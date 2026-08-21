/*
 * gd_options.h — the run options, described once.
 *
 * gd_analysis_tree_config and gd_analysis_pair_config are filled from more
 * than one place: the CLI parses flags into them, the WebAssembly shim parses
 * `key=value` lines into them, and the browser page builds those lines. Before
 * D4 each surface carried its own spelling of every option, its own copy of
 * every default, and its own idea of which values were acceptable — so
 * `chains` defaulted to 4 in the library, 4 in the CLI and 2 on the page, and
 * the page could send a rate_logsd of 0 into the sampler that the CLI would
 * have refused.
 *
 * This is the one description. A row names an option, says where it lives in
 * its config struct and what type it has, and states its bounds. Every surface
 * reads the same rows, so a new option is one row and a renamed option cannot
 * be renamed in only three of four places.
 *
 * DEFAULTS ARE NOT HERE. They stay in gd_analysis_tree_config_default and
 * gd_analysis_pair_config_default, which are already the one correct source.
 * A surface that wants to show a default fills a config with those and reads
 * the value back through gd_option_format — which is what stops the numbers
 * being retyped, and also avoids narrowing a uint64_t seed through a double
 * on the way to a table.
 */
#ifndef GD_OPTIONS_H
#define GD_OPTIONS_H

#include "gd_types.h"
#include "gd_error.h"

#include <stdbool.h>
#include <stddef.h>

/* The types an option field can have. Every one is a distinct C type, and the
 * table asserts each row's size against it at compile time. */
typedef enum gd_opt_type {
    GD_OPT_SIZE = 0,   /* size_t   */
    GD_OPT_U64,        /* uint64_t */
    GD_OPT_DOUBLE,     /* double   */
    GD_OPT_BOOL        /* bool     */
} gd_opt_type;

/*
 * One option.
 *
 * `name` is the canonical spelling — what the shim takes in a `key=value`
 * line. `flag` is the CLI spelling, including the leading dashes; it is
 * usually `name` with underscores turned into dashes, but not always
 * (`progress_every` is `--progress`), so it is stated rather than derived.
 *
 * `offset` locates the field inside its config struct. It is paired with
 * `type` by the macro that builds the table, which asserts that the named
 * field really has that type's size. That assert is the point of this
 * struct: the shim's previous table reached its fields through a void * that
 * no compiler could check, and it was wrong once — two fields became bool
 * while their slots still said int, so the parser wrote four bytes into one
 * and quietly corrupted whatever sat after them. Nothing in the C build
 * caught it; a headless-browser smoke test did, by noticing the draws had
 * changed.
 *
 * `lo` and `hi` bound the accepted value. `lo_open` and `hi_open` make the
 * respective bound exclusive: "must be positive" is a lo_open 0, and a
 * contamination weight in [0, 1) is a hi_open 1. An unbounded end is
 * -HUGE_VAL or HUGE_VAL.
 */
typedef struct gd_option {
    const char *name;
    const char *flag;
    gd_opt_type type;
    size_t      offset;
    double      lo, hi;
    bool        lo_open, hi_open;
    /* What to print where a help listing would show the default value, when
     * the number itself would say nothing -- GD_ANALYSIS_KEEP_ALL renders as
     * SIZE_MAX. Cosmetic only: the value still comes from the config's
     * _default function, and gd_option_format stays round-trippable through
     * gd_option_set. NULL means "just print the value". */
    const char *default_note;
    const char *help;
} gd_option;

/* The two tables. Borrowed static storage; the arrays outlive any caller. */
const gd_option *gd_options_tree(size_t *n);
const gd_option *gd_options_pair(size_t *n);

/* Look an option up by canonical name or by CLI flag. NULL when absent.
 * `name_len` may be GD_OPT_MEASURE to have it measured. */
#define GD_OPT_MEASURE ((size_t)-1)
const gd_option *gd_option_find(const gd_option *opts, size_t n,
                                const char *name, size_t name_len);
const gd_option *gd_option_find_flag(const gd_option *opts, size_t n,
                                     const char *flag);

/*
 * Parse `value` and store it in `config` at this option's field, refusing
 * anything outside its bounds. Returns GD_E_OK, GD_E_BAD_NUMERIC for a value
 * that is not a finite number, or GD_E_INVALID_ARG for one out of range; the
 * error names the option and what was wrong with it.
 *
 * `value_len` may be GD_OPT_MEASURE.
 */
int gd_option_set(const gd_option *opt, void *config,
                  const char *value, size_t value_len, gd_error *err);

/*
 * Write this option's current value from `config` into `buf` as text, in the
 * form its own parser would accept. Used to show defaults without retyping
 * them. Returns GD_E_OK, or GD_E_RANGE if `buf` is too small.
 */
int gd_option_format(const gd_option *opt, const void *config,
                     char *buf, size_t bufsz, gd_error *err);

#endif /* GD_OPTIONS_H */
