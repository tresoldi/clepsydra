/*
 * gd.h — public API for Clepsydra (libclepsydra).
 *
 * Include this one header to get the whole library. It pulls in the error
 * model, the shared types, the topology a result is reported against, and the
 * analysis pipelines that do the work:
 *
 *   gd_error.h     GD_E_* codes and gd_error
 *   gd_types.h     types shared across the API
 *   gd_tree.h      the phylogeny a chronology is reported against (read-only)
 *   gd_analysis.h  the pipelines: gd_analysis_tree_run, gd_analysis_pair_run
 *   gd_options.h   the run options, described once: name, type, bounds
 *   gd_report.h    the reported tables, and the renderers over them
 *
 * Strict C99; depends only on the C standard library and libm.
 *
 * Design constraints:
 *  - every pipeline works from in-memory buffers, never from paths, so the
 *    same engine runs natively, in a fuzzer, and in a browser with no
 *    filesystem;
 *  - the library never calls exit() or abort(); errors flow through gd_error,
 *    and every int-returning function returns a GD_E_* code (gd_error.h);
 *  - checked multiplication before every array allocation, because the
 *    element counts come from parsed input;
 *  - no allocation in the innermost likelihood loop;
 *  - iterative tree traversal, no deep recursion (the browser stack is small);
 *  - deterministic PCG32 RNG with independent streams, so a seeded run is
 *    bit-reproducible.
 *
 * Ownership vocabulary, used consistently across the API:
 *
 *   _create / _destroy     an owned handle; destroy exactly once
 *   _prepare / _run        a job is prepared, then run once
 *   _result_destroy        frees a result and everything it owns
 *   _read / _free          an enumerated array and its matching release
 *   const T *              borrowed; valid only as long as its owner is
 *
 * Anything returned as `const char *` inside a result — a node label, a model
 * name — is borrowed from that result and dies with it.
 *
 * Note on internal headers: src/gd_*.h are NOT part of this API. They are the
 * engine's own modules (the sampler, the likelihoods, the parsers), and the
 * test suite links against them deliberately in order to test units rather
 * than only the pipelines. Callers outside the repository should use this
 * header and nothing else.
 */
#ifndef GD_H
#define GD_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"
#include "gd_analysis.h"
#include "gd_options.h"
#include "gd_report.h"
#include "gd_treeset.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------- versioning */

/* The library version, e.g. "0.2.0". Borrowed static storage. */
const char *gd_version(void);

/* The git commit this library was built from, or "unknown" when it was built
 * outside a checkout. Borrowed static storage. Recorded in every run manifest,
 * so a result can always be traced back to the code that produced it. */
const char *gd_git_commit(void);

#ifdef __cplusplus
}
#endif
#endif /* GD_H */
