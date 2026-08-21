/*
 * commands.h — CLI command handlers for Clepsydra (M12).
 *
 * These wire the validated internal modules (TSV/Newick parsing, T1 likelihood,
 * the fixed-tree sampler, diagnostics) into runnable end-to-end commands over
 * the strict-TSV inputs of docs/input-formats.md. Errors are printed and a
 * non-zero status returned; the library itself never exits.
 */
#ifndef GD_CLI_COMMANDS_H
#define GD_CLI_COMMANDS_H

/* date-tree: fixed-topology T1 chronology from a Newick tree, a forms TSV, and
 * a calibrations TSV. Writes posterior_nodes.tsv, diagnostics.tsv and
 * run_manifest.tsv to the output directory. */
int cmd_date_tree(int argc, char **argv);

/* date-pair: P1 hierarchical pairwise divergence date from a forms TSV for two
 * languages (deterministic integration). Writes posterior_pair.tsv. */
int cmd_date_pair(int argc, char **argv);

/* calib-sensitivity: how much the T1 chronology depends on its calibrations.
 * --mode loco re-dates dropping each calibration in turn (does the model
 * recover the withheld age? how far does the rest of the tree move?);
 * --mode widen re-dates with every calibration window loosened by a factor.
 * Writes calibration_sensitivity.tsv to the output directory. */
int cmd_calib_sensitivity(int argc, char **argv);

#endif /* GD_CLI_COMMANDS_H */
