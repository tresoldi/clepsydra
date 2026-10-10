/*
 * main.c — CLI entry point for Clepsydra.
 *
 * Commands: version, date-tree (T1 fixed-topology chronology), date-pair (P1
 * pairwise divergence), calib-sensitivity (calibration-sensitivity workflows).
 * See cli/commands.c and docs/cli-manual.md.
 */
#include "gd.h"
#include "commands.h"

#include <stdio.h>
#include <string.h>

static void usage(const char *prog)
{
    fprintf(stderr, "clepsydra %s\n", gd_version());
    fprintf(stderr, "Usage: %s <command> [options]\n", prog);
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "  date-tree          fixed-topology T1 chronology (MCMC)\n");
    fprintf(stderr, "  date-pair          P1 pairwise divergence date\n");
    fprintf(stderr, "  calib-sensitivity  leave-one-calibration-out / widen sensitivity\n");
    fprintf(stderr, "  version\n");
    fprintf(stderr, "Run a command with no options for its usage.\n");
}

int main(int argc, char **argv)
{
    if (argc < 2) { usage(argv[0]); return 1; }
    if (strcmp(argv[1], "version") == 0) {
        printf("clepsydra %s\n", gd_version());
        return 0;
    }
    if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        usage(argv[0]); return 0;
    }
    if (strcmp(argv[1], "date-tree") == 0) return cmd_date_tree(argc, argv);
    if (strcmp(argv[1], "date-trees") == 0) return cmd_date_trees(argc, argv);
    if (strcmp(argv[1], "date-pair") == 0) return cmd_date_pair(argc, argv);
    if (strcmp(argv[1], "calib-sensitivity") == 0) return cmd_calib_sensitivity(argc, argv);
    fprintf(stderr, "Unknown command: %s\n", argv[1]);
    usage(argv[0]);
    return 1;
}
