/*
 * gd_version.c — build identity for Clepsydra.
 *
 * Split out of the retired gd_context.c: the version and the commit are the
 * only things in that module the library still needs, and they belong nowhere
 * near a dataset loader. Both are read into the run manifest and reported by
 * the browser page, so they are provenance, not decoration.
 *
 * The string is deliberately unchanged from the pre-promotion library. Every
 * committed example artifact carries it in `run_manifest.tsv`, and bumping it
 * here would rewrite those files for a reason that has nothing to do with the
 * science.
 */
#include "gd.h"

static const char GD_VERSION_STRING[] = "0.2.0";

const char *gd_version(void)
{
    return GD_VERSION_STRING;
}

const char *gd_git_commit(void)
{
    /* Populated at build time via -DGD_GIT_COMMIT="..."; "unknown" otherwise,
     * which is more honest for provenance than an empty string. */
#ifdef GD_GIT_COMMIT
    return GD_GIT_COMMIT;
#else
    return "unknown";
#endif
}
