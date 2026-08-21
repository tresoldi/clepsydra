#!/usr/bin/env bash
#
# Records which sources produced the committed WebAssembly build, and checks
# that they still match.
#
#   ./scripts/wasm-provenance.sh --write    # after web/build-wasm.sh
#   ./scripts/wasm-provenance.sh --check    # asserts web/clepsydra.wasm is current
#
# web/clepsydra.wasm is committed rather than built on deploy, which means the
# page can ship an engine older than src/ without anything going red. The stamp
# in web/BUILD_INFO closes that: it hashes every source that goes into the
# artifact, so a source edit without a rebuild is a failing check rather than a
# silently stale deploy.
#
# Hashes sources rather than the artifact: emcc is not reproducible byte for
# byte across toolchain versions, so comparing a rebuilt binary would fail for
# reasons that have nothing to do with staleness. It also keeps this runnable
# anywhere bash and git are — no emsdk, no compiler — which is what lets CI gate
# the deploy on it.
#
# cli/ is deliberately NOT hashed: it is not part of the browser build, and
# including it would fail the check for changes that cannot reach the page.
set -euo pipefail

repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
stamp="$repo_dir/web/BUILD_INFO"

# Tracked *and* untracked-but-not-ignored, so a new source file counts before it
# is committed; in a clean CI checkout the two sets are the same.
sources_hash() {
    cd "$repo_dir"
    git ls-files --cached --others --exclude-standard -z -- \
            src include web/clepsydra_wasm.c \
        | LC_ALL=C sort -z \
        | while IFS= read -r -d '' path; do
              # The name goes into the hash too, so a rename counts as a change.
              [ -f "$path" ] && printf '%s  %s\n' \
                  "$(sha256sum < "$path" | cut -d' ' -f1)" "$path"
          done \
        | sha256sum | cut -d' ' -f1
}

mode="${1:---check}"
current="$(sources_hash)"

case "$mode" in
--write)
    commit="$(git -C "$repo_dir" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    printf 'sources_sha256 %s\n' "$current" > "$stamp"
    printf 'git_commit %s\n' "$commit" >> "$stamp"
    echo "wasm-provenance: recorded ${current:0:12} (${commit})"
    ;;
--check)
    if [ ! -f "$stamp" ]; then
        echo "wasm-provenance: web/BUILD_INFO is missing; run web/build-wasm.sh" >&2
        exit 1
    fi
    recorded="$(awk '$1 == "sources_sha256" { print $2 }' "$stamp")"
    if [ "$recorded" != "$current" ]; then
        cat >&2 <<EOF
wasm-provenance: web/clepsydra.wasm was built from different sources.

  recorded  $recorded
  current   $current

Rebuild and commit the artifact:

  source ~/emsdk/emsdk_env.sh
  ./web/build-wasm.sh
EOF
        exit 1
    fi
    echo "wasm-provenance: web/clepsydra.wasm is current (${current:0:12})"
    ;;
*)
    echo "usage: wasm-provenance.sh [--write|--check]" >&2
    exit 2
    ;;
esac
