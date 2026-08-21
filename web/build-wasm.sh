#!/usr/bin/env bash
#
# Builds clepsydra for the browser. Requires the Emscripten SDK:
#
#   source ~/emsdk/emsdk_env.sh
#   ./web/build-wasm.sh
#
# The artifacts (web/clepsydra.js, web/clepsydra.wasm) are committed, so rerun
# this after any change under src/, include/ or web/clepsydra_wasm.c and commit
# the result. scripts/wasm-provenance.sh records what was built from, and the
# Pages deploy refuses a stale artifact rather than publishing an engine older
# than the sources it claims to be.
#
# Unlike the sibling projects there is no CMake step: clepsydra is a plain
# Makefile project, so emcc compiles the sources directly. It links with
# -sFILESYSTEM=0 — the engine reads only the buffers JavaScript hands it.
#
# STACK_SIZE is raised from emscripten's 64 KB default for headroom only. The
# datasets themselves must be passed as heap pointers, never as ccall "string"
# arguments: that type copies the argument onto the wasm stack, so a 100 KB
# forms file overflows it and fails as an opaque out-of-bounds trap. See the
# `call` helper in scripts/wasm_verify.js and web/worker.js.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"

if ! command -v emcc >/dev/null 2>&1; then
    echo "build-wasm: emcc not found; run 'source ~/emsdk/emsdk_env.sh' first" >&2
    exit 2
fi

version=$(sed -n 's/.*GD_VERSION_STRING\[\] = "\(.*\)".*/\1/p' \
          "$repo_dir/src/gd_version.c" | head -1)
: "${version:=0.0.0}"

# The same warning set the native build uses, and -Werror: the browser build is
# not a second-class citizen with looser standards.
emcc \
    -O3 \
    -std=c99 \
    -Wall -Wextra -Wpedantic -Wconversion \
    -Wshadow -Wstrict-prototypes -Wmissing-prototypes \
    -Werror \
    "$script_dir/clepsydra_wasm.c" \
    "$repo_dir"/src/*.c \
    -I"$repo_dir/include" \
    -I"$repo_dir/src" \
    -sFILESYSTEM=0 \
    -sMODULARIZE=1 \
    -sEXPORT_NAME=createClepsydra \
    -sALLOW_MEMORY_GROWTH=1 \
    -sENVIRONMENT=web,worker,node \
    -sSTACK_SIZE=262144 \
    -sEXPORTED_FUNCTIONS='["_clepsydra_date_tree_json","_clepsydra_date_pair_json","_clepsydra_defaults_json","_clepsydra_version","_clepsydra_git_commit","_clepsydra_free","_malloc","_free"]' \
    -sEXPORTED_RUNTIME_METHODS='["ccall","cwrap","UTF8ToString","stringToNewUTF8"]' \
    -o "$script_dir/clepsydra.js"

"$repo_dir/scripts/wasm-provenance.sh" --write

echo "built clepsydra $version:"
ls -lh "$script_dir/clepsydra.js" "$script_dir/clepsydra.wasm" | awk '{print "  " $9 "  " $5}'
echo
echo "check it against the CLI: ./scripts/wasm-verify.sh"
echo "preview:                  cd web && python3 -m http.server 8080"
