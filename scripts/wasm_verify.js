/*
 * wasm_verify.js — runs an example through the WebAssembly engine under node,
 * so scripts/wasm-verify.sh can compare it with the same run on the CLI.
 *
 * M10's acceptance gate asks for identical seeded native-vs-Wasm results within
 * a documented floating-point tolerance. This half produces the Wasm numbers;
 * the shell script produces the native ones and does the comparing.
 *
 * Usage:
 *   node scripts/wasm_verify.js tree <dir> <warmup> <samples> <chains> <seed> \
 *        [rate_logmean] [rate_logsd]
 *   node scripts/wasm_verify.js pair <forms.tsv> <langA> <langB> <t_lo> <t_hi>
 */
"use strict";

const fs = require("fs");
const path = require("path");

const repo = path.resolve(__dirname, "..");
const createClepsydra = require(path.join(repo, "web", "clepsydra.js"));

/*
 * Calls an entry point with every argument passed as a HEAP pointer.
 *
 * ccall's "string" type copies the argument onto the *wasm stack*, which is
 * 256 KB: a 100 KB forms file plus a tree and calibrations overflows it, and
 * the failure is an opaque "memory access out of bounds" far from the cause.
 * Datasets are user-sized, so they go on the heap and are freed here.
 */
function call(mod, name, args) {
  const ptrs = args.map((s) => mod.stringToNewUTF8(s === null ? "" : s));
  try {
    const out = mod.ccall(name, "number", ptrs.map(() => "number"), ptrs);
    const json = mod.UTF8ToString(out);
    mod.ccall("clepsydra_free", null, ["number"], [out]);
    return JSON.parse(json);
  } finally {
    for (const p of ptrs) mod._free(p);
  }
}

async function main() {
  const [mode] = process.argv.slice(2);
  let progressCalls = 0;

  const mod = await createClepsydra({
    onProgress() { progressCalls++; return 0; },
  });

  let out;
  if (mode === "tree") {
    const [dir, warmup, samples, chains, seed, logmean, logsd] =
      process.argv.slice(3);
    const opts = [
      `warmup=${warmup}`,
      `samples=${samples}`,
      `chains=${chains}`,
      `seed=${seed}`,
      "relaxed=0",
      logmean !== undefined ? `rate_logmean=${logmean}` : null,
      logsd !== undefined ? `rate_logsd=${logsd}` : null,
    ].filter(Boolean).join("\n");

    out = call(mod, "clepsydra_date_tree_json",
               [fs.readFileSync(path.join(dir, "tree.nwk"), "utf8"),
                fs.readFileSync(path.join(dir, "data", "forms.tsv"), "utf8"),
                fs.readFileSync(path.join(dir, "calibrations.tsv"), "utf8"),
                opts]);
    out.progress_calls = progressCalls;
  } else if (mode === "pair") {
    const [forms, a, b, tLo, tHi] = process.argv.slice(3);
    out = call(mod, "clepsydra_date_pair_json",
               [fs.readFileSync(forms, "utf8"), a, b, "",
                `t_lo=${tLo}\nt_hi=${tHi}`]);
  } else {
    process.stderr.write("usage: wasm_verify.js tree|pair ...\n");
    process.exit(2);
  }

  if (out.ok !== true) {
    process.stderr.write("wasm run failed: " + JSON.stringify(out) + "\n");
    process.exit(1);
  }
  process.stdout.write(JSON.stringify(out) + "\n");
}

main().catch((e) => {
  process.stderr.write("wasm_verify: " + (e && e.stack ? e.stack : e) + "\n");
  process.exit(1);
});
