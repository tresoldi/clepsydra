/*
 * Runs the chronology sampler off the main thread.
 *
 * The sampler is a single synchronous call into WebAssembly, so while it runs
 * this worker cannot receive messages. Cancelling therefore means terminating
 * the worker from the page and starting a fresh one; re-instantiating the
 * module costs a few milliseconds once the bytes are cached. The engine's own
 * cancellation is reachable from here only through the progress callback, which
 * this worker uses to report — the page's Cancel button is the terminate path.
 */

/* global createClepsydra */

importScripts("clepsydra.js");

let engine = null;
let pending = null;      /* a request that arrived before the module resolved */

function post(type, payload) {
  self.postMessage(Object.assign({ type }, payload));
}

/*
 * Every argument goes on the heap, never through ccall's "string" type: that
 * copies the argument onto the wasm stack, and a forms file of any real size
 * overflows it — reported as an opaque out-of-bounds trap far from the cause.
 */
function callJson(name, args) {
  const ptrs = args.map((s) => engine.stringToNewUTF8(s === null ? "" : s));
  try {
    const out = engine.ccall(name, "number", ptrs.map(() => "number"), ptrs);
    const json = engine.UTF8ToString(out);
    engine.ccall("clepsydra_free", null, ["number"], [out]);
    return json;
  } finally {
    for (const p of ptrs) engine._free(p);
  }
}

createClepsydra({
  onProgress(chain, iter, total) {
    post("progress", { chain, iter, total });
    /* Returning nonzero would cancel; the page cancels by terminating this
     * worker instead, because a blocked worker cannot hear a message. */
    return 0;
  },
})
  .then((instance) => {
    engine = instance;
    post("ready", {
      version: instance.ccall("clepsydra_version", "string", [], []),
    });
    /* A visitor who presses Run immediately is not making a mistake; the
     * request waits for the engine instead of being refused. */
    if (pending !== null) {
      const req = pending;
      pending = null;
      runTree(req);
    }
  })
  .catch((error) => {
    post("fatal", {
      message: "the engine could not start: " + (error && error.message),
    });
  });

function runTree(req) {
  const started = Date.now();
  try {
    const json = callJson("clepsydra_date_tree_json",
                          [req.newick, req.forms, req.calibrations, req.opts]);
    post("result", { json, elapsedMs: Date.now() - started, tag: req.tag || "main" });
  } catch (error) {
    /* An exception here is the engine itself failing; a bad dataset comes back
     * as an ok:false payload instead. */
    post("fatal", { message: (error && error.message) || String(error) });
  }
}

self.onmessage = (event) => {
  const { type } = event.data;
  if (type !== "tree") {
    post("fatal", { message: "unknown request: " + type });
    return;
  }
  if (engine === null) { pending = event.data; return; }
  runTree(event.data);
};
