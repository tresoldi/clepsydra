# Clepsydra — hierarchical Bayesian lexical dating in portable C99.
# The browser engine is built separately by web/build-wasm.sh (needs emsdk);
# `make shim-check` compiles its adapter with the native compiler so the shim
# cannot rot unnoticed on a machine without Emscripten.

CC      ?= gcc

# The warning set is the same in every configuration; only the optimisation and
# the assertion policy change. It is not -Werror by default: a newer compiler
# should not break a user's build. CI turns warnings into errors instead, where
# the compiler version is pinned by the runner image.
WARNINGS  = -std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow \
            -Wstrict-prototypes -Wmissing-prototypes

# Three configurations (see docs/STYLE.md):
#
#   default   -O2 -g, assertions ON. What `make` builds, and what the examples
#             and the sanitizer tests run: an assertion that never runs is not
#             a check, so the everyday build keeps them.
#   debug     -O0 -g, assertions ON. For stepping through a sampler.
#   release   -O2 -DNDEBUG, assertions OFF. What a distributor would ship;
#             build it at least once before a release, because NDEBUG is the
#             one configuration where an assertion's *side effects* would go
#             missing — there are none, and this is how that stays true.
#
# Each writes to its own object directory, so switching does not silently mix
# objects built under different assertion policies.
CFLAGS  ?= $(WARNINGS) -O2 -g
SANFLAGS  ?= -fsanitize=address -fno-omit-frame-pointer
INCLUDE   = -Iinclude
LDLIBS    = -lm

# Build-time git commit stamped into the run manifest for provenance. Falls
# back to "unknown" outside a git checkout.
GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
CPPFLAGS += -DGD_GIT_COMMIT='"$(GIT_COMMIT)"'

SRC_DIR   = src
OBJ_DIR   = obj
SAN_OBJ   = obj_san
BUILD     = build
TEST_DIR  = tests/unit
LIK_DIR   = tests/likelihood
REG_DIR   = tests/regression
SIM_DIR   = tests/simulation

# ---- Release identity ----
#
# The single source of the version is src/gd_version.c; this is read from it so
# the soname and the pkg-config file cannot drift from what gd_version()
# reports. VERSION_MAJOR drives the soname: a 0.x series makes no ABI promise
# across minor bumps, so the soname carries the minor too.
VERSION       := $(shell sed -n 's/.*GD_VERSION_STRING\[\] = "\(.*\)".*/\1/p' $(SRC_DIR)/gd_version.c)
VERSION_ABI   := $(shell echo $(VERSION) | cut -d. -f1-2)

PREFIX      ?= /usr/local
DESTDIR     ?=
LIBDIR      ?= $(PREFIX)/lib
INCLUDEDIR  ?= $(PREFIX)/include
BINDIR      ?= $(PREFIX)/bin
PKGCONFIGDIR ?= $(LIBDIR)/pkgconfig

SHARED_SO     = libclepsydra.so
SHARED_SONAME = $(SHARED_SO).$(VERSION_ABI)
SHARED_REAL   = $(SHARED_SO).$(VERSION)

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS     = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))
PIC_OBJ  = obj_pic
PIC_OBJS = $(patsubst $(SRC_DIR)/%.c,$(PIC_OBJ)/%.o,$(SRCS))
SAN_OBJS = $(patsubst $(SRC_DIR)/%.c,$(SAN_OBJ)/%.o,$(SRCS))

TEST_SRCS = $(wildcard $(TEST_DIR)/test_*.c)
TEST_BINS = $(patsubst $(TEST_DIR)/test_%.c,$(BUILD)/test_%,$(TEST_SRCS))

LIK_SRCS = $(wildcard $(LIK_DIR)/test_*.c)
LIK_BINS = $(patsubst $(LIK_DIR)/test_%.c,$(BUILD)/lik_%,$(LIK_SRCS))

REG_SRCS = $(wildcard $(REG_DIR)/test_*.c)
REG_BINS = $(patsubst $(REG_DIR)/test_%.c,$(BUILD)/reg_%,$(REG_SRCS))

SIM_SRCS = $(wildcard $(SIM_DIR)/test_*.c)
SIM_BINS = $(patsubst $(SIM_DIR)/test_%.c,$(BUILD)/sim_%,$(SIM_SRCS))

.PHONY: all native test clean benchmark shim-check analyze fuzz shared install uninstall

all: native

native: $(BUILD)/libclepsydra.a $(BUILD)/clepsydra

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(SAN_OBJ):
	mkdir -p $(SAN_OBJ)

$(BUILD):
	mkdir -p $(BUILD)

# ---- Provenance stamp ----
#
# gd_version.c embeds GIT_COMMIT through -D, and make cannot see that value
# change: nothing on disk moved, so gd_version.o is never rebuilt and every
# run_manifest.tsv keeps reporting the commit of whenever it was last
# compiled. That is precisely the thing the manifest exists to get right, and
# it was silently wrong -- a stale stamp showed up as a phantom example diff
# on every commit of the D-series.
#
# The stamp file holds the current value. It is rewritten only when the value
# differs, so an unchanged HEAD does not force a rebuild.
.PHONY: force-stamp
force-stamp:

%/.git-commit: force-stamp
	@mkdir -p $(@D)
	@printf '%s' '$(GIT_COMMIT)' | cmp -s - $@ 2>/dev/null || \
	    printf '%s' '$(GIT_COMMIT)' > $@

# ---- Header dependencies ----
#
# Without these, editing a header rebuilds nothing. The result is not a stale
# binary but a MIXED one: the library compiled against the old struct and the
# test against the new, so a field written at one offset is read at another.
# That is how it surfaced -- a test asserting on a bool two fields after an
# inserted size_t, failing for reasons that had nothing to do with the change.
# A full rebuild "fixed" it, which is the tell.
#
# -MMD -MP writes a .d beside each .o listing the headers it used; -MP adds
# phony targets so a deleted header does not wedge the build.
DEPFLAGS = -MMD -MP

# ---- Native library ----
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(DEPFLAGS) $(INCLUDE) -c $< -o $@

$(OBJ_DIR)/gd_version.o: $(OBJ_DIR)/.git-commit

$(BUILD)/libclepsydra.a: $(OBJS) | $(BUILD)
	ar rcs $@ $(OBJS)

# ---- Shared library ----
#
# Separate objects because a shared library needs -fPIC and the static one does
# not; mixing them is how a static archive silently acquires relocations it
# cannot have. The soname carries major.minor: this is a 0.x series and the
# structs in include/ are young, so a minor bump is allowed to break ABI and
# says so rather than pretending otherwise.
$(PIC_OBJ):
	@mkdir -p $@

$(PIC_OBJ)/%.o: $(SRC_DIR)/%.c | $(PIC_OBJ)
	$(CC) $(CFLAGS) -fPIC $(CPPFLAGS) $(DEPFLAGS) $(INCLUDE) -c $< -o $@

$(PIC_OBJ)/gd_version.o: $(PIC_OBJ)/.git-commit

$(BUILD)/$(SHARED_REAL): $(PIC_OBJS) | $(BUILD)
	$(CC) -shared -Wl,-soname,$(SHARED_SONAME) $(PIC_OBJS) $(LDLIBS) -o $@
	@ln -sf $(SHARED_REAL) $(BUILD)/$(SHARED_SONAME)
	@ln -sf $(SHARED_SONAME) $(BUILD)/$(SHARED_SO)

shared: $(BUILD)/$(SHARED_REAL)

# No -Isrc: the CLI is an ordinary consumer of the public API and must stay
# buildable against include/ alone. If this rule ever needs src/ again, the
# thing to fix is the header that should have been public, not this line.
#
# Linked against the static archive by NAME, not by -lclepsydra: once a shared
# library exists in the same directory the linker prefers it, and the CLI
# quietly acquires a runtime dependency on a .so that is not on the loader path
# until ldconfig has run. An installed command that cannot start is a poor way
# to discover that.
$(BUILD)/clepsydra: cli/main.c cli/commands.c $(BUILD)/libclepsydra.a | $(BUILD)
	$(CC) $(CFLAGS) $(INCLUDE) cli/main.c cli/commands.c \
	    $(BUILD)/libclepsydra.a $(LDLIBS) -o $@

# ---- Sanitizer library (for tests) ----
$(SAN_OBJ)/%.o: $(SRC_DIR)/%.c | $(SAN_OBJ)
	$(CC) $(CFLAGS) $(SANFLAGS) $(CPPFLAGS) $(DEPFLAGS) $(INCLUDE) -c $< -o $@

$(SAN_OBJ)/gd_version.o: $(SAN_OBJ)/.git-commit

$(BUILD)/libclepsydra_san.a: $(SAN_OBJS) | $(BUILD)
	ar rcs $@ $(SAN_OBJS)

# ---- Tests (sanitizer-clean) ----
test: shim-check $(TEST_BINS) $(LIK_BINS) $(REG_BINS) $(SIM_BINS)
	@for t in $(TEST_BINS) $(LIK_BINS) $(REG_BINS) $(SIM_BINS); do \
		echo "=== $$t ==="; \
		$$t || exit 1; \
	done
	@echo "All tests passed."

# The WebAssembly adapter is only compiled by web/build-wasm.sh, which needs
# emsdk. Compiling it here with the native toolchain (its non-Emscripten branch)
# means a change that breaks the browser surface fails on any machine, not just
# on the one that happens to have Emscripten installed.
shim-check: | $(BUILD)
	@$(CC) $(CFLAGS) $(INCLUDE) -fsyntax-only web/clepsydra_wasm.c
	@echo "=== web/clepsydra_wasm.c: compiles clean ==="

# Tests are the one consumer allowed into src/. That is deliberate: a unit test
# exists to exercise a module directly, not only through the pipelines that
# happen to call it. Every other consumer builds against include/ alone.
$(BUILD)/test_%: $(TEST_DIR)/test_%.c $(BUILD)/libclepsydra_san.a | $(BUILD)
	$(CC) $(CFLAGS) $(SANFLAGS) $(INCLUDE) -Isrc $< -L$(BUILD) -lclepsydra_san $(LDLIBS) -o $@

$(BUILD)/lik_%: $(LIK_DIR)/test_%.c $(BUILD)/libclepsydra_san.a | $(BUILD)
	$(CC) $(CFLAGS) $(SANFLAGS) $(INCLUDE) -Isrc -I$(TEST_DIR) $< -L$(BUILD) -lclepsydra_san $(LDLIBS) -o $@

$(BUILD)/reg_%: $(REG_DIR)/test_%.c $(BUILD)/libclepsydra_san.a | $(BUILD)
	$(CC) $(CFLAGS) $(SANFLAGS) $(INCLUDE) -Isrc -I$(TEST_DIR) $< -L$(BUILD) -lclepsydra_san $(LDLIBS) -o $@

$(BUILD)/sim_%: $(SIM_DIR)/test_%.c $(BUILD)/libclepsydra_san.a | $(BUILD)
	$(CC) $(CFLAGS) $(SANFLAGS) $(INCLUDE) -Isrc -I$(TEST_DIR) $< -L$(BUILD) -lclepsydra_san $(LDLIBS) -o $@

# ---- Static analysis ----
# Two analyzers, both already part of the toolchains CI pins, so this adds no
# dependency. They are separate targets because they find different things:
# gcc's is stronger on allocation and lifetime, clang's on path-sensitive null
# and dead-value reasoning. Neither runs in `make test` — analysis is slower
# than compilation and belongs in CI and in a deliberate local run.
#
# `analyze` must stay at zero findings; see docs/static-analysis.md for the
# baseline, what was fixed to reach it, and the policy on suppressions.
ANALYZE_SRCS = $(SRCS) cli/main.c cli/commands.c web/clepsydra_wasm.c bench/bench_mcmc.c

# Distributions that default to a hardened toolchain (Ubuntu, and so CI) turn
# on _FORTIFY_SOURCE, which enables -Wformat-truncation among others. Without
# this, `make analyze` reports clean on a distribution that does not, and the
# same tree fails in CI -- which is how a 192-byte buffer holding a 256-byte
# error message got as far as a pull request.
ANALYZE_HARDENING = -Wall -Wextra -D_FORTIFY_SOURCE=2

.PHONY: analyze analyze-gcc analyze-clang

analyze: analyze-gcc analyze-clang
	@echo "=== static analysis: clean ==="

# Neither analyzer exits non-zero on a finding, so the targets check the output
# rather than the status — a target that reports "clean" while the tool printed
# a warning is worse than no target at all.
analyze-gcc:
	@echo "=== gcc -fanalyzer ==="
	@out=$$(for f in $(ANALYZE_SRCS); do \
		gcc -std=c99 $(INCLUDE) -Isrc -fanalyzer -O2 $(ANALYZE_HARDENING) \
		    -c $$f -o /dev/null 2>&1; \
	done); \
	if echo "$$out" | grep -q 'warning:'; then \
		echo "$$out"; \
		echo "gcc -fanalyzer: $$(echo "$$out" | grep -c 'warning:') finding(s)"; \
		exit 1; \
	fi

analyze-clang:
	@echo "=== clang --analyze ==="
	@out=$$(for f in $(ANALYZE_SRCS); do \
		clang -std=c99 $(INCLUDE) -Isrc --analyze \
		    -Xanalyzer -analyzer-output=text $$f -o /dev/null 2>&1; \
	done); \
	if echo "$$out" | grep -q 'warning:'; then \
		echo "$$out"; \
		echo "clang --analyze: $$(echo "$$out" | grep -c 'warning:') finding(s)"; \
		exit 1; \
	fi

# ---- Fuzzing ----
# The parsers consume untrusted text: the browser page lets a visitor paste a
# Newick tree and upload a forms file, so a stranger's bytes reach this C. Each
# harness calls the production parser directly — a harness that pre-validates
# its input only proves the validator works.
#
#   make fuzz                 build the harnesses
#   make fuzz-run             a short bounded run over the seed corpora (CI)
#   build/fuzz_newick corpus/newick   an open-ended run, locally
#
# clang only: libFuzzer is a clang feature. ASan and UBSan are on, because a
# fuzzer without sanitizers finds only the crashes that were going to happen
# anyway. Findings that stick go in tests/regression/ as ordinary tests.
FUZZ_DIR    = fuzz
FUZZ_SRCS   = $(wildcard $(FUZZ_DIR)/fuzz_*.c)
FUZZ_BINS   = $(patsubst $(FUZZ_DIR)/fuzz_%.c,$(BUILD)/fuzz_%,$(FUZZ_SRCS))
FUZZ_CC    ?= clang
# -fno-sanitize-recover is not optional here. Without it UBSan *prints* a
# runtime error and the fuzzer carries on and exits 0, so a UB finding would
# scroll past a green run — the flags matched the CI UBSan job's, which does
# set it. Unsigned wrap is excluded for the same reason CI excludes it: the
# PCG32 stream relies on it deliberately.
FUZZ_FLAGS ?= -fsanitize=fuzzer,address,undefined,integer \
              -fno-sanitize=unsigned-integer-overflow,unsigned-shift-base \
              -fno-sanitize-recover=all \
              -fno-omit-frame-pointer
# Seconds per harness in `fuzz-run`. Short by design: this is a regression
# gate, not a campaign. Raise it for a real hunt.
FUZZ_SECONDS ?= 30

.PHONY: fuzz fuzz-run

fuzz: $(FUZZ_BINS)

$(BUILD)/fuzz_%: $(FUZZ_DIR)/fuzz_%.c $(SRCS) | $(BUILD)
	$(FUZZ_CC) -std=c99 $(FUZZ_FLAGS) -O1 -g $(INCLUDE) -Isrc \
	    $< $(SRCS) $(LDLIBS) -o $@

# New units go to a scratch directory under build/ (gitignored), not into the
# committed corpus: libFuzzer writes discoveries to the FIRST corpus argument,
# and a run that quietly adds a few hundred files to the repository is not a
# check, it is a mess. The committed corpus is passed second, read-only, so a
# run still starts from real parser states.
#
# To grow the committed corpus deliberately, run a harness directly and
# re-minimise — see docs/fuzzing.md.
fuzz-run: fuzz
	@for b in $(FUZZ_BINS); do \
		name=$$(basename $$b); \
		corpus=$(FUZZ_DIR)/corpus/$${name#fuzz_}; \
		scratch=$(BUILD)/fuzz-scratch/$${name#fuzz_}; \
		mkdir -p $$scratch; \
		echo "=== $$name ($(FUZZ_SECONDS)s) ==="; \
		$$b $$scratch $$corpus -max_total_time=$(FUZZ_SECONDS) \
		    -print_final_stats=1 -artifact_prefix=$(BUILD)/ || exit 1; \
	done
	@echo "All fuzz targets survived."

# ---- Alternative configurations ----
# Recursive rather than a variable switch, so `make debug` and `make release`
# cannot pick up objects the default build left behind.
.PHONY: debug release

debug:
	@$(MAKE) --no-print-directory native \
	    CFLAGS="$(WARNINGS) -O0 -g" OBJ_DIR=obj_debug BUILD=build_debug
	@echo "=== debug build in build_debug/ (assertions on, -O0) ==="

release:
	@$(MAKE) --no-print-directory native \
	    CFLAGS="$(WARNINGS) -O2 -DNDEBUG" OBJ_DIR=obj_release BUILD=build_release
	@echo "=== release build in build_release/ (assertions off, -O2) ==="

# ---- Benchmark ----
#
# Deliberately NOT part of `make test`: it takes minutes, and a slow suite is a
# suite people stop running. The numbers are CPU seconds, not wall clock --
# see bench/bench_mcmc.c on why.
#
# Override the workload to match a question:
#   make benchmark BENCH_WARMUP=4000 BENCH_SAMPLES=4000 BENCH_CHAINS=4
BENCH_WARMUP  ?= 1500
BENCH_SAMPLES ?= 1500
BENCH_CHAINS  ?= 2
BENCH_SEED    ?= 42
BENCH_REPEATS ?= 3

$(BUILD)/bench_mcmc: bench/bench_mcmc.c $(BUILD)/libclepsydra.a | $(BUILD)
	$(CC) $(CFLAGS) $(INCLUDE) bench/bench_mcmc.c \
	    -L$(BUILD) -lclepsydra $(LDLIBS) -o $@

benchmark: $(BUILD)/bench_mcmc
	@echo "clepsydra benchmark — CPU seconds, seeded, $(BENCH_REPEATS) repeats"
	@$(CC) --version | head -1 | sed 's/^/  compiler: /'
	@echo "  workload: $(BENCH_WARMUP) warmup + $(BENCH_SAMPLES) samples x $(BENCH_CHAINS) chains, seed $(BENCH_SEED)"
	@echo
	@printf 'dataset\trepeat\twarmup\tsamples\tchains\tprepare_cpu_s\trun_cpu_s\titers\titers_per_s\n'
	@$(BUILD)/bench_mcmc romance examples/romance_tree \
	    $(BENCH_WARMUP) $(BENCH_SAMPLES) $(BENCH_CHAINS) $(BENCH_SEED) $(BENCH_REPEATS)
	@$(BUILD)/bench_mcmc five_ie examples/five_ie_reproduction \
	    $(BENCH_WARMUP) $(BENCH_SAMPLES) $(BENCH_CHAINS) $(BENCH_SEED) $(BENCH_REPEATS)
	@echo
	@echo "Compare against docs/benchmark-baseline.md. An effect smaller than the"
	@echo "spread above is not distinguishable from scheduling noise." 

-include $(PIC_OBJS:.o=.d)
-include $(OBJS:.o=.d)
-include $(SAN_OBJS:.o=.d)

# ---- Install ----
#
# DESTDIR is honoured so a packager can stage into a build root. Headers go in
# their own directory: <clepsydra/gd.h> rather than <gd.h>, so a name as
# generic as gd_types.h cannot collide with whatever else is on the include
# path.
install: native shared $(BUILD)/clepsydra.pc
	@install -d $(DESTDIR)$(LIBDIR) $(DESTDIR)$(INCLUDEDIR)/clepsydra \
	            $(DESTDIR)$(BINDIR) $(DESTDIR)$(PKGCONFIGDIR)
	install -m 644 $(BUILD)/libclepsydra.a $(DESTDIR)$(LIBDIR)/
	install -m 755 $(BUILD)/$(SHARED_REAL) $(DESTDIR)$(LIBDIR)/
	ln -sf $(SHARED_REAL) $(DESTDIR)$(LIBDIR)/$(SHARED_SONAME)
	ln -sf $(SHARED_SONAME) $(DESTDIR)$(LIBDIR)/$(SHARED_SO)
	install -m 644 include/*.h $(DESTDIR)$(INCLUDEDIR)/clepsydra/
	install -m 755 $(BUILD)/clepsydra $(DESTDIR)$(BINDIR)/
	install -m 644 $(BUILD)/clepsydra.pc $(DESTDIR)$(PKGCONFIGDIR)/
	@echo "installed clepsydra $(VERSION) under $(DESTDIR)$(PREFIX)"

uninstall:
	rm -f $(DESTDIR)$(LIBDIR)/libclepsydra.a
	rm -f $(DESTDIR)$(LIBDIR)/$(SHARED_REAL) \
	      $(DESTDIR)$(LIBDIR)/$(SHARED_SONAME) \
	      $(DESTDIR)$(LIBDIR)/$(SHARED_SO)
	rm -rf $(DESTDIR)$(INCLUDEDIR)/clepsydra
	rm -f $(DESTDIR)$(BINDIR)/clepsydra
	rm -f $(DESTDIR)$(PKGCONFIGDIR)/clepsydra.pc

$(BUILD)/clepsydra.pc: | $(BUILD)
	@printf 'prefix=%s\n' '$(PREFIX)'                        > $@
	@printf 'libdir=%s\n' '$(LIBDIR)'                       >> $@
	@printf 'includedir=%s\n\n' '$(INCLUDEDIR)'             >> $@
	@printf 'Name: clepsydra\n'                             >> $@
	@printf 'Description: Hierarchical Bayesian lexical dating in portable C99\n' >> $@
	@printf 'URL: https://github.com/tresoldi/clepsydra\n'   >> $@
	@printf 'Version: %s\n' '$(VERSION)'                    >> $@
	@printf 'Libs: -L$${libdir} -lclepsydra\n'               >> $@
	@printf 'Libs.private: -lm\n'                            >> $@
	@printf 'Cflags: -I$${includedir}\n'                      >> $@

clean:
	rm -rf $(OBJ_DIR) $(SAN_OBJ) $(PIC_OBJ) $(BUILD) obj_debug build_debug obj_release build_release
