# Changelog

All notable changes to this project should be documented in this file.

The format is based on Keep a Changelog, and this project aims to use Semantic Versioning.

## [Unreleased]

### Added
- **Export bundle schema 1.2.0.** The per-sequence constrained refold
  (`*_refold.out`) is exported — it was computed for every candidate and never
  carried into a bundle. `input_alignment_type` now reports the format the run
  was actually given.
- **`RNAcs --export-overwrite`**, and `threshold_sweep --representative`.
- `split_stockholm` is configurable through `tools:` like every other tool; it
  was the one hardcoded invocation in the workflow.

### Changed
- **`{candidate}.refold.json` is now `generated_files/consensus/{candidate}.consensus.json`**,
  and the `extract-refold` subcommand is `extract-consensus`. The file holds the
  RNAalifold consensus structure, not the refold; it was named for the leg that
  runs beside it. The export artifact is `consensus_structure` /
  "Consensus secondary structure" rather than "Refold summary".
- `results/versions.yaml` is written with `yaml.safe_dump`.
- PyYAML is a declared dependency; it was imported directly but only ever
  arrived transitively through Snakemake.
- The `reciprocal_overlap` column of the benchmark recovery table moved from
  first to beside `overlap_fraction`, so the table leads with the element.
- `remove_gaponly --infmt` defaults to `stockholm`, the only value it accepts.
- The test suite is split by area under `tests/`, with shared fakes and builders
  in `tests/helpers.py`.

### Fixed
- **`RNAcs --export-bundle` deleted the directory it was pointed at.** The
  wrapper forced `overwrite=True` on every run, so an existing directory was
  removed without asking. Overwriting is now opt-in (`--export-overwrite`) and
  refuses any directory holding files that are not part of an export bundle.
- **`threshold_sweep` ignored the representative rule**, always clustering with
  `widest`, so a sweep of a run configured otherwise described different scores.
- **`sensitivity_envelope` crashed after writing its outputs** when an alignment
  yielded no subsets — `min()` over an empty sequence, raised at the last line.
- **`--config null=<non-mapping>` raised an AttributeError traceback** instead of
  saying what was wrong.
- `run_with_progress` built progress bars in its `finally` block only to close
  them, and could wait on a child whose pipe nobody was reading.
- **`consensus_mfe` recorded RNAz's mean single-sequence MFE, not the consensus
  MFE.** Both labels sat in one regex alternation, and `re.search` returns the
  leftmost match — RNAz prints `Mean single sequence MFE` first, so on every
  real output the wrong quantity was recorded (on an RNAz 2.1.1 example, -20.93
  where the consensus MFE was -18.74). It reached `RNAConSnake.log.csv`,
  `RNAConSnake.md`, every `*.summary.json`, and the `consensus_mfe` column of
  the exported `candidates.csv`. **Re-export any bundle whose `consensus_mfe`
  values matter; the other columns are unaffected.**
- **`verify_run_consistency` reported two calibrated runs as identical without
  comparing anything.** With the null-model arm enabled every output lives under
  `arms/real/`, where the tool never looked; it found no window lengths, compared
  nothing, and exited 0. It now resolves the analysis root the way the export
  bundle does, and refuses to call an empty comparison agreement.
- **`alifold_maxcovar` exited non-zero on success**, returning the covariation
  count as the process exit status — so any caller checking the status read a
  candidate carrying covariation as a failed tool. The count was always on
  stdout; the status is now 0.
- **`#=GR` annotation was not merged across interleaved Stockholm blocks.** The
  sequences and `#=GC` lines beside it were, so a round trip emitted several
  short `#=GR` rows against a full-length alignment.
- **`split_stockholm --accession` decided nothing.** Both branches returned the
  accession, so an accession always won and the default never reached `#=GF ID`.
  The default now prefers `#=GF ID` and the flag selects the accession, as its
  help always claimed. RNALalifold writes no `#=GF AC`, so pipeline candidate
  names are unchanged.

## [0.3.0] - 2026-08-28

### Added
- **Null-model calibration arm** (`--null-arm`). Optionally reruns the *same*
  rule DAG over simulated alignments carrying no structural signal, and turns
  the score distributions into empirical FDR / q-values in
  `results/calibration/`. SISSIz (default) or `rnazRandomizeAln.pl`. Counting is
  on de-replicated loci, and two-stage mode records the FDR as conditional.
  See [`docs/usage.md`](docs/usage.md).
- **De-replication of overlapping windows** (`--dereplicate`). RNALalifold
  reports one element many times over; `RNAConSnake.nr.csv` now gives one row
  per locus alongside the untouched per-window table, and `RNAConSnake.md`
  leads with it. Methods: `containment` (default), `substructure`, `overlap`,
  `none`. See [`docs/dereplication.md`](docs/dereplication.md).
- **Positive-control benchmark** (`--benchmark`). Scores recovery of curated
  elements from `resources/benchmark/` against the calibrated q-values, with a
  null-arm baseline (`--null-loci`) — without which a recovery count is
  uninterpretable. `benchmark_scaffold` drafts truth files from an
  `#=GC SS_cons` reference, filling in coordinates but never element labels.
- **Alignment screenability report** (`generated_files/alignment_screenability.tsv`).
  Per-window occupancy and identity, distinguishing regions where nothing *can*
  be reported from regions where power is merely reduced.
- **Export bundle schema 1.1.0.** Carries locus and q-value columns plus
  `dereplication`/`calibration` manifest blocks. All new fields are optional, so
  1.0.0 consumers keep working.
- **Analysis aids**, outside the DAG: `fold_region` (fold and plot an arbitrary
  span), `threshold_sweep` (thresholds vs FDR), `sensitivity_envelope`
  (alignment subsets for measuring the detection floor).
- **`tools/refold`**, a reimplementation of ViennaRNA's `refold.pl` that needs
  only the ViennaRNA Python bindings: no Perl, and no `RNAfold` subprocess,
  since the constrained fold runs in-process. Verified byte-identical to
  `refold.pl | RNAfold --noPS -C` over 38 windows from two flavivirus
  alignments, in both constraint modes.
- **The AlifoldZ zero-variance guard ships as a patch**
  (`container/patches/alifoldz-zero-variance.patch`), applied by
  `prepare-context.sh`. It was previously a local modification of
  `alifoldz.pl`, so the image could not be rebuilt from published sources
  alone. Without it, a candidate whose shuffles all fold to the same energy
  divides by zero and dies; with it, that candidate is recorded as unscored.
- **Reproducibility artefacts:** `results/versions.yaml`, `CITATION.cff`,
  `environment.yaml`, a `profiles/test/` CI profile, and a
  [container](container/README.md) pinning the whole toolchain — including
  SISSIz (<https://github.com/mtw/SISSIz>) and `alifoldz.pl`, neither of which
  is on any package index.
- Input format detection for the initial `RNALalifold` run: `.stk` uses `-f S`,
  `.aln` uses `-f C`.

### Fixed
- **Lower-case alignments hung the workflow indefinitely.** `refold.pl` matches
  only upper-case Clustal rows, so a MAFFT-style alignment parsed as nothing and
  its `<>` fell through to STDIN, parking a scheduler slot forever. Case is now
  normalised once in `strip_aln`, and every workflow subprocess pins `stdin`.
- **Lower-case alignments also produced garbage AlifoldZ z-scores.**
  `alifoldz.pl` read `0 sequences of 0 columns`, exited 0, and printed `9999` —
  its running-minimum initialiser — which was recorded as a z-score.
  `extract-alifoldz` now refuses non-results, and a crashed AlifoldZ is `NA`
  rather than a plausible-looking `0.0`.
- **A calibrated run was not reproducible, and said it was.** SISSIz seeds from
  the clock, so `null.seed` never reproduced a pool, yet the metadata claimed
  determinism. Corrected; `--null-pool` pins a previously generated pool.
- **AlifoldZ made q-values irreproducible** — 100 unseeded internal shuffles.
  `alifoldz_seed` (on by default) pins them per candidate; two runs with a
  pinned pool are now byte-identical.
- **A stray `alirna.ps` was written to the run root**, once per shuffle from
  every candidate in parallel. AlifoldZ now runs in a private scratch directory.
- **`fold_region` dropped a stray `alirna.ps` into the caller's directory.** It
  ran AlifoldZ with no working directory of its own; AlifoldZ shells out to
  RNAalifold without `--noPS`. It now runs in the output directory, like the
  workflow rule does.
- `parse_stockholm_records` silently truncated interleaved Stockholm alignments
  to their last block.
- The export bundle now resolves `arms/real/` when the null-model arm is enabled.
- **The dependency preflight checked the wrong commands.** It resolved hard-coded
  program names, so a tool redirected through `tools:` in the config file was
  never checked at the path the run would actually use, and branches switched on
  in the config file (`do_rscape`, `null.method`) contributed no dependencies
  unless they were also passed as flags.
- **`--benchmark` produced a recovery count with no null baseline.** The rule now
  passes every null arm's locus tables to the benchmark, so the report states how
  many elements the null arms "recover" too. `--benchmark` without the null arm
  now fails immediately with an explanation, rather than on a missing rule for
  the calibrated table it scores.
- **AlifoldZ was ranked in opposite directions in different places.**
  De-replication treats a more negative z-score as the stronger one; the summary
  reports and the export bundle sorted the raw score, ranking `-2` above `-3`.
- **`results/calibration/summary.json` omitted clustering parameters that change
  the q-values** (`max_container_width`, `container_min_coverage`,
  `representative_rule` under `thresholds`), and the export manifest copied that
  incomplete block downstream.
- **Only the default curated truth file was packaged.** `benchmark_truth` can
  name any file in `resources/benchmark/`, but the build hook copied one of
  them, so an installed run could not resolve the others.
- **`--benchmark` replaced the default targets instead of adding to them.**
  Naming a target makes Snakemake ignore `rule all`, so `RNAcs --benchmark`
  built the recovery table and whatever it depended on, silently skipping the
  rest of the run -- `results/versions.yaml` among it.

### Changed
- **`refold.pl` and the `RNAfold` binary are no longer dependencies.** The
  refold leg runs `rnaconsnake.tools.refold`, so neither has to be found on
  `PATH`, vendored into the container, or copied out of the ViennaRNA source
  tarball. It needs the ViennaRNA **Python module** instead, which must be the
  **same build as the ViennaRNA binaries**: `RNAcs --check-deps` requires the
  module, compares its version against every ViennaRNA binary, and refuses a
  mixed toolchain -- the consensus comes from the binaries and the refold from
  the module, so two builds would mean two sets of energy parameters in one
  run. The supported minimum is ViennaRNA 2.7. `results/versions.yaml` records
  the module's version beside the binaries' and states whether they match.
- **The command-line tool is now `RNAcs`.** `rnaconsnake-run` remains a
  deprecated alias.
- **The minimum Python is now 3.11.** `>=3.10` was never installable — snakemake
  8 requires 3.11, so `pip install` failed outright there. A test compares the
  declared floor against what the dependencies support.
- **Default behaviour changed in three places**, each altering results relative
  to earlier runs: `alifoldz_seed` is now set; `dereplicate.max_container_width`
  defaults to 120; `strip_aln` upper-cases sequences. The cascade thresholds and
  representative rule are locked in `config.yaml`, and a test pins them.
- **Ruff is now the linter and formatter**, enforced in CI (`ruff check`,
  `ruff format --check`) and configured in `pyproject.toml` at 110 columns,
  with pycodestyle, pyflakes, import order, pyupgrade, bugbear,
  comprehensions, simplifications, ruff's own checks, pytest style (tests
  only), naming, performance, return shape, modern idioms, tidy imports and a
  complexity ceiling set to today's worst function. The tree was formatted in
  one pass; the lint fixes it required were unused imports, import order,
  `zip(..., strict=)` on alignment rows, and a file handle left unclosed in the
  tests.
- **Tests and CI.** 196 tests at 92% line coverage (measured with subprocess
  tracing; the naive figure understates it, because most code runs through
  spawned tools). CI lints, runs the suite on 3.11-3.13, and verifies that the
  built sdist and wheel carry the workflow and every curated truth file. The
  workflow tests drive the real Snakemake DAG against stubbed external tools,
  including both arms of the calibration and the positive control. A separate
  CI job builds the container and screens the JEV 3'UTR example through the
  *real* toolchain -- the only tier where a change in what RNALalifold, RNAz or
  AlifoldZ prints would surface -- and fails if the run comes back without
  genuine scores. It is skipped on pull requests, where it would cost minutes
  rather than seconds.
- Internal cleanups: `SUMMARY_FIELDS` defined once and shared; `run_checked`
  uses `contextlib.ExitStack`; `normalize_rnaalifold_side_output` raises rather
  than silently no-oping; manifest rules share a helper; `refold_firstseq` and
  `refold_firststruc` removed from summaries and export bundles.

## [0.2.0] - 2026-03-24

### Added
- Installable Python package via `pyproject.toml`
- AGPL-3.0-or-later `LICENSE` file
- Packaged workflow runner CLI: `rnaconsnake-run`
- Python replacements for the former Perl helper tools
- Python helper module for legacy post-processing metric extraction and summary generation
- Python replacement for `alifold_maxcovar.pl`
- Helper tool CLIs: `rnaconsnake-split-stockholm`, `rnaconsnake-remove-gaponly`, and `rnaconsnake-strip-aln`
- Pytest-based regression tests for helper tools and CLI behavior
- End-to-end workflow smoke test with a stubbed `RNALalifold`
- Curated real-derived fixtures for Stockholm regression coverage
- GitHub Actions CI
- Developer docs for local virtualenv usage and workflow entry points

### Changed
- Workflow helper tools now run via packaged Python modules
- Snakemake workflow now implements the full legacy post-processing chain
- Repository workflow file is now tracked as lowercase `snakefile`
- Root `snakefile` and `config.yaml` are now the workflow source of truth for packaging
- Optional shell-based post-processing helper retained outside the default packaged runtime interface
- ViennaRNA and legacy post-processing commands remain external runtime dependencies

## [0.1.0] - 2026-03-19

### Added
- Initial packaged release of RNAConSnake
- Snakemake workflow distribution with Python helper tools
- Editable-install developer workflow
- Release checklist and local developer command guidance
