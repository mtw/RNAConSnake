# Development Notes

## Local Environment

Use a local virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

This project packages the Python-side workflow logic, but still expects the external RNA analysis toolchain to be installed separately and available on `PATH`.

Two of those tools are on no package index and must be built or copied by
hand: **SISSIz** (<https://github.com/mtw/SISSIz>) and **`alifoldz.pl`** (RNAz
source tarball). Refolding additionally needs the ViennaRNA **Python module**
(`import RNA`) in the same environment. The
[container](../container/README.md) packages all of them and is the least
painful way to get a complete environment:

```bash
cd container && ./prepare-context.sh
docker build --platform linux/arm64 -t rnacs:0.3.0 .
```

Note that `--use-conda` does **not** work: most rules use Snakemake `run:`
directives, which execute in the Snakemake process and cannot take per-rule
conda environments. Use [`environment.yaml`](../environment.yaml) or the
container instead.

## Workflow Entry Points

There are two workflow entry points in the repository:

- [`snakefile`](../snakefile): workflow source of truth
- [`config.yaml`](../config.yaml): default workflow config source of truth

The packaged CLI `RNAcs` invokes workflow files copied from those root files at build time. In editable development installs, the CLI falls back to the root workflow files directly.

For repository-local use without shell-activating a virtual environment, the root launcher [`RNAcs`](../RNAcs) executes `python -m rnaconsnake.cli` with `PYTHONPATH=src` and automatically prefers `./.venv/bin/python` when it exists.

Important preprocessing knobs in [`config.yaml`](../config.yaml):

- `remove_gaponly_gapratio`
  Maximum allowed gap fraction per sequence.
- `remove_gaponly_max_n`
  Maximum allowed count of `N`/`n` characters per sequence.
  The default is `0`, which removes any sequence containing ambiguous `N` bases.
- `do_rscape`
  Enable optional per-hit `R-scape` analysis and summary integration.

## Common Commands

```bash
pytest -v
ruff check .
ruff format .
RNAcs --input-alignment /path/to/input.stk --output-dir /path/to/run_dir
RNAcs --check-deps
RNAcs --version
```

Lint and formatting are enforced in CI (`ruff check .`, `ruff format --check .`),
configured in `pyproject.toml` under `[tool.ruff]`.

Example CLI override for window sizes:

```bash
RNAcs --input-alignment /path/to/input.stk --maxbpspan 150 --maxbpspan 250 --cores all
```

Optional `R-scape` branch:

```bash
RNAcs --input-alignment /path/to/input.stk --rscape --cores all
RNAcs --check-deps --rscape
```

RNAz defaults to `-d -n`, so alignment shuffling is disabled unless explicitly requested. To allow shuffling for a run:

```bash
RNAcs --input-alignment /path/to/input.stk --rnaz-shuffle --cores all
```

To inspect which external tool commands RNAConSnake is configured to use, and where their executables resolve on the current machine, run:

```bash
RNAcs --show-tool-paths
```

Structured export can be triggered from the main CLI after a successful workflow run:

```bash
RNAcs \
  --input-alignment /path/to/input.stk \
  --output-dir /path/to/run_dir \
  --export-bundle /path/to/export_bundle \
  --cores all
```

## Tests And Fixtures

The test suite lives in [`tests/test_python_tools.py`](../tests/test_python_tools.py).

Fixtures are split into two groups:

- [`tests/fixtures`](../tests/fixtures): small synthetic fixtures for focused tool behavior
- [`tests/fixtures/real`](../tests/fixtures/real): small real-derived fixtures curated from repository input and workflow output data

The real-derived fixtures are intended to stay small enough for CI while still guarding against formatting and behavior regressions in the Stockholm helper tools.

The suite also includes an end-to-end Snakemake smoke test that stubs the external workflow toolchain with temporary fake executables on `PATH`. That keeps CI independent of the real binaries while still exercising the packaged workflow DAG.

To verify that two completed runs generated and processed the same deterministic candidate set from the same input alignment, run:

```bash
python -m rnaconsnake.tools.verify_run_consistency /path/to/run_a /path/to/run_b
```

This checks:
- `Lalifold/len_*/RC_*_0001.stk`
- split manifests and per-candidate split Stockholm files
- cleaned alignment manifests and per-candidate files under `generated_files/stk`

It deliberately does not compare stochastic downstream metrics such as `alifoldz` z-scores.

RNAConSnake should remain the scientific computation and export layer only. Any future HTML browser UI should remain outside the scope of this repository.

## Tool inventory

Everything under `src/rnaconsnake/tools/` is a standalone module with a `main()`,
runnable as `python -m rnaconsnake.tools.<name>`.

| module | role |
| --- | --- |
| `split_stockholm`, `remove_gaponly`, `strip_aln` | candidate preprocessing |
| `legacy_postprocess`, `alifold_maxcovar` | metric extraction, summary rendering |
| `stockholm_utils`, `alignment_io` | Stockholm/Clustal parsing, composition diagnostics |
| `null_model` | null-alignment generation (SISSIz, rnazRandomizeAln), pool pinning |
| `loci` | alignment-coordinate primitives |
| `dereplicate` | grouping overlapping windows into loci; representative selection |
| `calibration` | empirical FDR, q-values, filter funnel |
| `benchmark` | recovery scoring against a curated truth file, with null baseline |
| `benchmark_scaffold` | drafts a truth file from an `#=GC SS_cons` reference |
| `fold_region` | folds and plots an arbitrary alignment span |
| `refold` | `refold.pl` in Python: per-sequence constraints from the consensus, folded in-process |
| `alignment_report` | per-window screenability: where the screen is blind |
| `threshold_sweep` | cascade thresholds vs FDR, reusing a calibrated run |
| `sensitivity_envelope` | alignment subsets for measuring the detection floor |
| `versions` | toolchain provenance (`results/versions.yaml`) |
| `verify_run_consistency` | compares two completed runs |

Two are analysis aids rather than pipeline steps, and are not wired into the
DAG: `fold_region` and `sensitivity_envelope`. `threshold_sweep` likewise runs
against a finished calibration.

### refold

`refold` replaces ViennaRNA's `refold.pl` and the `RNAfold -C` it feeds:

```bash
python -m rnaconsnake.tools.refold \
  --alignment run_dir/.../RNAalifold_results.cleaned.aln \
  --consensus run_dir/.../<candidate>_dp.ps \
  --output refold.out
```

`--consensus` takes either an `RNAalifold -p` dot plot (pairs above
`--threshold`, default 0.9) or `RNAalifold`'s own output. `--constraints-only`
writes the stream `refold.pl` produced, for piping into `RNAfold -C`; the
default folds in-process through the ViennaRNA Python bindings and writes what
`RNAfold --noPS -C` would have written.

`rule run_refold_file` runs this instead of `refold.pl | RNAfold --noPS -C`,
which is why neither is a dependency any more. Output was checked
byte-for-byte against that pipe over 38 candidate windows from two flavivirus
alignments, in both constraint modes. Where the two differ, `refold.pl` fails
rather than disagreeing — see the module docstring. The `RNA` module is
imported only when folding; the constraint logic imports nothing.

ViennaRNA must be 2.7 or newer, and the module and the binaries must report
the same version:
`RNAcs --check-deps` compares them and refuses a mixed toolchain, since the
consensus structure comes from `RNAalifold` and the refold from the module.
Running `snakemake` directly bypasses that preflight — `results/versions.yaml`
still records both versions and whether they match.

## Shipping Model

Current release model:

- Python package via [`pyproject.toml`](../pyproject.toml)
- Snakemake orchestration
- root workflow files copied into the package at build time
- ViennaRNA and legacy post-processing commands as external runtime dependencies

Refolding already runs through the ViennaRNA Python bindings rather than
`refold.pl` and `RNAfold` (see [refold](#refold) above).

Planned future direction:

- loosen the hard external dependency on `RNALalifold`
- `alifoldz.pl` and `rnazRandomizeAln.pl` are the remaining Perl scripts;
  replacing either means reimplementing a statistic, not just text handling
- an integration tier that runs the container's real toolchain, since every
  workflow test here stubs it

## Null-model calibration arm

The calibration arm adds an `{arm}` dimension over the whole pipeline. The
invariant that makes it valid is that **both arms traverse the same rule DAG**,
so the code is structured to make divergence unrepresentable rather than merely
discouraged:

- Rules never hardcode an output prefix. They call `A(path)`, which prepends
  `arms/{arm}/`, or the empty string when the arm is disabled.
- `CandidatePaths` carries the same `arm_prefix`, so the Python body of a rule
  cannot drift from its declared outputs.
- Checkpoint lookups go through `checkpoint_wildcards()`, which adds the `arm`
  wildcard only when the arm is enabled.

If you add a rule downstream of `make_arm_alignment`, wrap its paths in `A()`
and take its `CandidatePaths` from `paths_for(wildcards)`. Do not add an
arm-specific variant of an existing rule; `tests/test_python_tools.py` asserts
that no rule name or rule body is duplicated.

Behaviour when the arm is disabled (`null.method: none`, the default) must stay
byte-identical to the pre-calibration pipeline. Verify with:

```bash
python -m rnaconsnake.tools.verify_run_consistency /path/to/run_a /path/to/run_b
```

Note that a plain `diff -r` between any two runs is expected to show
differences in `alifoldz` outputs (the tool shuffles internally) and in
PostScript/PDF/PNG artifacts (embedded `%%CreationDate`). Compare a run of the
changed code against *two* runs of the unchanged code before concluding that a
difference is yours.

The fast feedback loop for the arm is the CI toy target:

```bash
pytest -v -k null_arm_smoke
```

It runs the whole DAG, including two null replicates, against fake external
tools, and asserts that `funnel.tsv` has one row per filter stage per arm.
