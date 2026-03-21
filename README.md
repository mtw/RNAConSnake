# RNAConSnake

RNAConSnake is a Snakemake-based RNA conservation workflow packaged with Python helper tools and a small CLI. The current release still expects the ViennaRNA and legacy post-processing command-line tools to be available on `PATH`.

## Install

Recommended: use a local virtual environment.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install .[dev]
```

This keeps Python dependencies local to the project while leaving the external RNA analysis tools managed separately.

If you want a runtime-only install without test dependencies:

```bash
pip install .
```

Check that the external runtime dependencies are available:

```bash
rnaconsnake-run --check-deps
```

If you want to include optional `R-scape` analysis in a run, check that dependency set with:

```bash
rnaconsnake-run --check-deps --rscape
```

The full workflow currently expects these external commands on `PATH`:

- `RNALalifold`
- `esl-reformat`
- `RNAz`
- `alifoldz.pl`
- `RNAalifold`
- `ps2eps`
- `epstopdf`
- `refold.pl`
- `RNAfold`
- `magick`

Show the installed package version:

```bash
rnaconsnake-run --version
```

## Run

From a working directory that contains your input files and desired outputs:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --cores all
```

To write all workflow outputs into a separate working directory:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --output-dir /path/to/run_dir --cores all
```

To override the default `maxbpspan` values from [`config.yaml`](config.yaml), repeat `--maxbpspan`:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --maxbpspan 150 --maxbpspan 250 --cores all
```

To enable optional per-hit `R-scape` analysis:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --rscape --cores all
```

If you want Snakemake to trust existing outputs more conservatively and avoid reruns caused only by workflow code changes, use:

```bash
rnaconsnake-run --conservative --cores all
```

This launches Snakemake with the packaged workflow. The root [`snakefile`](snakefile) and [`config.yaml`](config.yaml) are the source of truth; build artifacts copy them into the installed package.

For repository-local development, the root workflow file is [`snakefile`](snakefile) and the convenience wrapper is [`XFILE.sh`](XFILE.sh).

## Configuration

The default workflow configuration lives in [`config.yaml`](config.yaml).

Two useful preprocessing settings are:

```yaml
remove_gaponly_gapratio: 0.5
remove_gaponly_max_n: 0
do_rscape: false
```

They control the `remove_gaponly` filtering step as follows:

- `remove_gaponly_gapratio`
  Maximum allowed gap fraction per sequence in the split Stockholm alignment.
- `remove_gaponly_max_n`
  Maximum allowed count of `N` or `n` characters per sequence.
  `0` means sequences containing any `N`/`n` are removed.

Example: allow up to two ambiguous bases per sequence:

```yaml
remove_gaponly_max_n: 2
```

## Package Contents

The package currently ships:

- the workflow runner CLI `rnaconsnake-run`
- helper CLIs:
  - `rnaconsnake-split-stockholm`
  - `rnaconsnake-remove-gaponly`
  - `rnaconsnake-strip-aln`
  - `rnaconsnake-alifold-maxcovar`
- Python helper module for legacy post-processing metric extraction
- the root workflow files copied into the installed package at build time
- pytest-based regression tests in the repository

## Development

Install in editable mode:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

Typical local development loop:

```bash
source .venv/bin/activate
pip install -e .[dev]
pytest -v
sh XFILE.sh /path/to/input.stk --output-dir /path/to/run_dir
```

You can also use the included [`Makefile`](Makefile):

```bash
make venv
make install
make test
make check-deps
make run INPUT=/path/to/input.stk OUTPUT=/path/to/run_dir
make clean-data
```

The repository CI uses the same basic model: create a local virtual environment, install the package in editable mode, run the pytest suite, and smoke-test the CLI.

More development details are in [`docs/development.md`](docs/development.md).
For a publication-style description of the workflow and outputs, see [`docs/pipeline_summary.md`](docs/pipeline_summary.md).
For a ready-to-use manuscript companion-data bundle template, see [`publication_companion`](publication_companion).
To assemble that bundle automatically from a multi-gene result tree, use [`scripts/build_publication_companion.py`](scripts/build_publication_companion.py).

## Tests

Run the test suite with:

```bash
pytest -v
```

Current coverage includes:

- CLI `--version`
- CLI `--check-deps`
- minimal end-to-end workflow smoke test with a stubbed `RNALalifold`
- workflow smoke test covering the full legacy post-processing chain with fake external tools
- Stockholm splitting on synthetic and real-derived fixtures
- `remove_gaponly` regression behavior
- `strip_aln --nosingle` regression behavior

The real-derived fixtures under [`tests/fixtures/real`](tests/fixtures/real) are intentionally small curated samples derived from `data/` and `Lalifold/` outputs.

## Optional Scripts

The repository also contains an optional post-processing helper:

- [`scripts/pp_RNALalifold.sh`](scripts/pp_RNALalifold.sh)

This is not part of the default packaged runtime interface. It documents the original shell-based implementation that the current Snakemake workflow now aims to match.

## Cleaning Outputs

To remove generated workflow outputs, run:

```bash
snakemake --cores 1 clean
```

or use:

```bash
make clean-data
```

## Release Checklist

1. Create a fresh virtual environment and install the package.
2. Run `pytest -v`.
3. Run `rnaconsnake-run --version`.
4. Run `rnaconsnake-run --check-deps` in an environment with the required external tools available.
5. Dry-run the workflow with representative inputs.
6. Confirm the workflow smoke test passes in `pytest -v`.
7. Update the version in [`pyproject.toml`](pyproject.toml) and [`src/rnaconsnake/__init__.py`](src/rnaconsnake/__init__.py).
8. Update [`CHANGELOG.md`](CHANGELOG.md).
