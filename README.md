# RNAConSnake

RNAConSnake is a Snakemake-based RNA conservation workflow packaged with Python helper tools and a small CLI. The current release still expects `RNALalifold` as an external runtime dependency on `PATH`.

## Install

Recommended: use a local virtual environment.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install .[dev]
```

This keeps Python dependencies local to the project while leaving `RNALalifold` managed separately.

If you want a runtime-only install without test dependencies:

```bash
pip install .
```

Check that the external runtime dependencies are available:

```bash
rnaconsnake-run --check-deps
```

Show the installed package version:

```bash
rnaconsnake-run --version
```

## Run

From a working directory that contains your input files and desired outputs:

```bash
rnaconsnake-run --cores all
```

This launches Snakemake with the packaged workflow. The root [`snakefile`](snakefile) and [`config.yaml`](config.yaml) are the source of truth; build artifacts copy them into the installed package.

For repository-local development, the root workflow file is [`snakefile`](snakefile) and the convenience wrapper is [`XFILE.sh`](XFILE.sh).

## Package Contents

The package currently ships:

- the workflow runner CLI `rnaconsnake-run`
- helper CLIs:
  - `rnaconsnake-split-stockholm`
  - `rnaconsnake-remove-gaponly`
  - `rnaconsnake-strip-aln`
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
sh XFILE.sh
```

You can also use the included [`Makefile`](Makefile):

```bash
make venv
make install
make test
make check-deps
make run
```

The repository CI uses the same basic model: create a local virtual environment, install the package in editable mode, run the pytest suite, and smoke-test the CLI.

More development details are in [`docs/development.md`](docs/development.md).

## Tests

Run the test suite with:

```bash
pytest -v
```

Current coverage includes:

- CLI `--version`
- CLI `--check-deps`
- minimal end-to-end workflow smoke test with a stubbed `RNALalifold`
- Stockholm splitting on synthetic and real-derived fixtures
- `remove_gaponly` regression behavior
- `strip_aln --nosingle` regression behavior

The real-derived fixtures under [`tests/fixtures/real`](tests/fixtures/real) are intentionally small curated samples derived from `data/` and `Lalifold/` outputs.

## Optional Scripts

The repository also contains an optional post-processing helper:

- [`scripts/pp_RNALalifold.sh`](scripts/pp_RNALalifold.sh)

This is not part of the default packaged runtime interface. It has additional external dependencies and environment assumptions beyond the main workflow.

## Release Checklist

1. Create a fresh virtual environment and install the package.
2. Run `pytest -v`.
3. Run `rnaconsnake-run --version`.
4. Run `rnaconsnake-run --check-deps` in an environment with `RNALalifold` available.
5. Dry-run the workflow with representative inputs.
6. Confirm the workflow smoke test passes in `pytest -v`.
7. Update the version in [`pyproject.toml`](pyproject.toml) and [`src/rnaconsnake/__init__.py`](src/rnaconsnake/__init__.py).
8. Update [`CHANGELOG.md`](CHANGELOG.md).
