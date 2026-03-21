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

## Workflow Entry Points

There are two workflow entry points in the repository:

- [`snakefile`](../snakefile): workflow source of truth
- [`config.yaml`](../config.yaml): default workflow config source of truth

The packaged CLI `rnaconsnake-run` invokes workflow files copied from those root files at build time. In editable development installs, the CLI falls back to the root workflow files directly.

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
sh XFILE.sh /path/to/input.stk --output-dir /path/to/run_dir
rnaconsnake-run --check-deps
rnaconsnake-run --version
```

Example CLI override for window sizes:

```bash
rnaconsnake-run --input-alignment /path/to/input.stk --maxbpspan 150 --maxbpspan 250 --cores all
```

Optional `R-scape` branch:

```bash
rnaconsnake-run --input-alignment /path/to/input.stk --rscape --cores all
rnaconsnake-run --check-deps --rscape
```

If you prefer the convenience targets:

```bash
make venv
make install
make test
make check-deps
make run INPUT=/path/to/input.stk OUTPUT=/path/to/run_dir
```

## Tests And Fixtures

The test suite lives in [`tests/test_python_tools.py`](../tests/test_python_tools.py).

Fixtures are split into two groups:

- [`tests/fixtures`](../tests/fixtures): small synthetic fixtures for focused tool behavior
- [`tests/fixtures/real`](../tests/fixtures/real): small real-derived fixtures curated from repository input and workflow output data

The real-derived fixtures are intended to stay small enough for CI while still guarding against formatting and behavior regressions in the Stockholm helper tools.

The suite also includes an end-to-end Snakemake smoke test that stubs the external workflow toolchain with temporary fake executables on `PATH`. That keeps CI independent of the real binaries while still exercising the packaged workflow DAG.

## Shipping Model

Current release model:

- Python package via [`pyproject.toml`](../pyproject.toml)
- Snakemake orchestration
- root workflow files copied into the package at build time
- ViennaRNA and legacy post-processing commands as external runtime dependencies

Planned future direction:

- loosen the hard external dependency on `RNALalifold`
- evaluate wrapping ViennaRNA Python bindings directly from Python
