# Development Notes

## Local Environment

Use a local virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

This project packages the Python-side workflow logic, but still expects `RNALalifold` to be installed separately and available on `PATH`.

## Workflow Entry Points

There are two workflow entry points in the repository:

- [`snakefile`](../snakefile): workflow source of truth
- [`config.yaml`](../config.yaml): default workflow config source of truth

The packaged CLI `rnaconsnake-run` invokes workflow files copied from those root files at build time. In editable development installs, the CLI falls back to the root workflow files directly.

## Common Commands

```bash
pytest -v
sh XFILE.sh
rnaconsnake-run --check-deps
rnaconsnake-run --version
```

If you prefer the convenience targets:

```bash
make venv
make install
make test
make check-deps
make run
```

## Tests And Fixtures

The test suite lives in [`tests/test_python_tools.py`](../tests/test_python_tools.py).

Fixtures are split into two groups:

- [`tests/fixtures`](../tests/fixtures): small synthetic fixtures for focused tool behavior
- [`tests/fixtures/real`](../tests/fixtures/real): small real-derived fixtures curated from repository input and workflow output data

The real-derived fixtures are intended to stay small enough for CI while still guarding against formatting and behavior regressions in the Stockholm helper tools.

The suite also includes a minimal end-to-end Snakemake smoke test that stubs `RNALalifold` with a temporary fake executable on `PATH`. That keeps CI independent of the real ViennaRNA binary while still exercising the packaged workflow DAG.

## Shipping Model

Current release model:

- Python package via [`pyproject.toml`](../pyproject.toml)
- Snakemake orchestration
- root workflow files copied into the package at build time
- `RNALalifold` as an external runtime dependency

Planned future direction:

- loosen the hard external dependency on `RNALalifold`
- evaluate wrapping ViennaRNA Python bindings directly from Python
