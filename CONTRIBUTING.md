# Contributing

## Local Setup

Use a local virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

## Common Development Loop

```bash
source .venv/bin/activate
pytest -v
rnaconsnake-run --check-deps
```

If `RNALalifold` is available on `PATH`, you can also run the local workflow:

```bash
sh XFILE.sh
```

## Make Targets

```bash
make venv
make install
make test
make check-deps
make dist
make run
```

## Notes

- The root [`snakefile`](snakefile) and [`config.yaml`](config.yaml) are the workflow source of truth.
- The packaged workflow files are copied from those root files at build time.
- Keep tests green before opening a pull request.
