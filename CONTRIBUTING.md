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
rnaconsnake-run --input-alignment /path/to/input.stk --output-dir /path/to/run_dir
```

You can also point the packaged CLI at an arbitrary Stockholm input:

```bash
rnaconsnake-run --input-alignment /path/to/input.stk --cores all
```

To write outputs into a dedicated run directory:

```bash
rnaconsnake-run --input-alignment /path/to/input.stk --output-dir /path/to/run_dir --cores all
```

To override `maxbpspan` from the CLI:

```bash
rnaconsnake-run --input-alignment /path/to/input.stk --maxbpspan 150 --maxbpspan 250 --cores all
```

To enable optional `R-scape` analysis:

```bash
rnaconsnake-run --input-alignment /path/to/input.stk --rscape --cores all
rnaconsnake-run --check-deps --rscape
```

If you want to reuse existing outputs without rerunning rules just because the workflow code changed:

```bash
rnaconsnake-run --conservative --cores all
```

## Packaging

```bash
source .venv/bin/activate
python -m build
```

## Notes

- The root [`snakefile`](snakefile) and [`config.yaml`](config.yaml) are the workflow source of truth.
- The packaged workflow files are copied from those root files at build time.
- Keep tests green before opening a pull request.
