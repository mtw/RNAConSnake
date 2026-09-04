# Contributing

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
```

## Development loop

```bash
pytest -v
ruff check . && ruff format .     # CI enforces both; 110 char line length
RNAcs --check-deps
```

Run workflow (requires RNALalifold on PATH):
```bash
RNAcs --input-alignment aln.stk --output-dir run_dir --cores all
RNAcs --input-alignment aln.stk --maxbpspan 150 --maxbpspan 250 --cores all
RNAcs --input-alignment aln.stk --rscape --cores all
RNAcs --input-alignment aln.stk --null-arm sissiz --null-replicates 10 --cores all
```

Full DAG with two null replicates on toy set:
```bash
pytest -v -k null_arm_smoke
```

## Invariants (test-enforced)

- **Both arms traverse the same rule DAG** (no arm-specific rules; `A()` helper prepends `arms/{arm}/` path prefix only).
- **No HTML, templates, or web assets.** Export bundle is the downstream interface.
- **Calibration counts on de-replicated loci** using the same grouping as the reported table.
- **Configuration is locked** (representative rule and cascade thresholds in `config.yaml` are pinned; changing them changes all published numbers).
- **Non-numeric sentinels stay non-numeric** (failed AlifoldZ is `NA`, not `0.0`).

## Packaging

```bash
python -m build
```

Root `snakefile` and `config.yaml` are source of truth; packaged copies
generated at build time. Keep tests green before opening a PR.
