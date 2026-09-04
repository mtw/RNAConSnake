# Development

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
```

## Dependencies

**On PATH:** ViennaRNA 2.7+, RNAz, Infernal, R-scape (optional)

**Not on package index:** SISSIz 0.2.0 (build or copy by hand), alifoldz.pl (from RNAz source tarball)

For a complete environment, use the container (see [`container/README.md`](../container/README.md)).

## Workflow entry points

- [`snakefile`](../snakefile): workflow DAG
- [`config.yaml`](../config.yaml): default config
- Root launcher [`RNAcs`](../RNAcs): uses editable install or `./.venv/bin/python`

## Common commands

```bash
pytest -v
ruff check .
ruff format .
RNAcs --input-alignment aln.stk --output-dir run_dir --cores all
RNAcs --check-deps
RNAcs --show-tool-paths
```

## Tests

Split by area in [`tests/`](../tests). Fixtures in `tests/fixtures/` (synthetic) and `tests/fixtures/real/` (curated from real data). Smoke test stubs the external toolchain; CI includes real-toolchain tier.

Verify run consistency across two completed runs:

```bash
python -m rnaconsnake.tools.verify_run_consistency run_a run_b
```

## Tool inventory

Modules under `src/rnaconsnake/tools/`, each runnable as `python -m rnaconsnake.tools.<name>`:

- `split_stockholm`, `remove_gaponly`, `strip_aln`: preprocessing
- `legacy_postprocess`, `alifold_maxcovar`: metrics
- `stockholm_utils`, `alignment_io`: parsing
- `null_model`: null-alignment generation, pool pinning
- `dereplicate`, `calibration`, `benchmark`: de-replication, FDR, recovery
- `fold_region`, `threshold_sweep`, `sensitivity_envelope`: analysis aids
- `refold`: Python replacement for `refold.pl` (ViennaRNA bindings only, no Perl)
- `verify_run_consistency`: compare two runs

## Null-model arm

Rules never hardcode output prefixes; all use `A()` which prepends `arms/{arm}/` when the arm is enabled. Checkpoint lookups via `checkpoint_wildcards()`. No rule duplicates.

Verify byte-identical behavior when disabled:

```bash
python -m rnaconsnake.tools.verify_run_consistency run_a run_b
```

Fast feedback loop:

```bash
pytest -v -k null_arm_smoke
```
