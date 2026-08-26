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
RNAcs --check-deps
```

If `RNALalifold` is available on `PATH`, you can also run the local workflow:

```bash
RNAcs --input-alignment /path/to/input.stk --output-dir /path/to/run_dir
```

You can also point the packaged CLI at an arbitrary Stockholm input:

```bash
RNAcs --input-alignment /path/to/input.stk --cores all
```

To write outputs into a dedicated run directory:

```bash
RNAcs --input-alignment /path/to/input.stk --output-dir /path/to/run_dir --cores all
```

To override `maxbpspan` from the CLI:

```bash
RNAcs --input-alignment /path/to/input.stk --maxbpspan 150 --maxbpspan 250 --cores all
```

To enable optional `R-scape` analysis:

```bash
RNAcs --input-alignment /path/to/input.stk --rscape --cores all
RNAcs --check-deps --rscape
```

If you want to reuse existing outputs without rerunning rules just because the workflow code changed:

```bash
RNAcs --conservative --cores all
```

To run the null-model calibration arm (see [`docs/usage.md`](docs/usage.md)):

```bash
RNAcs --input-alignment /path/to/input.stk --null-arm sissiz --null-replicates 10 --cores all
```

The full DAG including two null replicates runs on a toy set in a few minutes:

```bash
pytest -v -k null_arm_smoke
snakemake --profile profiles/test --config input_alignment=examples/all_DENVG_3UTR.relabel.stk
```

## Invariants a change must not break

These are enforced by tests; if one fails, the fix is almost never to relax the
test.

- **Both arms traverse the same rule DAG.** There are no arm-specific rules.
  Every rule downstream of alignment generation takes its paths from the `A()`
  helper, which prepends `arms/{arm}/`; when the arm is disabled the prefix is
  empty and the legacy layout is preserved byte-for-byte. Any divergence
  between arms invalidates the calibration.
- **No HTML reporting, templates, or web assets.** RNAConSnake is the
  scientific computation layer: it produces alignments, figures, metrics and
  structured export bundles. Rendering or browsing those results is out of
  scope, and the export bundle is the interface downstream consumers read.
- **Calibration counts on de-replicated loci, never on raw windows**, and using
  the same grouping the reported table uses.
- **Configuration is locked.** `docs/configuration_lock.md` records the
  representative rule and cascade thresholds, fixed on one clade before another
  was evaluated. Changing those values changes every published number.
- **Non-numeric sentinels stay non-numeric.** A failed AlifoldZ is `NA`, never
  `0.0`; a plausible-looking number is indistinguishable from a real score.

## Packaging

```bash
source .venv/bin/activate
python -m build
```

## Notes

- The root [`snakefile`](snakefile) and [`config.yaml`](config.yaml) are the workflow source of truth.
- The packaged workflow files are copied from those root files at build time.
- Keep tests green before opening a pull request.
