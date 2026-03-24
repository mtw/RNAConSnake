# Usage Guide

## Installation

Recommended setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

If you want a runtime-only install without test dependencies:

```bash
pip install .
```

## Repository-local launcher

From a repository checkout, you can use:

```bash
./rnaconsnake-run --version
```

It automatically prefers `./.venv/bin/python` when present and otherwise falls back to `python3`.

## Dependency checks

Check the standard runtime dependencies:

```bash
rnaconsnake-run --check-deps
```

Check the dependency set including optional `R-scape` support:

```bash
rnaconsnake-run --check-deps --rscape
```

To inspect which external tool commands RNAConSnake is configured to use, and where their executables resolve on the current machine:

```bash
rnaconsnake-run --show-tool-paths
```

## Basic run

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --cores all
```

To write all outputs into a dedicated working directory:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --output-dir /path/to/run_dir --cores all
```

## Common options

Override the default `maxbpspan` values:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --maxbpspan 150 --maxbpspan 250 --cores all
```

Enable optional per-hit `R-scape` analysis:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --rscape --cores all
```

RNAz defaults to `-d -n`, which disables alignment shuffling. To explicitly allow shuffling:

```bash
rnaconsnake-run --input-alignment /path/to/my_alignment.stk --rnaz-shuffle --cores all
```

If you want Snakemake to reuse existing outputs more conservatively and avoid reruns caused only by workflow code changes:

```bash
rnaconsnake-run --conservative --cores all
```

## Structured export

RNAConSnake can generate a structured export bundle after a successful run:

```bash
rnaconsnake-run \
  --input-alignment /path/to/my_alignment.stk \
  --output-dir /path/to/run_dir \
  --export-bundle /path/to/export_bundle \
  --cores all
```

You can also override export metadata, for example:

```bash
rnaconsnake-run \
  --input-alignment /path/to/my_alignment.stk \
  --output-dir /path/to/run_dir \
  --export-bundle /path/to/export_bundle \
  --export-feature-type other \
  --export-feature-label "Example aligned feature" \
  --export-dataset-label "Example dataset" \
  --cores all
```

## Configuration

The default workflow configuration lives in [`../config.yaml`](../config.yaml).

Useful preprocessing settings include:

```yaml
remove_gaponly_gapratio: 0.5
remove_gaponly_max_n: 0
do_rscape: false
rnaz_no_shuffle: true
```

Meaning:

- `remove_gaponly_gapratio`
  Maximum allowed gap fraction per sequence in the split Stockholm alignment.
- `remove_gaponly_max_n`
  Maximum allowed count of `N` or `n` characters per sequence.
  `0` means sequences containing any `N`/`n` are removed.
- `rnaz_no_shuffle`
  Controls whether RNAz is called with `-n`.
  The default is `true`, which disables alignment shuffling.

Example:

```yaml
remove_gaponly_max_n: 2
```

## Cleaning outputs

To remove generated workflow outputs from a run directory:

```bash
snakemake --cores 1 clean
```

or:

```bash
rm -rf Lalifold generated_files
```
