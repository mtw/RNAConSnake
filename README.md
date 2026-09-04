[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19210864.svg)](https://doi.org/10.5281/zenodo.19210864)

# RNAConSnake

Workflow for discovering and prioritizing conserved RNA secondary structure in multiple sequence alignments. Combines RNALalifold discovery with filtering (RNAz, AlifoldZ, covariation), de-replication, and optional empirical q-values via null-model calibration.

## Quick start

```bash
# Install
python3 -m venv .venv
source .venv/bin/activate
pip install -e .

# Run (accepts .stk or .aln)
RNAcs --input-alignment aln.stk --output-dir run_dir --cores all
RNAcs --input-alignment aln.aln --output-dir run_dir --cores all

# With calibration (100 null replicates)
RNAcs --input-alignment aln.stk --output-dir run_dir \
  --null-arm sissiz --null-replicates 100 --cores all
```

Check dependencies:
```bash
RNAcs --check-deps
```

## What you get

- `RNAConSnake.nr.csv` — one row per locus (de-replicated via de-replicating overlapping windows)
- `RNAConSnake.log.csv` — all RNALalifold windows
- Structure/alignment figures and metrics per candidate
- `alignment_screenability.tsv` — blind-region reporting
- With `--null-arm`: empirical q-values, FDR, and null score distributions
- With `--benchmark`: recovery of known elements (Positive control) against q-values

## Dependencies

External toolchain on `PATH`: ViennaRNA 2.7+, RNAz, Infernal, R-scape (optional).

Two tools are not on package indices and must be built by hand:
- **SISSIz 0.2.0** (<https://github.com/mtw/SISSIz>) — for null-model arm
- **alifoldz.pl** — from RNAz source tarball

Note: Snakemake's `--use-conda` flag is not supported in this workflow (tools must be pre-installed on PATH). See `environment.yaml` for conda packages, or use the container (see [`container/README.md`](container/README.md)).

## Documentation

- **Usage:** [`docs/usage.md`](docs/usage.md)
- **De-replication:** [`docs/dereplication.md`](docs/dereplication.md)
- **Pipeline overview:** [`docs/pipeline_summary.md`](docs/pipeline_summary.md)
- **Development:** [`docs/development.md`](docs/development.md)
- **Container:** [`container/README.md`](container/README.md)

## Container

Pre-built x86_64 image published to ghcr.io on every push to `dev`:

```bash
apptainer pull rnacs.sif docker://ghcr.io/mtw/rnaconsnake:dev-<sha>
apptainer run --cleanenv --bind data:/data --bind run:/work rnacs.sif \
  --input-alignment /data/aln.stk --output-dir /work --cores all
```

See [`container/README.md`](container/README.md) for build and Apptainer details.

## License

AGPL-3.0-or-later
