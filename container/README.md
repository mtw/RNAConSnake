# Container

Reproducible image pinning SISSIz 0.2.0, alifoldz.pl, and Debian tools not in conda-forge. Multi-arch native builds (x86_64 or ARM64); SISSIz compiled in-image. On clusters with Apptainer, pull CI's tested x86_64 image instead of building.

## Build

Requires: SISSIz source at tag `0.2.0` (`https://github.com/mtw/SISSIz`), `alifoldz.pl` from RNAz source tarball.

```bash
cd container
./prepare-context.sh        # gathers SISSIz + alifoldz.pl
docker build -t rnacs:0.3.0 .
```

Or override locations:
```bash
SISSIZ_SRC=~/src/SISSIz ALIFOLDZ=~/.local/share/RNAz/perl/alifoldz.pl ./prepare-context.sh
```

Native build is default (what you want—emulated builds compile SISSIz under QEMU). Cross-build with `docker buildx build --platform linux/amd64` (slow).

Sources are **not committed**; image rebuilt from your copies, so nothing third-party redistributed.

## Docker

```bash
docker run --rm \
  -v "$PWD/examples:/data:ro" -v "$PWD/runs/myrun:/work" \
  rnacs:0.3.0 \
  --input-alignment /data/all_JEVG_3UTR.relabel.stk \
  --output-dir /work --maxbpspan 100 --cores 10 --null-arm sissiz --null-replicates 100
```

Entrypoint is `RNAcs`; arguments passed through. Runs as root (outputs owned by root on host).

## Apptainer

CI builds x86_64 image, tests, and pushes to ghcr.io. Pull that:

```bash
apptainer pull rnacs_0.3.0.sif docker://ghcr.io/mtw/rnaconsnake:0.3.0
apptainer run --cleanenv rnacs_0.3.0.sif --check-deps
```

Version tags publish version; dev tags (`dev-<sha>`) exist for pre-release tries:

```bash
apptainer pull rnacs_dev.sif docker://ghcr.io/mtw/rnaconsnake:dev-1c94c92
```

Private ghcr packages: `pull` fails with `unauthorized` until visibility set to public. To keep private, authenticate with `read:packages` token:

```bash
export APPTAINER_DOCKER_USERNAME=<user> APPTAINER_DOCKER_PASSWORD=<token>
```

Run with `--cleanenv` (prevents host environment shadowing toolchain):

```bash
apptainer run --cleanenv \
  --bind "$PWD/examples:/data:ro" --bind "$PWD/runs/ebov:/work" \
  rnacs_0.3.0.sif \
  --input-alignment /data/aln.stk --output-dir /work --cores 10
```

Key differences from Docker:
- **`--cleanenv` required** (module loads, conda prefixes shadow image)
- **Run as yourself** (not root; outputs owned by you)
- **Image FS read-only** (writes go to bound `/work`)
- **Working dir is host's** (pass `--output-dir` explicitly)

Create run directory before binding. Keep `--cores` ≤ job allocation on schedulers.

## Reproducibility across machines

- **SISSIz unseeded:** splitting null replicates across hosts gives different alignments (statistically fine, but then only reproducible with `--null-pool`).
- **Check versions:** Record `results/versions.yaml` from each host and verify they agree before pooling arms.
- **Safest split:** one input alignment per machine, not arms across machines.
