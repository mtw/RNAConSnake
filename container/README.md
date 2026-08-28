# RNAcs container

Native **linux/arm64** image for running the calibrated screen on Apple Silicon.
Built for distributing long calibration runs across machines.

## Why a container

Most of the toolchain installs from conda-forge/bioconda, but three pieces do
not, and getting them right by hand on each machine is the failure mode this
image removes:

| piece | why it needs special handling |
| --- | --- |
| **SISSIz** | not in bioconda on any platform; compiled from source ([https://github.com/mtw/SISSIz](https://github.com/mtw/SISSIz)) |
| **alifoldz.pl** | ships in the RNAz *source* tarball, not the conda package |
| **ps2eps**, **epstopdf** | Debian packages; not in conda-forge |

Refolding needs no vendored tool: it runs through the ViennaRNA Python module,
which must be importable in the image's environment.

R-scape is also absent from bioconda and is **not** included: the default
configuration has `do_rscape: false`. Enabling R-scape means extending the
image.

## Build

```bash
cd container
./prepare-context.sh        # gathers SISSIz source + the two perl helpers
docker build --platform linux/arm64 -t rnacs:0.3.0 .
```

`prepare-context.sh` collects three things that are not on any package index
and must already be present on the build machine:

| piece | where to get it |
| --- | --- |
| SISSIz source | clone [https://github.com/mtw/SISSIz](https://github.com/mtw/SISSIz); the build tree, not a binary |
| `alifoldz.pl` | the **RNAz source tarball** (`perl/` directory) — absent from the conda package |

`ps2eps` and `epstopdf` come from Debian inside the image and need no action.

Override the locations if yours differ:

```bash
git clone https://github.com/mtw/SISSIz ~/src/SISSIz

SISSIZ_SRC=~/src/SISSIz \
ALIFOLDZ=~/.local/share/RNAz/perl/alifoldz.pl \
./prepare-context.sh
```

The vendored sources are **not committed** — the image is rebuilt from your own
copies, so nothing third-party is redistributed through this repository.

## Run

```bash
docker run --rm --platform linux/arm64 \
  -v "$PWD/examples:/data:ro" \
  -v "$PWD/runs/myrun:/work" \
  rnacs:0.3.0 \
  --input-alignment /data/all_JEVG_3UTR.relabel.stk \
  --output-dir /work \
  --maxbpspan 100 --cores 10 --no-progress \
  --null-arm sissiz --null-replicates 100
```

Mount inputs read-only at `/data` and the run directory at `/work`. The
entrypoint is `RNAcs`, so arguments are passed straight through.

Run it as root, which is the image's default. Passing `--user` with a uid that
has no entry in the image's `/etc/passwd` fails: snakemake looks the user up
(`getpwuid(): uid not found`). Outputs in the mounted run directory are then
owned by root on the host.

## Distributing to other machines

`docker save` / `docker load` moves the built image without rebuilding:

```bash
docker save rnacs:0.3.0 | gzip -1 | ssh <host> 'gunzip | docker load'
ssh <host> 'docker run --rm --platform linux/arm64 rnacs:0.3.0 --check-deps'
```

## Reproducibility caveats when splitting work across machines

The image pins the toolchain, so results are comparable across hosts in a way
that separately-installed toolchains would not be. Two things still need care:

- **SISSIz cannot be seeded.** Splitting null replicates across machines gives
  different simulated alignments per host. That is statistically fine — they
  are all draws from the same null — but the run is then only reproducible by
  keeping every `null_pool/pool.stk`. To reproduce a calibration exactly, pin
  the pool with `--null-pool`.
- **Do not merge arms from different images.** If you rebuild the image, the
  conda solve may pick up newer versions. Record `results/versions.yaml` from
  each host and check they agree before pooling arms into one calibration.

For a single manuscript figure, the safest split is **one input alignment per
machine**, not arms of one alignment across machines: each run is then
internally consistent and carries its own `versions.yaml`.
