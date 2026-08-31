# RNAcs container

Image pinning the whole toolchain, for running the calibrated screen off the
development machine. Nothing in the Dockerfile is architecture-specific: the
base image is multi-arch and the one tool built from source, SISSIz, is
compiled in-image, so the same context builds a native **linux/amd64** image on
an x86_64 host and a native **linux/arm64** one on Apple Silicon.

On a cluster with Apptainer but no Docker daemon, do not build at all — pull
the x86_64 image CI publishes. See [Running under Apptainer](#running-under-apptainer).

## Why a container

Most of the toolchain installs from conda-forge/bioconda, but three pieces do
not, and getting them right by hand on each machine is the failure mode this
image removes:

| piece | why it needs special handling |
| --- | --- |
| **SISSIz** | not in bioconda on any platform; compiled from source, pinned to release **0.2.0** ([https://github.com/mtw/SISSIz](https://github.com/mtw/SISSIz)) |
| **alifoldz.pl** | ships in the RNAz *source* tarball, not the conda package |
| **ps2eps**, **epstopdf** | Debian packages; not in conda-forge |

Refolding needs no vendored tool: it runs through the ViennaRNA Python module,
which must be importable in the image's environment.

R-scape **is** included (`rscape=2.0.4.a`), so `do_rscape: true` works in the
image without extending it, though the default configuration leaves it off.

## Build

```bash
cd container
./prepare-context.sh        # gathers SISSIz source + the two perl helpers
docker build -t rnacs:0.3.0 .
```

With no `--platform`, this builds natively for the host architecture, which is
what you want: an emulated build compiles SISSIz under QEMU, and SISSIz is what
simulates the null alignments. To cross-build anyway — an x86_64 image on
Apple Silicon, say — use `docker buildx build --platform linux/amd64` and
expect tens of minutes.

`prepare-context.sh` collects three things that are not on any package index
and must already be present on the build machine:

| piece | where to get it |
| --- | --- |
| SISSIz source | clone [https://github.com/mtw/SISSIz](https://github.com/mtw/SISSIz) at tag `0.2.0`; the build tree, not a binary |
| `alifoldz.pl` | the **RNAz source tarball** (`perl/` directory) — absent from the conda package |

`prepare-context.sh` applies
[`patches/alifoldz-zero-variance.patch`](patches/alifoldz-zero-variance.patch)
to whichever `alifoldz.pl` it is given, and says which it did. Without the
guard, a candidate whose shuffles all fold to the same energy makes the script
divide by zero and die with an unexplained Perl error; with it, the candidate
is recorded as unscored. A copy that already carries the guard is left alone,
so pointing `ALIFOLDZ` at either an upstream or an already-patched script
produces the same image.

`ps2eps` and `epstopdf` come from Debian inside the image and need no action.

Override the locations if yours differ:

```bash
git clone --branch 0.2.0 https://github.com/mtw/SISSIz ~/src/SISSIz

SISSIZ_SRC=~/src/SISSIz \
ALIFOLDZ=~/.local/share/RNAz/perl/alifoldz.pl \
./prepare-context.sh
```

The vendored sources are **not committed** — the image is rebuilt from your own
copies, so nothing third-party is redistributed through this repository.

## Run

```bash
docker run --rm \
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

`docker save` / `docker load` moves the built image without rebuilding, as long
as both hosts are the same architecture:

```bash
docker save rnacs:0.3.0 | gzip -1 | ssh <host> 'gunzip | docker load'
ssh <host> 'docker run --rm rnacs:0.3.0 --check-deps'
```

## Running under Apptainer

Cluster nodes typically have no Docker daemon and no root, so the image is not
built there — CI builds an x86_64 image, runs the real toolchain through it,
and pushes the *tested* image to ghcr.io. Pull that:

```bash
apptainer pull rnacs_0.3.0.sif docker://ghcr.io/mtw/rnaconsnake:0.3.0
apptainer run --cleanenv rnacs_0.3.0.sif --check-deps
```

A version tag publishes its version (`0.3.0`); every push to `dev` publishes
`dev-<short-sha>`. Untagged builds exist so an image can be tried on a cluster
before a release is tagged — pull the `dev-` tag matching the commit you want:

```bash
apptainer pull rnacs_dev.sif docker://ghcr.io/mtw/rnaconsnake:dev-1c94c92
```

A ghcr.io package is **private on first push**, so that `pull` fails with an
unhelpful `unauthorized` until the package's visibility is set to public (GitHub
→ the package page → Package settings). To keep it private instead,
authenticate on the cluster with a personal access token carrying
`read:packages`:

```bash
export APPTAINER_DOCKER_USERNAME=<github-user>
export APPTAINER_DOCKER_PASSWORD=<token>
```

`--check-deps` is worth running before anything long: it reports the ViennaRNA
and SISSIz versions the image resolved and verifies that the ViennaRNA binaries
and the Python module are one build.

```bash
apptainer run --cleanenv \
  --bind "$PWD/examples:/data:ro" \
  --bind "$PWD/runs/ebov:/work" \
  rnacs_0.3.0.sif \
  --input-alignment /data/all_JEVG_3UTR.relabel.stk \
  --output-dir /work \
  --maxbpspan 100 --cores 10 --no-progress \
  --null-arm sissiz --null-replicates 100
```

Four differences from the `docker run` above, all of which matter here:

- **`--cleanenv` is not optional.** Apptainer exports the host environment into
  the container by default, so a `module load` that set `PERL5LIB`,
  `PYTHONPATH` or a conda prefix will shadow the image's own toolchain — the
  image would no longer be pinning what it claims to pin. `--cleanenv` is what
  makes the run reproducible.
- **You run as yourself**, not root. The `getpwuid()` failure that `--user`
  causes under Docker does not arise: Apptainer binds the host `/etc/passwd`.
  Outputs in the bound run directory are owned by you rather than by root.
- **The image filesystem is read-only.** Everything the workflow writes goes to
  the bound `/work`, which is fine; if some tool insists on writing elsewhere
  in the image, add `--writable-tmpfs`.
- **The working directory is the host's, not the image's `/work`.** Always pass
  `--output-dir` explicitly, and bind the run directory rather than relying on
  the `WORKDIR`.

Create the run directory before binding it — Apptainer will not create a bind
target that does not exist. Under a scheduler, keep `--cores` at or below the
cores the job was allocated.

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
