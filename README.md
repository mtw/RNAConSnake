[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19210864.svg)](https://doi.org/10.5281/zenodo.19210864)

# RNAConSnake

RNAConSnake is a workflow for discovering and prioritizing locally conserved RNA secondary structure candidates from multiple-sequence alignments. The command-line tool is called `RNAcs`.

It combines `RNALalifold` candidate discovery with per-candidate filtering, conservation scoring, covariation analysis, structure prediction, and structured result export. RNAConSnake is intended as the scientific computation layer: it produces candidate alignments, figures, metrics, and summary tables, but it does not generate browser-style HTML reports.

RNAConSnake is released under the GNU Affero General Public License v3 or later (AGPL-3.0-or-later).

## What Goes In

- one input multiple-sequence alignment in Stockholm (`.stk`) or Clustal (`.aln`) format
- one or more `maxbpspan` window sizes
- an installed external RNA analysis toolchain available on `PATH`

An example input alignment is provided at [`examples/aln1.stk`](examples/aln1.stk).

## What Comes Out

RNAConSnake writes all outputs into a user-selected run directory. The main result classes are:

- split and cleaned candidate alignments in Stockholm and Clustal format
- per-candidate `RNAalifold` outputs and structure/alignment figures
- `RNAz`, `alifoldz`, max-covariation, and refold-derived metrics
- optional `R-scape` outputs
- per-candidate JSON summary records
- per-run summary files:
  - `RNAConSnake.log`
  - `RNAConSnake.log.csv` - every RNALalifold window
  - `RNAConSnake.nr.csv` - one row per locus after de-replication
  - `RNAConSnake.md` - both tables, non-redundant block first
- optional structured export bundles for downstream consumers

With the null-model calibration arm enabled, the run additionally produces:

- `results/calibration/funnel.tsv` - survivor counts per filter stage, per arm
- `results/calibration/qvalues.tsv` - per real candidate locus: scores plus an empirical q-value
- `results/calibration/score_dists.tsv` - null score distributions, for plotting
- `results/calibration/summary.json` - method, replicates, thresholds, and whether the FDR is conditional
- `results/versions.yaml` - the exact toolchain the run used

For a workflow-oriented overview, see [`docs/pipeline_summary.md`](docs/pipeline_summary.md).

## De-replicating overlapping windows

RNALalifold reports every locally stable structure, so one real element comes
back many times: once near its true extent, and again as shorter windows over
its individual stable helices. Those fragments score well in their own right,
which is what makes a raw ranked list hard to read.

Each run therefore also writes `RNAConSnake.nr.csv`, with one row per locus -
the best-scoring window, plus the span, the number of windows collapsed into
it, and a `members` column naming each one. The full per-window table is left
untouched. On a 20-sequence flavivirus 3'UTR alignment this collapses 22
windows to 7 loci. It is also the counting unit the calibration uses, so
q-values describe the candidates that actually get reported.

The exact procedure is specified in
[`docs/dereplication.md`](docs/dereplication.md).

## Null-model calibration

Without a null model, "stringent filtering" is not quantified anywhere: nothing
says how many candidates would survive the same cascade by chance. The
calibration arm runs the *identical* rule DAG over simulated alignments that
preserve phylogeny, base and dinucleotide composition, gap structure and mean
pairwise identity but carry no real structural signal, and turns the resulting
score distributions into empirical FDR / q-values.

```bash
RNAcs --input-alignment aln.stk --output-dir run_dir --null-arm --null-replicates 100 --cores all
```

Key properties:

- The arms are not parallel rule sets. Every rule downstream of alignment
  generation is traversed by both arms; the only difference is an
  `arms/<arm>/` path prefix, which is empty when the arm is off.
- Counting is on **merged loci**, not raw RNALalifold windows, because
  overlapping windows are not independent. The window-to-locus collapse ratio
  is compared across arms and a warning is recorded if they disagree.
- `null.method: none` (the shipped default) or `replicates: 0` reproduces the
  pre-calibration behaviour and output layout.

See [`docs/usage.md`](docs/usage.md) for the full option set.

## Is the alignment screenable?

A local structure screen is blind where most sequences are gaps: no consensus
can be folded, so no candidate is reported — indistinguishable in the output
from a genuine absence of structure. Every run therefore writes
`generated_files/alignment_screenability.tsv`, giving per-window occupancy and
mean pairwise identity, and warning where the screen could not have looked.

Read it before concluding that a region contains nothing. See
[`docs/usage.md`](docs/usage.md).

## Positive control

`resources/benchmark/` holds curated coordinates of known elements, and
`RNAcs --benchmark` reports which are recovered and at what q-value. Truth
files can be drafted from a structure-annotated alignment with
`rnaconsnake.tools.benchmark_scaffold`, but element labels are always assigned
by a curator, never by code.

Report recovery **with a null baseline** (`--null-loci`): reported loci cover
much of an alignment, so an overlap test can be satisfied by chance, and a
recovery count without that baseline is uninterpretable.

## Quickstart

Install in a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

Check external dependencies:

```bash
RNAcs --check-deps
```

Run the workflow:

```bash
RNAcs --input-alignment /path/to/input_alignment.{stk,aln} --output-dir /path/to/run_dir --cores all
```

With the null-model calibration arm, which attaches empirical q-values:

```bash
RNAcs --input-alignment aln.stk --output-dir run_dir \
      --null-arm sissiz --null-replicates 100 --cores all
```

100 replicates gives a q resolution of 0.01; 10 is an order-of-magnitude sanity
check. The arm is off by default, and `null.method: none` reproduces the
pre-calibration output layout exactly.

If you want to use the repository checkout without activating the virtual environment, the repository-local launcher is:

```bash
./RNAcs --version
```

`rnaconsnake-run` remains as a deprecated alias for `RNAcs`.

## Running in a container

A native linux/arm64 image pins the entire toolchain, including the three tools
that are not on any package index. This is the recommended way to run
reproducible calibration jobs.

```bash
# one-off: fetch the sources conda cannot supply
git clone https://github.com/mtw/SISSIz ~/src/SISSIz

cd container
SISSIZ_SRC=~/src/SISSIz ./prepare-context.sh
docker build --platform linux/arm64 -t rnacs:0.3.0 .
```

Run it, mounting inputs read-only and the run directory read-write:

```bash
docker run --rm --platform linux/arm64 \
  -v "$PWD/examples:/data:ro" \
  -v "$PWD/myrun:/work" \
  rnacs:0.3.0 \
  --input-alignment /data/all_JEVG_3UTR.relabel.stk \
  --output-dir /work \
  --maxbpspan 100 --cores 8 --no-progress \
  --null-arm sissiz --null-replicates 100
```

The entrypoint is `RNAcs`, so all options pass straight through. Full details,
including how to move the image to another machine, are in
[`container/README.md`](container/README.md).

## Documentation

- User guide: [`docs/usage.md`](docs/usage.md)
- Window de-replication: [`docs/dereplication.md`](docs/dereplication.md)
- Development notes: [`docs/development.md`](docs/development.md)
- Pipeline summary: [`docs/pipeline_summary.md`](docs/pipeline_summary.md)
- Contributing: [`CONTRIBUTING.md`](CONTRIBUTING.md)
- Container build: [`container/README.md`](container/README.md)

## Runtime Dependencies

RNAConSnake expects the external RNA analysis toolchain on `PATH`. Most of it
installs from conda-forge/bioconda (see [`environment.yaml`](environment.yaml));
three pieces do not and must be installed by hand:

| tool | source |
| --- | --- |
| `SISSIz` | <https://github.com/mtw/SISSIz> — build from source; required for the null-model arm |
| `alifoldz.pl` | the **RNAz source tarball** (`perl/`); absent from the conda package |
| `refold.pl` | the **ViennaRNA source tarball** (`src/Utils/`); absent from the conda package |

From conda: `RNALalifold`, `RNAalifold`, `RNAfold`, `RNAz`, `esl-reformat`,
`ps2eps`, `epstopdf`, `magick`. Optional branches additionally use `R-scape`,
`cmbuild` and `cmcalibrate`.

Check what is resolvable on the current machine:

```bash
RNAcs --check-deps
RNAcs --show-tool-paths
```

The [container](container/README.md) packages all of the above, including the
three that are not on any package index, and is the least painful way to get a
working installation.

**Note on `--use-conda`:** most rules use Snakemake `run:` directives, which
execute in the Snakemake process and therefore cannot take per-rule conda
environments. Use `environment.yaml` or the container instead.

## Development And Testing

The repository includes a pytest-based regression suite and a Snakemake workflow smoke test built around fake external tools. Run:

```bash
pytest -v
```

Detailed developer workflow notes are in [`docs/development.md`](docs/development.md).
