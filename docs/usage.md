# Usage

## Installation

Venv:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

From repo checkout:
```bash
./RNAcs --version
```

## Dependencies

```bash
RNAcs --check-deps
RNAcs --check-deps --rscape
RNAcs --show-tool-paths
```

## Basic run

```bash
RNAcs --input-alignment aln.{stk,aln} --output-dir run_dir --cores all
```

## Common options

Override window sizes:
```bash
RNAcs --input-alignment aln.stk --maxbpspan 100 --maxbpspan 200 --cores all
```

Enable R-scape:
```bash
RNAcs --input-alignment aln.stk --rscape --cores all
```

Allow RNAz shuffling (default: disabled):
```bash
RNAcs --input-alignment aln.stk --rnaz-shuffle --cores all
```

Conservative output reuse (only rebuild if inputs changed):
```bash
RNAcs --input-alignment aln.stk --conservative --cores all
```

## Output

Main results per `maxbpspan`:

```
results/versions.yaml                            # toolchain versions
generated_files/summary/len_*/RNAConSnake.log.csv        # all windows
generated_files/summary/len_*/RNAConSnake.nr.csv         # one per locus
generated_files/summary/len_*/RNAConSnake.md             # both tables
generated_files/alignment_screenability.tsv     # per-window occupancy
```

Calibration (with `--null-arm`):

```
results/calibration/qvalues.tsv          # per-locus q-values
results/calibration/summary.json         # FDR, method, parameters
results/calibration/funnel.tsv           # filter stages per arm
```

Export (with `--export-bundle`):

```
<bundle>/manifest.json                   # schema and provenance
<bundle>/candidates.csv                  # all loci with q-values
```

## De-replication

Each run writes both full and de-replicated tables. Change method:

```bash
RNAcs --input-alignment aln.stk --dereplicate substructure --cores all
```

Methods: `containment` (default), `substructure`, `overlap`, `none`.

Details: [`dereplication.md`](dereplication.md)

## Alignment screenability

Every run writes where the screen is blind (too gappy, too divergent):

```bash
column -ts$'\t' generated_files/alignment_screenability.tsv | less
```

Or on alignment alone:

```bash
python -m rnaconsnake.tools.alignment_report --alignment aln.stk \
  --output screenability.tsv --metadata screenability.json
```

## Null-model calibration

```bash
RNAcs --input-alignment aln.stk --null-arm sissiz \
  --null-replicates 100 --output-dir run_dir --cores all
```

Methods: `sissiz` (default), `rnazRandomizeAln`, `none`.

**Reproducibility:** SISSIz cannot be seeded. Pin a pool:

```bash
RNAcs --input-alignment aln.stk --null-arm sissiz \
  --null-pool /path/to/run/null_pool/pool.stk --output-dir run_dir \
  --cores all
```

Details: [`configuration_lock.md`](configuration_lock.md)

## Analysis aids

Fold and plot a span:

```bash
python -m rnaconsnake.tools.fold_region --alignment aln.stk \
  --start 100 --end 200 --output-dir review/span1 --label elem1
```

Sweep cascade thresholds:

```bash
python -m rnaconsnake.tools.threshold_sweep \
  --arm-input real:100:run/arms/real/generated_files/summary/len_100/RNAConSnake.log.csv \
  --arm-input null_000:100:run/arms/null_000/generated_files/summary/len_100/RNAConSnake.log.csv \
  ... one per arm ... \
  --output sweep.tsv
```

Measure detection floor:

```bash
python -m rnaconsnake.tools.sensitivity_envelope \
  --alignment aln_known.stk --manifest envelope/manifest.tsv \
  --output-dir envelope/subsets --sizes 3,4,5,6
```
