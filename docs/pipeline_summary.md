# Pipeline summary

RNAConSnake screens multiple sequence alignments for conserved RNA secondary structure.

## Inputs

- Multiple-sequence alignment (Stockholm `.stk` or Clustal `.aln`). At least 3 sequences.
- Window sizes (`--maxbpspan`, default 100 and 200 nt).

## Core pipeline

1. **RNALalifold:** consensus structure per window
2. **RNAz:** structure conservation; probability and SCI
3. **AlifoldZ:** shuffle-derived z-score
4. **Cascade filter:** RNAz probability ≥ 0.8 AND AlifoldZ z ≤ -1.5 (configurable)

## Key outputs

```
generated_files/summary/len_<wlen>/RNAConSnake.nr.csv      # one per locus
results/versions.yaml                                      # toolchain versions
```

With `--null-arm`:

```
results/calibration/qvalues.tsv                            # empirical q-values
results/calibration/summary.json                           # FDR, filter counts
```

With `--benchmark`:

```
results/calibration/flavivirus_recovery.tsv                # recovery vs q-value
```

## Pre-processing

- Sequences with >50% gaps removed
- Sequences with ambiguous `N` removed
- Case normalized to upper

## De-replication

RNALalifold reports overlapping windows. De-replication groups by containment
(default) and reports one per locus. Collapse ratio (real vs null) checked; >20%
divergence flags the summary. Details: [`dereplication.md`](dereplication.md).

## Null-model calibration

Runs the same DAG on shuffled alignments. Score distributions → empirical FDR
and q-values. Two-stage mode: AlifoldZ only on stage-one survivors (FDR then
conditional).

## Reproducibility

- All versions in `results/versions.yaml`
- SISSIz cannot be seeded; pin the pool with `--null-pool` to reproduce
- Different SISSIz builds = different null model

## Options

| flag | what |
| --- | --- |
| `--maxbpspan` | window size; repeat for multiple |
| `--dereplicate` | grouping: containment, substructure, overlap, none |
| `--null-arm` | calibration; sissiz, rnazRandomizeAln, none |
| `--null-replicates` | 10 for sanity check, 100 for publication |
| `--null-pool` | reuse previous pool for reproducibility |
| `--benchmark` | score recovery of curated elements |
| `--rscape` | per-candidate covariation |
| `--rnaz-shuffle` | allow RNAz shuffling (default: disabled) |
| `--export-bundle` | structured export at path |
| `--conservative` | only rebuild if inputs changed |
