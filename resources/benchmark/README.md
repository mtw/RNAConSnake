# Benchmark truth files

Repository data (not generated): curated coordinates of known RNA elements for scoring the screen against ground truth.

**Alignment-specific**: coordinates meaningful only relative to alignments in `examples/`. Useless without them.

| file | alignment | status |
| --- | --- | --- |
| `jevg_3utr_elements.tsv` | `examples/all_JEVG_3UTR.relabel.stk` | curated, 9 elements; **training set** — config locked on this |
| `denvg_3utr_elements.tsv` | `examples/all_DENVG_3UTR.relabel.stk` | curated, 4 elements; **held-out test** — evaluated once; do not retune |
| `flavivirus_elements.tsv` | (none yet) | stub; broad alignment 30% too gappy, 5' third unscreenable |

## Format

| column | meaning |
| --- | --- |
| `element_id` | unique identifier |
| `element_class` | xrRNA, DB, sHP, 3SL, SL, hairpin, … |
| `alignment` | row scored if matches `benchmark_alignment` arg |
| `start`, `end` | 1-based alignment column range |
| `notes` | provenance, caveats |

Lines beginning with `#` are curator comments.

## Running

```bash
RNAcs --input-alignment examples/all_JEVG_3UTR.relabel.stk \
  --output-dir run_dir --null-arm sissiz --null-replicates 100 --benchmark --cores all
```

**Always report recovery with null baseline** (overlap easily satisfied by chance):

```bash
python -m rnaconsnake.tools.benchmark \
  --truth resources/benchmark/jevg_3utr_elements.tsv \
  --qvalues run_dir/results/calibration/qvalues.tsv \
  --null-loci run_dir/arms/null_000/generated_files/summary/len_100/RNAConSnake.nr.csv \
  ... one --null-loci per arm ... \
  --min-overlap-fraction 1.0 --output recovery.tsv
```

At 0.5 overlap: JEV real arm recovers 9/9, null arms recover 8.5/9 (test close to vacuous). Requiring full containment: 8/9 vs 5.1/9 null baseline.

## Curating

`benchmark_scaffold` mechanical step: structure-annotated alignment (`#=GC SS_cons` from mLocARNA) → consensus domains with coordinates filled, labels blank:

```bash
python -m rnaconsnake.tools.benchmark_scaffold \
  --alignment examples/all_JEVG_3UTR.relabel.stk \
  --alignment-id jevg_3utr --loci run_dir/generated_files/summary/len_100/RNAConSnake.nr.csv \
  --output my_scaffold.tsv
```

Flags curator decisions: unreported domains, wide domains (may bracket multiple elements), tiny domains (easy to miss/invent). Never guesses labels (e.g., which is xrRNA1); all rows emerge as `TBD`. Refuses files with placeholders.

**Curate against plots, not coordinates:** fold each domain full-span and inspect structure (best-scoring window often one lobe of two-lobed element; half dumbbell ≠ dumbbell):

```bash
python -m rnaconsnake.tools.fold_region \
  --alignment examples/all_JEVG_3UTR.relabel.stk --start 425 --end 497 \
  --output-dir review/DB1 --label DB1
```

**Known limitations:**
- Pseudoknots unmodelled (ViennaRNA folds xrRNAs as two hairpins, lower scores; truncate pseudoknot region to recover three-way junction)
- Very short windows unreliable z-scores (JEV sHP: −21.1 at 18 nt, −6.3 folded full-domain; prefer domain-span)
