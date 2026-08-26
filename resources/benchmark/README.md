# Benchmark truth files

Curated coordinates of known RNA elements, so the calibrated screen can be
scored on a positive control: does it recover known structure, and at what
q-value?

These are **repository data, not generated artefacts**. Coordinates are
alignment-specific — they mean something only relative to the exact alignment
in `examples/`, and are useless without it.

## Files

| file | alignment | status |
| --- | --- | --- |
| `jevg_3utr_elements.tsv` | `examples/all_JEVG_3UTR.relabel.stk` | **curated, 9 elements** — the training set |
| `denvg_3utr_elements.tsv` | `examples/all_DENVG_3UTR.relabel.stk` | **curated, 4 elements** — the held-out test set |
| `flavivirus_elements.tsv` | (none yet) | uncurated stub for the broad mosquito-borne set |

### `jevg_3utr_elements.tsv` — training set

Nine elements of the JEV-group 3'UTR: xrRNA1, HP1, SL-III, xrRNA2, HP2, DB1,
DB2, sHP, 3SL. Two are reported here for the first time — **HP1 and HP2**,
small hairpins immediately 3' of xrRNA1 and xrRNA2 — and are flagged as such in
the file.

The screen configuration (representative rule and cascade thresholds) was
chosen against this set; see
[`../../docs/configuration_lock.md`](../../docs/configuration_lock.md).

### `denvg_3utr_elements.tsv` — held-out test set

Four elements of the dengue-group 3'UTR: DB1, DB2, sHP, 3SL. Curated **after**
the configuration was locked, and evaluated once without tuning. That is the
only out-of-sample number the study has, so it must stay held out: do not
adjust thresholds against results from this set.

Its header records what was deliberately **excluded** — three scaffold domains
examined and judged not to be xrRNAs. Omitting a row you do not believe is
correct; including it would inflate recovery.

### `flavivirus_elements.tsv` — stub

Element names and classes for the broad mosquito-borne set, with `TBD`
coordinates. Retained as a template. The broad alignment is a poor benchmark
substrate: 30% of its columns are too gappy to fold a consensus, and its 5'
third — where xrRNA1 sits — is unscreenable. See the screenability report
written by every run.

## File format

| column | meaning |
| --- | --- |
| `element_id` | unique identifier |
| `element_class` | `xrRNA`, `DB`, `sHP`, `3SL`, `SL`, `hairpin`, … |
| `alignment` | identifier passed as `benchmark_alignment`, so only matching rows are scored |
| `start` | 1-based **alignment column** of the first position |
| `end` | 1-based **alignment column** of the last position |
| `notes` | provenance, caveats, references |

Lines beginning with `#` are comments and carry the curation record.

## Running the benchmark

```bash
RNAcs --input-alignment examples/all_JEVG_3UTR.relabel.stk \
      --output-dir run_dir --null-arm sissiz --null-replicates 100 \
      --benchmark --cores all
```

**Always report recovery with a null baseline.** Reported loci cover much of an
alignment, so an overlap test is easily satisfied by chance:

```bash
python -m rnaconsnake.tools.benchmark \
  --truth resources/benchmark/jevg_3utr_elements.tsv \
  --qvalues run_dir/results/calibration/qvalues.tsv \
  --null-loci run_dir/arms/null_000/generated_files/summary/len_100/RNAConSnake.nr.csv \
  ... one --null-loci per arm ... \
  --min-overlap-fraction 1.0 \
  --output recovery.tsv
```

On the JEV group this matters a great deal. At the default
`--min-overlap-fraction 0.5` the real arm recovers 9/9 and the **null arms
recover 8.5/9** — the test is close to vacuous, and the report says so.
Requiring full containment gives 8/9 against a null baseline of 5.1/9.

## Curating a new truth file

`rnaconsnake.tools.benchmark_scaffold` does the mechanical half: given a
structure-annotated alignment (`#=GC SS_cons`, as produced by mLocARNA) it
decomposes the consensus into element-scale domains and writes a truth file
with **coordinates filled in and labels blank**.

```bash
python -m rnaconsnake.tools.benchmark_scaffold \
  --alignment    examples/all_JEVG_3UTR.relabel.stk \
  --alignment-id jevg_3utr \
  --loci run_dir/generated_files/summary/len_100/RNAConSnake.nr.csv \
  --output my_scaffold.tsv
```

It flags what a curator must decide about: consensus domains the screen did
**not** report, domains wider than 100 nt that may bracket several elements,
and domains under 20 nt that are easy to miss and easy to invent.

It will not guess which domain is xrRNA1 — that is a judgement about biology.
Every row emerges as `TBD`, and the benchmark **refuses a file containing
placeholder labels**, so an unedited scaffold cannot be mistaken for a curated
one.

### Curate against plots, not coordinates

Fold each domain over its **full span** and inspect the structure:

```bash
python -m rnaconsnake.tools.fold_region \
  --alignment examples/all_JEVG_3UTR.relabel.stk \
  --start 425 --end 497 --output-dir review/DB1 --label DB1
```

This matters. The screen reports the best-*scoring* window of a locus, which
for a two-lobed element is often one lobe — and half a dumbbell does not look
like a dumbbell. DB1 was rejected by the cascade as a 42 nt fragment scoring
z = −1.8; folded across its true 73 nt it scores −2.4 and passes.

### Known limitations worth recording in `notes`

- **Pseudoknots are unmodelled.** ViennaRNA cannot represent them, so xrRNAs
  fold as two hairpins and score below what their true structure warrants.
  Truncating xrRNA1's pseudoknot partner (columns 210–219) recovers the
  three-way junction, which is how the xrRNA1 note was established.
- **Very short windows give unreliable z-scores.** The JEV sHP scored −21.1 as
  an 18 nt window and −6.3 folded over its domain. Prefer the domain-span
  value.
