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
./RNAcs --version
```

`rnaconsnake-run` is kept as a deprecated alias. It automatically prefers `./.venv/bin/python` when present and otherwise falls back to `python3`.

## Dependency checks

Check the standard runtime dependencies:

```bash
RNAcs --check-deps
```

Check the dependency set including optional `R-scape` support:

```bash
RNAcs --check-deps --rscape
```

The check follows the run's own configuration: a tool redirected through
`tools:` is checked at the command configured there, and the branches a run
will actually take — `do_rscape`, `null.method`, whether set by flag, by
`--config`, or in the config file — bring in their own dependencies. Pass the
same `--configfile` you will run with to check that run:

```bash
RNAcs --check-deps -- --configfile /path/to/config.yaml
```

To inspect which external tool commands RNAConSnake is configured to use, and where their executables resolve on the current machine:

```bash
RNAcs --show-tool-paths
```

## Basic run

```bash
RNAcs --input-alignment /path/to/my_alignment.{stk,aln} --cores all
```

To write all outputs into a dedicated working directory:

```bash
RNAcs --input-alignment /path/to/my_alignment.{stk,aln} --output-dir /path/to/run_dir --cores all
```

## Common options

Override the default `maxbpspan` values:

```bash
RNAcs --input-alignment /path/to/my_alignment.{stk,aln} --maxbpspan 150 --maxbpspan 250 --cores all
```

Enable optional per-hit `R-scape` analysis:

```bash
RNAcs --input-alignment /path/to/my_alignment.{stk,aln} --rscape --cores all
```

RNAz defaults to `-d -n`, which disables alignment shuffling. To explicitly allow shuffling:

```bash
RNAcs --input-alignment /path/to/my_alignment.{stk,aln} --rnaz-shuffle --cores all
```

If you want Snakemake to reuse existing outputs more conservatively and avoid reruns caused only by workflow code changes:

```bash
RNAcs --conservative --cores all
```

## De-replicating overlapping windows

RNALalifold reports every locally stable structure it finds, so a single real
element comes back many times over: once near its true extent, and again as a
series of shorter windows covering its individual stable helices. Those
fragments are themselves stable and conserved, so they score well too, and a
raw ranked list ends up dominated by pieces of the same element.

Each run therefore writes a non-redundant table alongside the full one:

```
generated_files/summary/len_<wlen>/RNAConSnake.log.csv   # every window
generated_files/summary/len_<wlen>/RNAConSnake.nr.csv    # one row per locus
generated_files/summary/len_<wlen>/RNAConSnake.md        # both, loci block first
```

`RNAConSnake.nr.csv` carries the best-scoring window of each locus, plus
`locus_id`, `locus_start`, `locus_end`, `n_windows` and a `members` column
naming every window that was collapsed into it. Nothing is discarded: the full
per-window table is untouched.

It runs on every run; there is no flag to switch it on. To change how windows
are grouped:

```bash
RNAcs --input-alignment aln.stk --dereplicate substructure --cores all
```

| `--dereplicate` | rule |
| --- | --- |
| `containment` (default) | a window nested inside another is a fragment of it; merely overlapping windows stay separate |
| `substructure` | stricter: the fragment's consensus base pairs must also be a subset of the other's |
| `overlap` | classic RNAz-style merge on any shared column; chains |
| `none` | report every window |

The same de-replication is the counting unit for the null-model calibration, so
q-values describe the candidates the ranked table reports.

**The exact procedure — grouping predicates, representative selection, output
columns and caveats — is specified in
[`dereplication.md`](dereplication.md).** Read it before interpreting
`n_windows` or a locus span.

## Is the alignment screenable?

Every run writes `generated_files/alignment_screenability.tsv`. A local
structure screen is blind in two different ways, and neither announces itself
in the results:

- **Too gappy.** Where most sequences are gaps there is no consensus to fold,
  so no candidate can be reported. This is indistinguishable, in the output,
  from a genuine absence of structure.
- **Too divergent.** Below roughly 55% mean pairwise identity RNAz loses
  discriminating power. Strong elements are still found — this is a caveat, not
  a verdict — so these regions are reported as `reduced_power`, not as blind.

The report gives per-window occupancy and identity plus explicit warnings. Read
it before concluding that a region contains no structure.

This matters most for alignments spanning a broad clade. A 20-sequence
mosquito-borne-flavivirus 3'UTR alignment has 30% of its columns too gappy to
fold, all of them at the 5' end, and the screen returns nothing there. Per-clade
alignments of the same region (JEV group, dengue group) are far more tractable —
11% and 8% unscreenable respectively — and recover elements across the whole
UTR. If elements are missing, check this report before tuning thresholds.

The report can also be produced for an alignment on its own, before committing
to a screen — useful when deciding whether a new alignment is worth the compute:

```bash
python -m rnaconsnake.tools.alignment_report \
  --alignment my_alignment.stk \
  --output screenability.tsv --metadata screenability.json
```

## Inspecting one region

The screen reports the best-*scoring* window of a locus, which need not be the
element's true extent. To fold and plot an arbitrary span with the workflow's
own RNAalifold parameters:

```bash
python -m rnaconsnake.tools.fold_region \
  --alignment examples/all_JEVG_3UTR.relabel.stk \
  --start 425 --end 497 \
  --output-dir review/DB1 --label DB1
```

It writes the consensus structure, the structure and alignment plots, and
RNAz/AlifoldZ recomputed on that exact span, so the two are directly
comparable. This is the tool to reach for when a reported candidate looks like
a fragment of something larger, or when checking whether a short window's
z-score survives folding over the full element.

## Null-model calibration arm

The calibration arm reruns the *same* rule DAG over simulated alignments that
carry no real structural signal, and uses the resulting score distributions to
attach an empirical FDR / q-value to every real candidate.

```bash
RNAcs --input-alignment aln.stk --output-dir run_dir --null-arm --cores all
```

Options:

| option | meaning |
| --- | --- |
| `--null-arm [METHOD]` | Enable the arm. Bare `--null-arm` uses `sissiz`. `rnazRandomizeAln` is faster and weaker; `none` disables it. |
| `--null-replicates N` | Replicate count. 10 gives an order-of-magnitude sanity check; 100 gives usable q-values with resolution floored at 1/100. |
| `--null-seed N` | Base seed. Per-arm seeds are derived deterministically from it. |
| `--no-two-stage` | Run AlifoldZ on every candidate instead of only on stage-one survivors. Much slower, but the FDR is then unconditional. |
| `--benchmark` | Also build the positive-control recovery table, with the null arms as its baseline. Requires the arm. |
| `--emit-versions` | Write `results/versions.yaml` even without the calibration arm. |

### What changes when the arm is on

Every pipeline output moves under `arms/<arm>/`, where `<arm>` is `real` or
`null_000`, `null_001`, ... The real arm is *not* special-cased: it traverses
exactly the same rules, starting from a copy (never a symlink) of the input
alignment made by `make_arm_alignment`. Any divergence between arms would
invalidate the calibration, so there is deliberately no way to express one.

Outputs:

- `results/calibration/funnel.tsv` - `n_in`/`n_out` per filter stage, per arm,
  plus a `null_mean` aggregate. The header records whether the FDR is
  conditional on stage-one survival.
- `results/calibration/qvalues.tsv` - one row per real **locus** with its
  scores and q-values.
- `results/calibration/score_dists.tsv` - the null score distributions.
- `results/calibration/summary.json` - method, replicates, seed, thresholds
  (including the clustering parameters that decide how many loci each arm
  reports, and therefore the q-values), collapse ratios, warnings, and the
  composite cascade FDR.

### Statistics

For each score `s` and threshold `t`:

```
FDR(t) = mean_over_replicates( N_null(s >= t) ) / N_real(s >= t)
```

The reported q-value is the monotone envelope `q(c) = min over t <= s(c) of
FDR(t)`, clamped to `[0, 1]`. It is computed separately for the RNAz class
probability and the AlifoldZ z-score, and once more for the composite cascade
(candidates surviving every filter) - the cascade number is the headline one.

### Things to watch

- **Overlapping windows are not independent.** RNALalifold emits overlapping
  candidates, so all counting is done on merged loci. Each locus is represented
  by one deterministically chosen window (best RNAz probability, then best
  AlifoldZ, then name), and that window's scores are used for every filter, in
  both arms. The window-to-locus collapse ratio is compared across arms; if it
  differs by more than `calibration.collapse_ratio_tolerance` (default 20%), a
  warning is written to `summary.json` and into the funnel header, because a
  divergent collapse would inflate the null count and bias the FDR in the
  flattering direction.
- **Two-stage mode makes the FDR conditional.** AlifoldZ shuffles internally
  and is the expensive branch, so by default it runs only on candidates that
  clear `calibration.stage1_rnaz_prob`. The resulting FDR is conditional on
  passing stage one, and is labelled as such in `summary.json` and in the
  funnel header. It must not be reported as an unconditional FDR. Use
  `--no-two-stage` for an unconditional estimate.
- **Locus merging is transitive.** See "De-replicating overlapping windows"
  above: the counting unit is the de-replicated locus, and with the default
  `containment` method a single long window can absorb everything nested inside
  it. The collapse ratio is reported per arm so this is visible.
- **The collapse-ratio guard is noisy at small locus counts.** It compares a
  ratio of counts, so on short alignments that yield only a handful of loci it
  can fire on sampling noise alone. Read it as a prompt to look at the funnel,
  not as an automatic verdict.
- **R-scape not run.** With `do_rscape: false` every candidate carries `NA`, so
  the R-scape filter would reject everything for reasons unrelated to the data.
  In that case it is reported in the funnel but excluded from the cascade, and
  `summary.json` records `rscape_evaluated: false` together with the
  `cascade_filters` actually applied.
- **R-scape power at small n.** At around 6 sequences and ~35% mean pairwise
  identity, R-scape returns almost no significant pairs in *either* arm. That
  is expected, not a bug. The filter stays in the cascade but is necessarily
  permissive at this sample size, and the funnel reports the counts as they are.
- **`null.seed` does not reproduce a SISSIz pool.** SISSIz exposes no seed
  option and draws its own from the clock, so two runs a second apart simulate
  different alignments. Only the `rnazRandomizeAln` backend is genuinely
  seeded. To reproduce a calibration, **pin the pool that produced it**:

  ```yaml
  "null":
    pool_file: /path/to/previous_run/null_pool/pool.stk
  ```

  A pinned pool is re-validated against the current alignment, so one taken
  from a different alignment is rejected rather than silently miscalibrating.
  The pool's SHA-256 is recorded in `summary.json` either way, so a rerun can
  be verified rather than assumed, and an unseeded backend adds a warning
  saying so.

## Choosing the cascade thresholds

Threshold choices are usually inherited rather than measured. Once a calibrated
run exists the question is answerable directly: at each candidate cutoff, how
many real loci survive, how many survive by chance, and what FDR results?

```bash
python -m rnaconsnake.tools.threshold_sweep \
  --arm-input real:100:run_dir/arms/real/generated_files/summary/len_100/RNAConSnake.log.csv \
  --arm-input null_000:100:run_dir/arms/null_000/generated_files/summary/len_100/RNAConSnake.log.csv \
  ... one --arm-input per arm ... \
  --reference resources/benchmark/<truth>.tsv \
  --output sweep.tsv
```

It reuses the per-arm tables a calibrated run already wrote, so a sweep costs
seconds rather than another screen. With `--reference` it adds how many
reference spans a survivor covers.

Read the recovery column as **sensitivity only**: a screen reporting one locus
across the whole alignment would recover everything. The FDR column is what
constrains the other side, and the two have to be read together. A threshold
that loses real elements *and* worsens the FDR is simply dominated.

Resolution is bounded by the replicate count: with 10 replicates the null mean
moves in steps of 0.1, so small FDR differences are noise. Sweep on the
replicate count you intend to report.

## Sensitivity envelope: is "nothing" evidence of absence?

A screen returning no candidates has two explanations — nothing is there, or
the alignment could not have shown it. The null arm bounds false positives but
says nothing about this floor. Before claiming that sparsity in one organism is
biological rather than methodological, measure where detection fails:

```bash
python -m rnaconsnake.tools.sensitivity_envelope \
  --alignment examples/<alignment with curated elements>.stk \
  --output-dir envelope/subsets --manifest envelope/manifest.tsv \
  --sizes 3,4,5,6 --identities 0.35,0.45,0.55,0.65
```

It writes subsets of the alignment spanning a range of sequence count and mean
pairwise identity. Subsetting drops sequences and never columns, so curated
element coordinates stay valid. Screen each subset and record where the known
elements stop being recovered; that boundary is the envelope.

**A caution on what is achievable.** Subsetting a small alignment barely lowers
identity, and often raises it — dropping the most divergent sequence makes the
remainder more similar. The JEV-group alignment (5 sequences, MPI 0.61) yields
subsets spanning only 0.54-0.68, which is nowhere near the divergence of a
broad filovirus set. Measuring an envelope that reaches low identity needs a
source alignment that is already diverse: the 20-sequence mosquito-borne set
yields subsets from 0.35 to 0.68, but its elements have to be curated first.
Check the achievable range before designing the experiment.

## Positive control benchmark

```bash
RNAcs --input-alignment benchmark.stk --output-dir run_dir --null-arm --benchmark --cores all
```

This scores `results/calibration/qvalues.tsv` against the curated truth file at
`resources/benchmark/flavivirus_elements.tsv` and writes
`results/benchmark/flavivirus_recovery.tsv`: one row per known element, whether
it was recovered, and at what q-value.

`--benchmark` requires the null arm, because the table it scores is the
calibrated one. Without `--null-arm` (or `null.method` in the config file) the
run stops before starting and says so.

### Recovery needs a null baseline

Reported loci typically cover most of an alignment, so an overlap test can be
satisfied by chance. The workflow therefore passes every null arm's locus
tables to the benchmark, and the report states how many elements those arms
also "recover". Run by hand, it looks like this:

```bash
python -m rnaconsnake.tools.benchmark \
  --truth resources/benchmark/<truth>.tsv \
  --qvalues run_dir/results/calibration/qvalues.tsv \
  --null-loci run_dir/arms/null_000/generated_files/summary/len_{100,200}/RNAConSnake.nr.csv \
  ... one --null-loci per arm ... \
  --min-overlap-fraction 1.0 \
  --output recovery.tsv
```

One `--null-loci` per arm, listing all of that arm's window lengths after it:
the tables under a single `--null-loci` pool into one baseline sample. Giving
each window length its own `--null-loci` would count it as a separate
replicate, and understate the baseline.

On a JEV-group screen this matters a great deal. At the default
`--min-overlap-fraction 0.5` the real arm recovers 9/9 and the *null* arms
recover 8.5/9 — the test is close to vacuous, and the report says so. Requiring
full containment gives 8/9 against a null baseline of 5.1/9, which is a claim
worth making.

`reciprocal_overlap` reports how much of the *locus* the element accounts for.
A locus far larger than the element still scores 1.0 on `overlap_fraction`,
which is how a screen reporting one huge locus can look perfectly sensitive.

The truth file ships with element names and classes but **uncurated
coordinates**, because alignment-column coordinates only mean something
relative to the exact benchmark alignment used. The benchmark refuses to run
until they are filled in; see `resources/benchmark/README.md`. Set
`benchmark_allow_uncurated: true` to emit `uncurated` rows instead while wiring
things up.

## Structured export

RNAConSnake can generate a structured export bundle after a successful run:

```bash
RNAcs \
  --input-alignment /path/to/my_alignment.stk \
  --output-dir /path/to/run_dir \
  --export-bundle /path/to/export_bundle \
  --cores all
```

You can also override export metadata, for example:

```bash
RNAcs \
  --input-alignment /path/to/my_alignment.stk \
  --output-dir /path/to/run_dir \
  --export-bundle /path/to/export_bundle \
  --export-feature-type other \
  --export-feature-label "Example aligned feature" \
  --export-dataset-label "Example dataset" \
  --cores all
```

### What the export bundle carries

`--export-bundle` writes schema **1.1.0**, which adds:

- per-candidate `locus_id`, `locus_start`, `locus_end`, `locus_window_count`,
  `is_representative` and `redundant_to`, so a consumer can show one row per
  locus and still reach the collapsed windows;
- per-candidate `q_rnaz`, `q_alifoldz`, `q_cascade` and `cascade_pass` when the
  null-model arm was run;
- `dereplication` and `calibration` blocks in `manifest.json`, including
  `fdr_conditional_on_stage_one`, so the conditionality of the FDR is legible
  without reading a table;
- the calibration tables and `results/versions.yaml` as dataset artifacts.

Every new column is optional, so a consumer written against 1.0.0 keeps working
and simply does not see the new fields. When the null arm is enabled the export
reads the real arm from `arms/real/`.

## Configuration

The default workflow configuration lives in [`../config.yaml`](../config.yaml).

De-replication settings:

```yaml
dereplicate:
  method: containment      # containment | substructure | overlap | none
  pair_containment: 0.9    # substructure method only
```

Null-arm and calibration settings:

```yaml
"null":
  method: none            # sissiz | rnazRandomizeAln | none
  replicates: 10
  seed: 20261101
  two_stage: true

calibration:
  rnaz_prob_threshold: 0.9
  alifoldz_threshold: -2.0
  rscape_min_pairs: 1
  stage1_rnaz_prob: 0.5
  locus_min_overlap: 1
  collapse_ratio_tolerance: 0.2
```

Note the quotes around `"null"`. An unquoted `null:` key is the YAML null
scalar, so the section would land under a `None` key. RNAConSnake accepts both
spellings, but the quoted form is the one to write.

`stage1_rnaz_prob` must not exceed `rnaz_prob_threshold`; otherwise the
reported RNAz filter would select candidates whose AlifoldZ score was never
computed, and the calibration step refuses to run.

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
rm -rf Lalifold generated_files arms null_pool results
```
