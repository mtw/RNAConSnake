# De-replication Of Overlapping Windows

This document specifies exactly what RNAConSnake does about overlapping
RNALalifold windows, so that results can be interpreted and reproduced without
reading the source.

## The problem

`RNALalifold` reports every locally stable consensus structure it finds in a
sliding window. A single real element is therefore reported many times over:
once at close to its true extent, and again as a series of shorter windows
covering its individual stable helices. Those fragments are themselves
thermodynamically stable and phylogenetically conserved, so they score well in
their own right, and a raw ranked candidate list ends up dominated by pieces of
the same element.

Concretely, on a 20-sequence mosquito-borne flavivirus 3'UTR alignment at
`-L 100`, the terminal 3'SL region alone produced six windows, four of them
scoring well:

| window | width | `rnazprob` | `alifoldzscore` |
| --- | --- | --- | --- |
| `RC_100_0001_aln_712_799` | 88 | 0.999984 | −2.3 |
| `RC_100_0001_aln_719_794` | 76 | 0.999980 | −2.2 |
| `RC_100_0001_aln_739_778` | 40 | 0.969731 | −0.9 |
| `RC_100_0001_aln_745_759` | 15 | 0.356011 | −1.2 |
| `RC_100_0001_aln_754_777` | 24 | 0.000003 | 0.5 |
| `RC_100_0001_aln_755_772` | 18 | 0.000006 | NA |

De-replication groups windows that describe the same element and reports one
representative per group.

## What it produces

De-replication runs on every run; there is no flag to switch it on.

```
generated_files/summary/len_<wlen>/RNAConSnake.log.csv   # every window, unchanged
generated_files/summary/len_<wlen>/RNAConSnake.nr.csv    # one row per locus
generated_files/summary/len_<wlen>/RNAConSnake.md        # both, in one report
```

The per-window table is never modified. `RNAConSnake.nr.csv` has the columns

```
locus_id  locus_start  locus_end  n_windows  members  <then all per-window summary columns>
```

where the per-window columns are those of the chosen representative. `members`
is a `;`-separated, alphabetically sorted list of every window in the locus, so
nothing is silently discarded — anything collapsed can be looked up in the full
table.

Rows are sorted by descending representative `rnazprob`, then by `locus_start`.

`RNAConSnake.md` carries both tables in one report, non-redundant block first:

```markdown
# RNAConSnake Summary: len_200

- Loci (non-redundant): 6
- Windows (all): 18
- De-replication method: `containment`

## Non-redundant candidates
...one row per locus, members last...

## All windows
...every RNALalifold window, unchanged...
```

The de-replicated block leads because that is the list to read; the per-window
table follows so nothing is hidden. There is no separate non-redundant report.

## Procedure

Given all windows of one run for one `maxbpspan` value:

1. **Parse coordinates.** Each candidate basename ends in `_aln_<start>_<end>`;
   these are 1-based, inclusive **alignment column** coordinates.
2. **Parse the consensus structure** (`alifold_consstruc`, the `#=GC SS_cons`
   line from `RNAalifold`) into a set of base pairs, and shift it by the
   window's start so that all pairs are in absolute alignment columns.
   `(`, `<`, `[` and `{` are each matched against their own partner;
   unbalanced input yields fewer pairs rather than an error.
   Used by the `substructure` method only.
3. **Sort** candidates by `(start, end, name)`, so the grouping is
   deterministic and independent of input order.
4. **Group** by single-linkage clustering under the chosen method's pairwise
   predicate (below). Grouping is *transitive*: if A groups with B and B with
   C, then A, B and C are one locus.
5. **Order loci** by `(locus_start, locus_end)` and number them
   `<label>_0001`, `<label>_0002`, … where `<label>` is `len<wlen>`.
   `locus_start` and `locus_end` are the **union** of the members' spans.
6. **Pick a representative** per locus (below).

### Grouping predicates

Two windows are placed in the same locus when:

| `dereplicate.method` | predicate |
| --- | --- |
| `containment` (default) | one interval is nested inside the other: `a.start <= b.start and b.end <= a.end`, subject to the container guard below |
| `substructure` | at least `pair_containment` (default 0.9) of one window's base pairs are also formed by the other |
| `overlap` | the two intervals share at least `locus_min_overlap` (default 1) alignment columns |
| `none` | never — every window is its own locus |

Note what `containment` deliberately does **not** do: two windows that merely
overlap stay separate. In the example above, `RC_100_0001_aln_694_713`
(`rnazprob` 0.888, z −3.8) overlaps `..._712_799` by two columns but is not
nested in it, and remains its own candidate.

### The container guard

RNALalifold reports long windows as well as short ones, and a long window can
span several genuinely distinct elements. Without a guard, one 200-nt window
over three 60-nt hairpins absorbs all three and the run reports a single locus.

`containment` therefore restricts which windows may act as containers:

- a window **at most `max_container_width`** wide (default 120) is on the scale
  of a real element and may absorb anything nested inside it;
- a **wider** window may only absorb windows covering at least
  `container_min_coverage` (default 0.8) of it, so near-duplicate long windows
  still collapse into each other rather than each becoming its own locus;
- `max_container_width: 0` disables the guard.

Set `max_container_width` to the largest element you expect. Flavivirus 3'UTR
elements — xrRNAs, dumbbells, the 3'SL — are roughly 70–100 nt, hence the 120
default. On a JEV-group 3'UTR alignment at `-L 200` the guard takes the run
from 4 loci (7 RNAz-passing windows collapsed into 4) to 14 loci with 7
distinct RNAz-passing elements, while the long composite windows still collapse
into one locus instead of eight.

### Representative selection

The representative is the member with the largest key, comparing in order:

1. `rnazprob`, larger is better
2. `alifoldzscore`, **more negative** is better
3. `maxcovarval`, larger is better
4. `alilen`, larger is better — so a tie between a fragment and its parent
   resolves towards the parent
5. the candidate name, as a final deterministic tie-break

Non-numeric values (`NA`, empty) rank worst in every position, so a candidate
whose AlifoldZ crashed or was skipped never wins on that field.

Note that the representative's own span is in its `wbn`, and can be narrower
than `locus_start`–`locus_end`.

## Choosing a method

Measured on the flavivirus 3'UTR alignment above, at `-L 200` (18 windows):

| method | loci |
| --- | --- |
| `containment` | 6 |
| `substructure` | 8 |
| `overlap` | 3 |
| `none` | 18 |

**`containment` is the default** because it keys on coordinates and therefore
catches fragments even when a shorter window re-folds into slightly different
helices — which is common, and is exactly what defeats the structural test. In
this dataset the 76-nt window nested inside the 88-nt one shared only 18 of its
21 base pairs (86%), below the 0.9 subset threshold, so `substructure` keeps
both.

**`substructure`** is more conservative: it will not merge two genuinely
different structures that happen to be nested. Use it when you would rather
review a few fragments than risk collapsing distinct elements.

**`overlap`** is the classic RNAz-style locus merge. It chains: any two windows
sharing a single column join, so a run of pairwise-overlapping windows fuses
into one span longer than any real element. It is provided for comparability
with RNAz genome screens, not recommended as a default.

## Caveats

- **Grouping is transitive.** With `containment`, one very long RNALalifold
  window absorbs everything nested inside it. If two genuinely distinct
  elements both lie inside such a window, they become one locus. This is
  visible in the output: check `n_windows` and the `locus_start`–`locus_end`
  span. A span much wider than any plausible element is worth opening up, and
  `members` names what is inside it. `substructure` gives finer resolution.
- **`locus_start`/`locus_end` is the union of the members**, not the extent of
  the reported element.
- De-replication is a reporting and counting step. It does not re-fold, re-score
  or re-align anything; every value in `RNAConSnake.nr.csv` is copied verbatim
  from the representative's row in `RNAConSnake.log.csv`.
- Loci are formed **within** one `maxbpspan` value. Windows from different
  window sizes are never merged with each other.

## Relationship to the null-model calibration

The calibration arm counts the loci produced by this same procedure rather than
raw windows, so q-values describe the candidates the ranked table reports.
Overlapping windows are not independent, and counting them as if they were
would inflate both arms. `results/calibration/summary.json` records the method
used under `dereplication_method`, and reports the window-to-locus collapse
ratio per arm — if the arms de-replicate differently, the null count is
distorted and a warning is emitted.

## Configuration

```yaml
dereplicate:
  method: containment      # containment | substructure | overlap | none
  pair_containment: 0.9    # substructure only
calibration:
  locus_min_overlap: 1     # overlap only
```

Or on the command line:

```bash
RNAcs --input-alignment aln.stk --dereplicate substructure --cores all
```

## Running it standalone

The step is a plain command over the per-window table, so an existing run can
be re-grouped without repeating the pipeline:

```bash
python -m rnaconsnake.tools.dereplicate \
  --input  generated_files/summary/len_100/RNAConSnake.log.csv \
  --output regrouped.csv \
  --method substructure \
  --label  len100
```
