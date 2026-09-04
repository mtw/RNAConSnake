# De-replication of overlapping windows

RNALalifold reports every locally stable structure. One element may appear many
times: once at true extent, and as shorter windows over its
helices. De-replication groups these and reports one representative per locus.

## Output

Runs on every run; no flag to disable. Three files per `maxbpspan`:

```
generated_files/summary/len_<wlen>/RNAConSnake.log.csv   # all windows
generated_files/summary/len_<wlen>/RNAConSnake.nr.csv    # one per locus
generated_files/summary/len_<wlen>/RNAConSnake.md        # both tables
```

`RNAConSnake.nr.csv` columns:

```
locus_id  locus_start  locus_end  n_windows  members  <...per-window columns>
```

`members` lists every window in the locus (alphabetically sorted, `;`-separated).

## Grouping methods

Two windows group when:

| `--dereplicate` | rule |
| --- | --- |
| `containment` | one nested in the other: `a.start <= b.start and b.end <= a.end` |
| `substructure` | ≥90% of one's base pairs formed by the other |
| `overlap` | share ≥1 alignment column |
| `none` | never |

`containment` is the default. Substructure grouping uses default 0.9 base-pair
overlap threshold; **Container guard:** Long windows cannot absorb short ones
beyond `max_container_width` (default 120 nt). Prevents a 200-nt window over
three hairpins from collapsing them into one locus.

## Representative selection

Per locus, by: descending `rnazprob`, then `alifoldzscore` (more negative
better), then candidate name. Deterministic across runs.

## Caveats

- Grouping is **transitive**: A groups with B, B with C → A, B, C are one locus
- `locus_start`, `locus_end` are the union of all members
- Collapse ratio compared across arms (real vs null); >20% divergence flags the summary

## Standalone usage

Run de-replication outside the workflow:

```bash
python -m rnaconsnake.tools.dereplicate \
  --input RNAConSnake.log.csv \
  --output RNAConSnake.nr.csv \
  --method containment \
  --label my_run
```

Thresholds: default 0.9 pair-containment overlap for substructure method,
default 1 locus-minimum-overlap for grouping.
