# RNAConSnake Pipeline Summary

## Overview

RNAConSnake is a Snakemake-based workflow for detecting and prioritizing locally conserved RNA secondary structure candidates from multiple-sequence alignments in Stockholm or Clustal format. The pipeline combines a genome-wide screen with `RNALalifold` and a per-candidate post-processing stage that computes structure, covariation, and conservation metrics for each predicted local element. The workflow is intended for comparative RNA structure analysis in viral or other compact genomes, where candidate structural elements need to be identified systematically and ranked for downstream inspection.

## Input

The workflow takes a multiple-sequence alignment in Stockholm (`.stk`) or Clustal (`.aln`) format as input. For each user-specified `maxbpspan` window size, `RNALalifold` is executed on the full alignment to identify locally stable consensus structures and to emit a multi-record Stockholm output containing one candidate element per predicted local structure. RNAConSnake detects the input format from the filename suffix and calls `RNALalifold` with `-f S` for Stockholm input or `-f C` for Clustal input.

## Alignment Screenability

Before any candidate is reported, each run characterises the input alignment
itself and writes `generated_files/alignment_screenability.tsv`. A local
structure screen is blind in two distinct ways that the candidate list cannot
distinguish from a genuine absence of structure. Where most sequences are gaps
there is no consensus to fold, so nothing can be reported at all; where mean
pairwise identity falls below roughly 55%, `RNAz` retains sensitivity to strong
elements but loses discriminating power. The report separates these — the first
is disqualifying, the second a caveat — and warns explicitly, so that regions
returning no candidates are not silently read as regions containing no
structure.

## Candidate Generation And Preprocessing

The `RNALalifold` multi-record Stockholm output is split into one Stockholm file per candidate element. Each candidate alignment then undergoes a preprocessing step that preserves the original split alignment, removes sequences dominated by gaps, filters sequences containing more than a user-defined number of ambiguous `N`/`n` characters, and removes redundant sequences while retaining the cleaned Stockholm alignment for downstream analysis.

The two key preprocessing parameters are:

- `remove_gaponly_gapratio`: maximum tolerated fraction of gap characters per sequence
- `remove_gaponly_max_n`: maximum tolerated number of ambiguous `N`/`n` characters per sequence

## Per-candidate Analysis

Each cleaned candidate alignment is reformatted to Clustal and analyzed with `RNAz` and `alifoldz` to estimate structural conservation and thermodynamic support. In parallel, each candidate is passed to `RNAalifold` to compute a consensus secondary structure and associated alignment/structure plots. The workflow normalizes `RNAalifold` outputs across ViennaRNA versions, converts alignment PostScript plots to EPS and PDF, and retains the per-candidate Stockholm result produced by `RNAalifold`.

The `RNAalifold` result is subsequently reformatted and cleaned for constrained refolding. A packaged Python tool (`rnaconsnake.tools.refold`) then derives, for every sequence, its own constraint from the consensus structure and folds it under that constraint through the ViennaRNA Python bindings. It replaces the `refold.pl | RNAfold -C` pipe the legacy workflow used, and reproduces it byte-for-byte. In addition, covariation support is quantified from the `RNAalifold` alignment output using a Python reimplementation of the legacy `alifold-maxcovar.pl` logic, which reports the maximum covariation score and the corresponding count summary for each candidate.

## De-replication Of Overlapping Windows

Because `RNALalifold` reports every locally stable structure it finds, a single
real element is typically reported both near its true extent and as a series of
shorter windows covering its individual stable helices. Those fragments are
themselves stable and conserved and therefore score well, so the raw ranked
candidate list is dominated by pieces of the same element.

The workflow therefore emits a non-redundant candidate table alongside the full
per-window one. Windows describing the same element are grouped, and the
best-scoring member of each group is reported together with the locus span, the
number of windows collapsed, and the identity of every member. The Markdown
report contains both tables, the non-redundant block first and the full
per-window table after it. The default
grouping rule is coordinate containment: a window nested inside another is a
fragment of it, while merely overlapping windows remain separate candidates. A
stricter structure-aware rule, which additionally requires the fragment's
consensus base pairs to be a subset of the other's, is available, as is the
classic overlap-based locus merge. De-replicated loci are also the counting
unit used by the null-model calibration. The procedure is specified in
[`dereplication.md`](dereplication.md).

## Optional Null-Model Calibration

The screen can be calibrated against a null model. When enabled, the same rule
DAG is additionally traversed by one or more simulated replicates of the input
alignment, generated either by `SISSIz` in alignment-simulation mode (the
default, preserving the tree, dinucleotide content and gap pattern) or by
column shuffling with `rnazRandomizeAln.pl`. Every pipeline output then lives
under `arms/<arm>/`, where `<arm>` is `real` or `null_000`, `null_001`, and so
on. Setting `null.method: none` or `replicates: 0` restores the uncalibrated
behaviour and output layout.

Because `RNALalifold` emits overlapping windows, which are not independent, the
calibration step first de-redundifies windows into merged loci and does all
counting on loci. It then computes, for each score and threshold, the mean
number of null survivors divided by the number of real survivors, and reports
the monotone envelope of that ratio as a q-value. Separate q-values are
produced for the `RNAz` class probability and the `alifoldz` z-score, plus a
single composite value for candidates surviving the whole cascade.

`alifoldz` shuffles internally and is the expensive branch, so by default it is
evaluated only on candidates that clear a cheap first-stage `RNAz` threshold,
in both arms. The resulting FDR is then conditional on stage-one survival, and
is labelled as such in the funnel header and in `summary.json`.

## Optional R-scape Analysis

The workflow includes an optional `R-scape` branch that is disabled by default and can be enabled explicitly via configuration or CLI flag. When enabled, `R-scape` is run on each per-candidate `RNAalifold_results.stk` file. The pipeline parses the generated `.power` output and records the number of base pairs reported as significantly covarying. For candidates with positive covariation support, the corresponding `R-scape` structural PDF is retained as a scientific artifact alongside the structured summaries.

## Optional Positive Control

The screen can be scored against curated coordinates of known elements, held in
`resources/benchmark/`. The benchmark reports, per element, whether a reported
locus recovers it and at what q-value.

Recovery is reported alongside a **null-arm baseline**: because reported loci
cover much of an alignment, an overlap criterion can be satisfied by chance, and
a recovery count without that baseline is uninterpretable. On the JEV-group
benchmark, requiring only 50% overlap yields 9 of 9 elements recovered while the
null arms recover 8.5 of 9; requiring full containment yields 8 of 9 against a
null baseline of 5.1. The report warns when the criterion is close to vacuous.

Truth files are repository data, curated by hand. A scaffolding tool drafts them
from a structure-annotated alignment, filling in coordinates but never assigning
element labels, and the benchmark refuses a file that still contains placeholder
labels so that an unedited scaffold cannot be mistaken for a curated one.

## Optional Covariance Model Construction

The pipeline can optionally build and calibrate Infernal covariance models for candidates that satisfy at least one of two user-configurable support criteria: an `RNAz` probability threshold or an `alifoldz` threshold. This step is disabled by default. When enabled, candidates passing the thresholds are converted into calibrated covariance models with `cmbuild` and `cmcalibrate`.

## Summary Outputs

For each candidate element, the workflow aggregates the following classes of measurements into a structured summary record:

- alignment size metrics such as number of sequences and alignment length
- covariation metrics from the max-covariation analysis
- optional `R-scape` covariation support
- `RNAz` probability
- `alifoldz` score
- consensus and refold-derived structural descriptors

These per-candidate summaries are rendered into three structured summary formats for each `maxbpspan` run:

- plain text log
- CSV
- Markdown

RNAConSnake itself does not render HTML dashboards or websites. Visualization and browser-style presentation are intentionally outside the scope of this repository.

## Output Structure

The workflow writes all generated data into a user-selected output directory. Major result classes are organized into stage-specific subdirectories, including:

- `Lalifold/` for the initial `RNALalifold` outputs and split candidate files
- `generated_files/stk/` for cleaned per-candidate Stockholm alignments
- `generated_files/rnaalifold/` for per-candidate `RNAalifold` outputs and structure plots
- `generated_files/rscape/` for optional `R-scape` outputs
- `generated_files/summary/` for final summary tables and structured reports

## Intended Use

RNAConSnake is designed as a reproducible candidate-discovery and prioritization workflow rather than a single-score classifier. It is best used to generate a ranked set of local RNA structure candidates that can then be inspected manually, filtered by user-defined thresholds, or taken forward for more detailed comparative, structural, or experimental validation.
