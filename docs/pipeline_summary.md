# RNAConSnake Pipeline Summary

## Overview

RNAConSnake is a Snakemake-based workflow for detecting and prioritizing locally conserved RNA secondary structure candidates from multiple-sequence alignments in Stockholm format. The pipeline combines a genome-wide screen with `RNALalifold` and a per-candidate post-processing stage that computes structure, covariation, and conservation metrics for each predicted local element. The workflow is intended for comparative RNA structure analysis in viral or other compact genomes, where candidate structural elements need to be identified systematically and ranked for downstream inspection.

## Input

The workflow takes a multiple-sequence alignment in Stockholm format as input. For each user-specified `maxbpspan` window size, `RNALalifold` is executed on the full alignment to identify locally stable consensus structures and to emit a multi-record Stockholm output containing one candidate element per predicted local structure.

## Candidate Generation And Preprocessing

The `RNALalifold` multi-record Stockholm output is split into one Stockholm file per candidate element. Each candidate alignment then undergoes a preprocessing step that preserves the original split alignment, removes sequences dominated by gaps, filters sequences containing more than a user-defined number of ambiguous `N`/`n` characters, and removes redundant sequences while retaining the cleaned Stockholm alignment for downstream analysis.

The two key preprocessing parameters are:

- `remove_gaponly_gapratio`: maximum tolerated fraction of gap characters per sequence
- `remove_gaponly_max_n`: maximum tolerated number of ambiguous `N`/`n` characters per sequence

## Per-candidate Analysis

Each cleaned candidate alignment is reformatted to Clustal and analyzed with `RNAz` and `alifoldz` to estimate structural conservation and thermodynamic support. In parallel, each candidate is passed to `RNAalifold` to compute a consensus secondary structure and associated alignment/structure plots. The workflow normalizes `RNAalifold` outputs across ViennaRNA versions, converts alignment PostScript plots to EPS and PDF, and retains the per-candidate Stockholm result produced by `RNAalifold`.

The `RNAalifold` result is subsequently reformatted and cleaned for constrained refolding. `refold.pl` and `RNAfold` are then used to derive a representative refolded sequence and structure under the consensus constraints. In addition, covariation support is quantified from the `RNAalifold` alignment output using a Python reimplementation of the legacy `alifold-maxcovar.pl` logic, which reports the maximum covariation score and the corresponding count summary for each candidate.

## Optional R-scape Analysis

The workflow includes an optional `R-scape` branch that is disabled by default and can be enabled explicitly via configuration or CLI flag. When enabled, `R-scape` is run on each per-candidate `RNAalifold_results.stk` file. The pipeline parses the generated `.power` output and records the number of base pairs reported as significantly covarying. For candidates with positive covariation support, the corresponding `R-scape` structural PDF is retained for interactive inspection in the final dashboard.

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

These per-candidate summaries are rendered into three report formats for each `maxbpspan` run:

- CSV
- Markdown
- interactive HTML

The HTML report supports numeric sorting and filtering and links directly to the generated structure plots. When `R-scape` is enabled and a candidate has positive covariation support, the `R-scape` PDF can be previewed interactively from within the dashboard.

## Output Structure

The workflow writes all generated data into a user-selected output directory. Major result classes are organized into stage-specific subdirectories, including:

- `Lalifold/` for the initial `RNALalifold` outputs and split candidate files
- `generated_files/stk/` for cleaned per-candidate Stockholm alignments
- `generated_files/rnaalifold/` for per-candidate `RNAalifold` outputs and structure plots
- `generated_files/rscape/` for optional `R-scape` outputs
- `generated_files/summary/` for final summary tables and dashboards

## Intended Use

RNAConSnake is designed as a reproducible candidate-discovery and prioritization workflow rather than a single-score classifier. It is best used to generate a ranked set of local RNA structure candidates that can then be inspected manually, filtered by user-defined thresholds, or taken forward for more detailed comparative, structural, or experimental validation.
