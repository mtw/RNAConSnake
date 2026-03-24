# RNAConSnake Architecture

## Purpose

RNAConSnake discovers local RNA structure candidates from aligned nucleotide features, computes per-candidate conservation and covariation metrics, and exports structured scientific results for downstream consumption.

## System Boundary

RNAConSnake is responsible for:

- workflow orchestration
- external scientific tool invocation
- candidate preprocessing and analysis
- structured summary generation
- scientific artifact generation
- export bundle generation

RNAConSnake is not responsible for:

- downstream data processing
- interactive data browsing

## High-Level Flow

The core workflow is:

1. Accept one input Stockholm alignment.
2. Run `RNALalifold` for one or more `maxbpspan` values.
3. Split the multi-record output into one candidate alignment per predicted element.
4. Preprocess each candidate alignment:
   - preserve original split alignment
   - remove gap-dominated sequences
   - filter ambiguous `N`/`n` content
   - remove redundant sequences
5. Analyze each cleaned candidate:
   - reformat to Clustal
   - run `RNAz`
   - run `alifoldz`
   - run `RNAalifold`
   - normalize `RNAalifold` outputs across versions
   - derive EPS/PDF structure figures
   - run constrained refolding
   - extract covariation metrics
   - optionally run `R-scape`
   - optionally build covariance models
6. Aggregate per-candidate metrics into summary records.
7. Emit structured outputs:
   - per-candidate JSON
   - text summary log
   - CSV summary
   - Markdown summary
8. Optionally export a bundle conforming to the public RNAConSnake export contract.

## Main Components

### Workflow Layer

Files:

- [`snakefile`](snakefile)
- [`config.yaml`](config.yaml)

Responsibilities:

- define Snakemake rules and dependencies
- coordinate stage execution
- connect stage outputs to manifests and summaries
- keep active pipeline behavior explicit

The Snakefile should remain focused on workflow topology. Reusable config and path logic belongs in Python modules rather than inline workflow code.

### CLI Layer

Files:

- [`src/rnaconsnake/cli.py`](src/rnaconsnake/cli.py)
- [`rnaconsnake-run`](rnaconsnake-run)

Responsibilities:

- expose the main workflow runner
- validate high-level CLI arguments
- check external runtime dependencies
- resolve packaged or repository workflow files
- invoke Snakemake
- optionally trigger export-bundle generation after a successful run

The repository-local launcher exists so the workflow can be run from a checkout without activating a virtual environment.

### Workflow Helper Layer

File:

- [`src/rnaconsnake/workflow_helpers.py`](src/rnaconsnake/workflow_helpers.py)

Responsibilities:

- parse and normalize workflow settings
- define canonical candidate paths
- provide workflow-facing file and manifest helpers
- centralize output normalization logic shared by multiple rules

This module is the preferred home for reusable workflow plumbing that would otherwise be duplicated in the Snakefile.

### Scientific Tool Layer

Files under [`src/rnaconsnake/tools/`](src/rnaconsnake/tools):

- `split_stockholm.py`
- `remove_gaponly.py`
- `strip_aln.py`
- `alifold_maxcovar.py`
- `legacy_postprocess.py`
- `verify_run_consistency.py`

Responsibilities:

- provide Python implementations for legacy shell-script behavior
- keep deterministic preprocessing and metric extraction logic testable
- isolate data-oriented postprocessing from workflow orchestration

Important note:

- `legacy_postprocess.py` now writes structured summary outputs only
- it no longer generates HTML

### Export Layer

File:

- [`src/rnaconsnake/export_bundle.py`](src/rnaconsnake/export_bundle.py)

Responsibilities:

- transform a completed run into a stable export bundle
- emit:
  - `manifest.json`
  - `features.csv`
  - `candidates.csv`
  - `artifacts.csv`
  - optional Markdown summaries
- preserve relative paths and artifact registration required by downstream tools

This module defines the public structured export boundary of the package.

## Data Model

### Workflow-Internal Data

The workflow is organized around per-candidate files keyed by:

- `wlen`
- `file`

Canonical path construction is centralized in `CandidatePaths` within [`workflow_helpers.py`](src/rnaconsnake/workflow_helpers.py).

### Export Data

When optional export-bundle generation is enabled, RNAConSnake emits a dataset-centric
structured bundle built around:

- dataset
- feature
- candidate
- artifact

This external model is intentionally more stable than the workflow's internal file
layout and rule topology.

## Output Classes

RNAConSnake produces:

- raw and cleaned candidate alignments
- Clustal reformatted alignments
- `RNAalifold` text outputs
- alignment and structure figures in PS/EPS/PDF
- refold summaries
- covariation metrics
- optional `R-scape` outputs
- optional covariance model outputs
- per-candidate JSON summaries
- summary log / CSV / Markdown
- optional structured export bundles

## Output Formats

RNAConSnake stores its outputs as regular files in the workflow output directory.

The main output formats are:

- Stockholm alignments: `.stk`
- Clustal alignments: `.aln`
- plain-text logs and tool outputs: `.log`, `.out`, `.txt`
- structured per-candidate records: `.json`
- tabular summaries: `.csv`
- narrative summaries: `.md`
- scientific figures and structure plots: `.ps`, `.eps`, `.pdf`, and optionally `.png`
- optional covariance model artifacts: Infernal CM files and related calibration outputs

At the summary/reporting layer, the primary human- and machine-readable outputs are:

- per-candidate JSON summary files
- per-run CSV summary tables
- per-run Markdown summaries

These are the intended stable structured outputs for downstream consumption. HTML is not generated by RNAConSnake.

## Testing Strategy

Tests live in:

- [`tests/test_python_tools.py`](tests/test_python_tools.py)

Current coverage emphasizes:

- CLI behavior
- workflow smoke tests with fake external tools
- deterministic preprocessing regressions
- metric extraction behavior
- export bundle generation
- run-consistency checks for deterministic workflow stages

The preferred testing approach is:

- small deterministic fixtures
- fake external tools for workflow integration
- explicit regression tests for legacy behavior replication

## Design Principles

- keep scientific computation separate from presentation
- favor structured outputs over implicit filesystem interpretation
- centralize path and config logic
- preserve compatibility with legacy shell behavior when scientifically relevant
- treat external tool quirks as explicit normalization logic
- keep workflow orchestration conservative and maintainable

## Downstream Boundary

Downstream consumers should depend on the public export bundle and should not depend on RNAConSnake internal directories, rule names, or ad hoc parsing of workflow outputs.

The intended boundary is:

- RNAConSnake computes and exports
- downstream tools consume the exported bundle or structured outputs

## Near-Term Refactoring Priorities

- continue shrinking Python `run:` logic inside the Snakefile
- further separate summary aggregation from report formatting helpers
- keep export-bundle generation aligned with the formal export specification
- avoid introducing HTML or presentation logic into this repository
