# Architecture

RNAConSnake discovers local RNA structure candidates from aligned sequences, computes per-candidate metrics, and exports structured results.

## Scope

**In:** workflow orchestration, external tool invocation, preprocessing, analysis, structured output.  
**Out:** downstream processing, interactive browsing, HTML rendering.

## Pipeline

1. Accept Stockholm alignment; run `RNALalifold` per window size
2. Split output into one candidate per predicted element
3. Preprocess: remove gap-heavy and ambiguous sequences
4. Analyze: `RNAz`, `alifoldz`, `RNAalifold`, refold, covariation
5. Aggregate metrics and emit: JSON, CSV, Markdown, figures

## Components

| layer | file | purpose |
| --- | --- | --- |
| workflow | `snakefile`, `config.yaml` | rule topology, default config |
| CLI | `src/rnaconsnake/cli.py`, `RNAcs` | argument parsing, dependency check |
| helpers | `src/rnaconsnake/workflow_helpers.py` | paths, manifests, output normalization |
| tools | `src/rnaconsnake/tools/` | testable preprocessing, metrics, postprocessing |
| export | `src/rnaconsnake/export_bundle.py` | stable structured export contract |

## Outputs

Files: Stockholm (`.stk`), Clustal (`.aln`), structured records (`.json`, `.csv`, `.md`), figures (`.ps`, `.eps`, `.pdf`).

Primary outputs for consumption: per-candidate JSON, per-run CSV and Markdown. Optional structured export bundle. No HTML.

## Testing

Small synthetic fixtures, fake external tools, regression tests. Split by area in `tests/`.

## Design

- Separate scientific computation from presentation
- Favor structured outputs over filesystem interpretation
- Centralize path and config logic
- Preserve deterministic behavior for Stockholm parsing, manifests, candidates, export
- Keep export-bundle generation aligned with formal spec
