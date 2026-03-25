# Changelog

All notable changes to this project should be documented in this file.

The format is based on Keep a Changelog, and this project aims to use Semantic Versioning.

## [Unreleased]

### Changed
- `SUMMARY_FIELDS` is now defined once in `workflow_helpers` and imported by both
  `legacy_postprocess` and `export_bundle`, eliminating a previously diverged duplicate
- `refold_firstseq` and `refold_firststruc` fields removed from summary outputs and export bundles
- `run_checked` now uses `contextlib.ExitStack` for guaranteed file handle cleanup

## [0.2.0] - 2026-03-24

### Added
- Installable Python package via `pyproject.toml`
- AGPL-3.0-or-later `LICENSE` file
- Packaged workflow runner CLI: `rnaconsnake-run`
- Python replacements for the former Perl helper tools
- Python helper module for legacy post-processing metric extraction and summary generation
- Python replacement for `alifold_maxcovar.pl`
- Helper tool CLIs: `rnaconsnake-split-stockholm`, `rnaconsnake-remove-gaponly`, and `rnaconsnake-strip-aln`
- Pytest-based regression tests for helper tools and CLI behavior
- End-to-end workflow smoke test with a stubbed `RNALalifold`
- Curated real-derived fixtures for Stockholm regression coverage
- GitHub Actions CI
- Developer docs for local virtualenv usage and workflow entry points

### Changed
- Workflow helper tools now run via packaged Python modules
- Snakemake workflow now implements the full legacy post-processing chain
- Repository workflow file is now tracked as lowercase `snakefile`
- Root `snakefile` and `config.yaml` are now the workflow source of truth for packaging
- Optional shell-based post-processing helper retained outside the default packaged runtime interface
- ViennaRNA and legacy post-processing commands remain external runtime dependencies

## [0.1.0] - 2026-03-19

### Added
- Initial packaged release of RNAConSnake
- Snakemake workflow distribution with Python helper tools
- Editable-install developer workflow
- Release checklist and local developer command guidance
