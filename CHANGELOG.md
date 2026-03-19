# Changelog

All notable changes to this project should be documented in this file.

The format is based on Keep a Changelog, and this project aims to use Semantic Versioning.

## [Unreleased]

### Added
- Installable Python package via `pyproject.toml`
- MIT `LICENSE` file
- Packaged workflow runner CLI: `rnaconsnake-run`
- Python replacements for the former Perl helper tools
- Helper tool CLIs: `rnaconsnake-split-stockholm`, `rnaconsnake-remove-gaponly`, and `rnaconsnake-strip-aln`
- Pytest-based regression tests for helper tools and CLI behavior
- End-to-end workflow smoke test with a stubbed `RNALalifold`
- Curated real-derived fixtures for Stockholm regression coverage
- GitHub Actions CI
- Developer docs for local virtualenv usage and workflow entry points

### Changed
- Workflow helper tools now run via packaged Python modules
- Repository workflow file is now tracked as lowercase `snakefile`
- Root `snakefile` and `config.yaml` are now the workflow source of truth for packaging
- Optional post-processing helper moved to `scripts/pp_RNALalifold.sh`
- `RNALalifold` remains an external runtime dependency

## [0.1.0] - 2026-03-19

### Added
- Initial packaged release of RNAConSnake
- Snakemake workflow distribution with Python helper tools
- Editable-install developer workflow
- Release checklist and Makefile targets
