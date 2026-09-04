# Changelog

All notable changes to this project should be documented in this file.

The format is based on Keep a Changelog, and this project aims to use Semantic Versioning.

## [0.3.0] - 2026-08-28

### Added
- Null-model calibration arm: `--null-arm` with SISSIz (default) or `rnazRandomizeAln`; empirical q-values in `results/calibration/`
- De-replication of overlapping windows: `--dereplicate` with methods `containment`, `substructure`, `overlap`, `none`
- Positive-control benchmark: `--benchmark` scores recovery against curated elements; `--null-loci` provides null baseline
- Alignment screenability report: `generated_files/alignment_screenability.tsv` for per-window occupancy and identity
- Analysis aids outside the DAG: `fold_region`, `threshold_sweep`, `sensitivity_envelope`
- `tools/refold`: Python reimplementation of ViennaRNA's `refold.pl` using bindings only
- AlifoldZ zero-variance patch: `container/patches/alifoldz-zero-variance.patch` applied by `prepare-context.sh`
- Reproducibility artifacts: `results/versions.yaml`, `CITATION.cff`, container pinning SISSIz and `alifoldz.pl`
- Export bundle schema 1.1.0: locus and q-value columns; `dereplication`/`calibration` manifest blocks
- Input format detection: `.stk` uses `-f S`, `.aln` uses `-f C` in RNALalifold invocation
- Image published to `ghcr.io` on every `dev` push and version tag; Apptainer support documented
- R-scape included in container (`rscape=2.0.4.a`)
- SISSIz pinned to release 0.2.0; `--check-deps` enforces and reports version
- Export bundle schema 1.2.0: constrained refolds exported; `input_alignment_type` recorded
- `RNAcs --export-overwrite` for explicit directory replacement; `threshold_sweep --representative` flag

### Changed
- Refold: `refold.pl` and `RNAfold` no longer dependencies; Python module required, must match ViennaRNA binary version
- Command-line tool: `RNAcs` (alias `rnaconsnake-run` retained)
- Minimum Python: 3.11 (required by Snakemake 8)
- Defaults: `alifoldz_seed` set; `dereplicate.max_container_width` 120; `strip_aln` upper-cases sequences
- `{candidate}.refold.json` → `generated_files/consensus/{candidate}.consensus.json`; `extract-refold` → `extract-consensus`
- Linting and formatting: Ruff on CI and in `pyproject.toml` (110-column width)
- Tests: 196 tests at 92% coverage; CI includes container build with real toolchain
- Test suite split by area under `tests/` with shared helpers in `tests/helpers.py`
- `results/versions.yaml` written with `yaml.safe_dump`; PyYAML declared as dependency
- Benchmark recovery table: `reciprocal_overlap` moved beside `overlap_fraction`
- `remove_gaponly --infmt` defaults to `stockholm`

### Fixed
- Lower-case alignments: no longer hang workflow or produce garbage z-scores; normalized once in `strip_aln`
- Calibrated runs: SISSIz reproducibility via `--null-pool` pinning; `alifoldz_seed` pinning per-candidate shuffles
- AlifoldZ subprocess: no longer writes stray `alirna.ps` to run root; runs in private scratch directory
- `fold_region`: no longer writes `alirna.ps` to caller's directory
- Interleaved Stockholm blocks: `parse_stockholm_records` no longer truncates; `#=GR` annotation merged correctly
- `split_stockholm --accession`: now works as documented; defaults to `#=GF ID`
- `--config null=<non-mapping>`: reports error message instead of AttributeError
- `consensus_mfe`: now records consensus MFE, not mean single-sequence MFE (re-export affected bundles)
- `verify_run_consistency`: correctly compares calibrated runs when null-model arm enabled
- `alifold_maxcovar`: exits with status 0 on success
- `threshold_sweep`: respects configured representative rule
- `sensitivity_envelope`: no longer crashes on alignments with no subsets
- `RNAcs --export-bundle`: refuses to overwrite existing directories without `--export-overwrite`
- Dependency preflight: checks tools at paths resolved through `tools:` config and respects conditional flags
- `--benchmark` without null arm: fails immediately with explanation; with arm: includes null survivor counts
- AlifoldZ ranking: consistent direction in de-replication, summary reports, and export bundles
- `results/calibration/summary.json`: includes clustering parameters that affect q-values
- Curated truth files: all packaged, resolvable through `benchmark_truth` config
- `--benchmark`: adds to default targets instead of replacing them

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
