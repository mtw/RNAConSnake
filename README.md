[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19210864.svg)](https://doi.org/10.5281/zenodo.19210864)

# RNAConSnake

RNAConSnake is a workflow for discovering and prioritizing locally conserved RNA secondary structure candidates from multiple-sequence alignments.

It combines `RNALalifold` candidate discovery with per-candidate filtering, conservation scoring, covariation analysis, structure prediction, and structured result export. RNAConSnake is intended as the scientific computation layer: it produces candidate alignments, figures, metrics, and summary tables, but it does not generate browser-style HTML reports.

RNAConSnake is released under the GNU Affero General Public License v3 or later (AGPL-3.0-or-later).

## What Goes In

- one input multiple-sequence alignment in Stockholm format
- one or more `maxbpspan` window sizes
- an installed external RNA analysis toolchain available on `PATH`

An example input alignment is provided at [`examples/aln1.stk`](examples/aln1.stk).

## What Comes Out

RNAConSnake writes all outputs into a user-selected run directory. The main result classes are:

- split and cleaned candidate alignments in Stockholm and Clustal format
- per-candidate `RNAalifold` outputs and structure/alignment figures
- `RNAz`, `alifoldz`, max-covariation, and refold-derived metrics
- optional `R-scape` outputs
- per-candidate JSON summary records
- per-run summary files:
  - `RNAConSnake.log`
  - `RNAConSnake.log.csv`
  - `RNAConSnake.md`
- optional structured export bundles for downstream consumers

For a workflow-oriented overview, see [`docs/pipeline_summary.md`](docs/pipeline_summary.md).

## Quickstart

Install in a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

Check external dependencies:

```bash
rnaconsnake-run --check-deps
```

Run the workflow:

```bash
rnaconsnake-run --input-alignment /path/to/input.stk --output-dir /path/to/run_dir --cores all
```

If you want to use the repository checkout without activating the virtual environment, the repository-local launcher is:

```bash
./rnaconsnake-run --version
```

## Documentation

- User guide: [`docs/usage.md`](docs/usage.md)
- Development notes: [`docs/development.md`](docs/development.md)
- Pipeline summary: [`docs/pipeline_summary.md`](docs/pipeline_summary.md)
- Architecture: [`ARCHITECTURE.md`](ARCHITECTURE.md)
- Contributing: [`CONTRIBUTING.md`](CONTRIBUTING.md)

## Runtime Dependencies

RNAConSnake expects the external RNA analysis toolchain to be installed separately. The core workflow currently requires these commands on `PATH`:

- `RNALalifold`
- `esl-reformat`
- `RNAz`
- `alifoldz.pl`
- `RNAalifold`
- `ps2eps`
- `epstopdf`
- `refold.pl`
- `RNAfold`
- `magick`

Optional branches may additionally require tools such as `R-scape`, `cmbuild`, and `cmcalibrate`.

## Development And Testing

The repository includes a pytest-based regression suite and a Snakemake workflow smoke test built around fake external tools. Run:

```bash
pytest -v
```

Detailed developer workflow notes are in [`docs/development.md`](docs/development.md).
