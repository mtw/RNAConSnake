# Code Review Handoff

Guidance for comprehensive pipeline review focusing on logic, production safety, hidden bugs, and test/documentation completeness.

## Pipeline Logic & Correctness

**Workflow topology (snakefile)**
- [ ] Rule DAG structure: real arm + null arms traverse identical rules via `A()` helper. Any arm-specific branches invalidate calibration.
- [ ] Checkpoint lookups: `checkpoint_wildcards()` used correctly in downstream rules.
- [ ] Null-pool pinning: `--null-pool` correctly loads pool from prior run and reproduces null-arm output.
- [ ] `null.method: none` inert: output layout and file paths byte-identical to pre-calibration runs.

**De-replication logic (tools/dereplicate.py, tools/loci.py)**
- [ ] Containment rule: windows correctly grouped when `a.start <= b.start and b.end <= a.end`.
- [ ] Container guard: windows beyond `max_container_width` (120 nt) prevent absorption of short windows by long ones.
- [ ] Grouping transitivity: A+B, B+C → A,B,C same locus.
- [ ] Representative selection: stable tie-breaking (rnazprob desc, alifoldzscore desc, name alpha) produces identical representatives across runs.
- [ ] Locus coordinates: `start`/`end` union of all members' ranges.

**Calibration math (tools/calibration.py)**
- [ ] Empirical q-values: correct percentile-rank calculation against null score distribution (real score ≥ null scores).
- [ ] FDR calculation: two-stage mode correctly conditions AlifoldZ on RNAz survivors.
- [ ] Collapse-ratio divergence: >20% real/null divergence correctly flagged in summary.
- [ ] Filter funnel: counting logic applies de-replication before cascade filters (not after).

**Metric extraction (tools/legacy_postprocess.py, tools/alifold_maxcovar.py)**
- [ ] RNAz parsing: `rnazprob`, `rnazstructcons`, `rnazgc` extracted correctly.
- [ ] AlifoldZ parsing: z-score correctly extracted; NA handling for zero-variance shuffles (via patch).
- [ ] Covariation metrics: `rscape_covary_count` correctly parsed from R-scape output.
- [ ] Case normalization: `strip_aln` handles all sequences to upper; `alifoldz.pl` and `refold.pl` receive upper-case input.

**Stockholm/Clustal conversion (tools/split_stockholm.py, tools/strip_aln.py)**
- [ ] Stockholm parsing: correctly preserves `#=GC` and sequence metadata; handles comments and blank lines.
- [ ] Clustal reformatting: `esl-reformat` invocation correct; output format matches RNAz expectations.
- [ ] Gap removal: sequences >50% gaps correctly removed; ambiguous `N` content correctly filtered.
- [ ] Redundancy removal: duplicate sequences correctly identified and removed.

## Production Safety

**Error handling & edge cases**
- [ ] Missing input alignment: error message is clear and actionable.
- [ ] Empty alignments: handled gracefully (raises error, not silent failure).
- [ ] Single-sequence alignments: rejected (minimum 3 sequences enforced).
- [ ] Very short windows (< 20 nt): handled (may give unreliable z-scores but don't crash).
- [ ] Pseudoknots: correctly noted as unmodelled in benchmark docs; users warned in `fold_region` output.
- [ ] Zero-variance shuffles: handled by alifoldz.patch (records NA, not division-by-zero crash).

**File I/O safety**
- [ ] Overwrite guards: output directory collision detection (warn or error, don't silently overwrite).
- [ ] Temp files: created in `{output_dir}/temp` with cleanup on success; retained on error for debugging.
- [ ] Permission errors: propagate clearly (not swallowed in try/except).
- [ ] Disk space: no unbounded temp accumulation in null-arm runs.

**Configuration validation**
- [ ] `--maxbpspan`: rejected if ≤0 or >10000.
- [ ] `--null-replicates`: rejected if ≤0 when null-arm is enabled; warning if >1000 (slow).
- [ ] `--dereplicate`: rejected if not in {containment, substructure, overlap, none}.
- [ ] `--null-arm`: rejected if not in {sissiz, rnazRandomizeAln, none}.
- [ ] Version check: `--check-deps` correctly reports external tool versions and warns on incompatibility.

**Reproducibility**
- [ ] SISSIz unseeded: reproducibility documented; users directed to `--null-pool` for exact reproduction.
- [ ] Version pinning: `results/versions.yaml` written with all toolchain versions on all runs; users warned not to merge arms across images.
- [ ] Conda solver non-determinism: `environment.yaml` uses pinned versions; container locks entire environment.

## Hidden Bugs & Errors

**Common pitfalls (check all)**
- [ ] Off-by-one errors in coordinate extraction (alignment columns are 1-based; Python lists 0-based).
- [ ] Float comparison: any `==` comparisons with floats should use epsilon tolerance.
- [ ] Uninitialized variables: all variables assigned before use.
- [ ] Resource leaks: file handles closed (via `with` or explicit close); processes reaped.
- [ ] Regex escaping: special characters escaped if used in file path interpolation.
- [ ] Unicode: alignment sequences correctly handled as bytes or strings consistently.

**Workflow-specific pitfalls**
- [ ] Wildcard expansion: Snakemake wildcard constraints correct and complete (avoid infinite expansions).
- [ ] Rule inputs/outputs: all rules list all their inputs and outputs (no implicit file creation).
- [ ] Checkpoints: correct use of `checkpoint_output()` and `checkpoint_wildcards()` — no stale file assumptions.
- [ ] Symlinks: if using symlinks in output, Snakemake correctly tracks them (usually does; verify if output is mounted).
- [ ] Circular dependencies: no rule outputs feed back into its own inputs (should be obvious but check multi-arm setups).

**Data corruption risks**
- [ ] Alignment modification: any in-place modification to input alignment caught and documented.
- [ ] Score truncation: float scores retained to full precision (not rounded prematurely).
- [ ] Missing field handling: empty `rnaz_prob` or `alifoldz_zscore` correctly handled as NA, not 0 or ignored.
- [ ] Parallel I/O: if multiple Snakemake jobs write to same directory, serialization enforced (or each writes to unique subdir).

## Documentation Completeness

**User-facing docs (README.md, docs/usage.md, docs/pipeline_summary.md)**
- [ ] Quick start runnable without external research.
- [ ] All CLI flags documented with examples.
- [ ] Output files clearly described.
- [ ] Calibration workflow step-by-step (what is `--null-arm`, why use it, how to interpret q-values).
- [ ] De-replication methods explained (when to use each).
- [ ] Known limitations documented (pseudoknots, z-score reliability at small sample sizes, screenability).
- [ ] Container usage clear (Apptainer `--cleanenv` requirement explained).

**Developer docs (AGENTS.md, CONTRIBUTING.md, docs/development.md)**
- [ ] Setup instructions complete (venv, pip install, external tool requirements).
- [ ] Development loop clear (test → lint → check-deps cycle).
- [ ] Invariants documented (arm traversal, no HTML, de-replication counts, config locked).
- [ ] Git hygiene rules stated (no run directories, venv, caches committed).
- [ ] Packaging notes present (root snakefile/config.yaml copied at build time).

**Architecture docs (ARCHITECTURE.md)**
- [ ] Purpose and scope clear.
- [ ] Component boundaries defined.
- [ ] Data model (workflow-internal vs export) explained.
- [ ] Output formats and locations specified.
- [ ] Design principles stated.

**Inline code comments**
- [ ] WHY (non-obvious constraints, invariants, workarounds) explained; WHAT omitted (code is self-documenting).
- [ ] Workarounds linked to issues/PRs when relevant.
- [ ] Deprecated code paths marked with deadline or removal version.

## Test Completeness

**Test coverage (tests/ directory)**
- [ ] CLI: argument parsing, version output, dependency checking, error messages.
- [ ] Workflow: smoke test with fake external tools; minimal synthetic fixtures.
- [ ] De-replication: all grouping methods tested; tie-breaking determinism verified.
- [ ] Calibration: empirical q-value calculation, two-stage filtering, collapse-ratio flagging.
- [ ] Preprocessing: gap removal, case normalization, redundancy removal.
- [ ] Metric extraction: RNAz/AlifoldZ/covariation parsing against real tool outputs.
- [ ] Export bundle: structure, schema version, manifest correctness.
- [ ] Null-model arm: `null.method: none` produces byte-identical output to pre-calibration layout.
- [ ] Project hygiene: no HTML, templates, or rendering code in repo.

**Regression tests**
- [ ] Deterministic outputs on fixed inputs: runs against benchmark sets (`jevg_3utr_elements.tsv`, `denvg_3utr_elements.tsv`) produce expected recovery.
- [ ] Run consistency: `verify_run_consistency` compares two completed runs byte-for-byte (excluding stochastic metrics like alifoldz z-scores).
- [ ] Version output: `results/versions.yaml` schema stable; changes recorded in export schema version.

**Edge cases tested**
- [ ] Empty alignments, single-sequence alignments, very small alignments (3 sequences).
- [ ] Very wide/narrow windows (`--maxbpspan 20`, `--maxbpspan 10000`).
- [ ] All de-replication methods including `--dereplicate none`.
- [ ] All null-arm methods including `--null-arm none`.
- [ ] Large numbers of replicates (`--null-replicates 1000`).
- [ ] Alignments where cascade filter rejects all candidates.
- [ ] Zero-variance shuffle handling in AlifoldZ.

**Test maintainability**
- [ ] Fixtures small and focused (CI-friendly).
- [ ] No brittle hardcoded paths (use `Path()` and `conftest.py` fixtures).
- [ ] Test names descriptive (test purpose clear from name).
- [ ] Mocks/fakes minimal and realistic (fake external tools produce realistic output format).

## Checks to Perform

```bash
# CI-enforced checks (should all pass)
pytest -v                          # full test suite
ruff check .                       # linting
ruff format --check .              # formatting
RNAcs --check-deps                 # dependency verification

# Manual verification
python -m rnaconsnake.tools.verify_run_consistency run_a run_b  # determinism check
git log --oneline -20              # recent commits and messages
git diff master..HEAD              # full diff against main branch
```

## Sign-off Criteria

✓ All logic checks pass (pipeline, de-replication, calibration, metrics)  
✓ No production safety issues found  
✓ No hidden bugs detected in code review  
✓ User-facing documentation complete and accurate  
✓ Developer documentation guides new contributors  
✓ Tests cover core paths and edge cases  
✓ Deterministic outputs verified  
✓ CI passes (pytest, ruff, check-deps)
