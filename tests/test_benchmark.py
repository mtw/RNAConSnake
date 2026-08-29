"""The positive control: curated truth files, the scaffolding that
drafts them, and the null baseline that makes a recovery count readable."""

from __future__ import annotations

import itertools
import subprocess
import sys
from pathlib import Path

import pytest

from tests.helpers import (
    BENCHMARK_TRUTH,
    PYTHON,
    _locus,
    _two_elements,
    read_text,
    subprocess_env,
)


def test_benchmark_scores_recovery_against_calibrated_loci(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, read_qvalues, read_truth, write_recovery

    truth_path = tmp_path / "truth.tsv"
    truth_path.write_text(BENCHMARK_TRUTH, encoding="utf-8")
    qvalues_path = tmp_path / "qvalues.tsv"
    qvalues_path.write_text(
        "\n".join(
            [
                "# fdr_conditional_on_stage_one\ttrue",
                "locus_id\twlen\tstart\tend\trnazprob\tq_rnaz\tq_alifoldz\tq_cascade\tcascade_pass",
                "len100_0001\t100\t90\t210\t0.97\t0\t0.01\t0\tyes",
                "len100_0002\t100\t400\t500\t0.5\t0.4\t0.5\tNA\tno",
                "",
            ]
        ),
        encoding="utf-8",
    )

    results = evaluate(
        read_truth(truth_path, "flavivirus_3utr"),
        read_qvalues(qvalues_path),
        min_overlap_fraction=0.5,
        allow_uncurated=False,
    )
    by_id = {row["element_id"]: row for row in results}
    assert by_id["xrRNA1"]["recovered"] == "yes"
    assert by_id["xrRNA1"]["best_locus"] == "len100_0001"
    assert by_id["xrRNA1"]["q_cascade"] == "0"
    assert by_id["xrRNA2"]["recovered"] == "no"

    out = tmp_path / "recovery.tsv"
    write_recovery(results, out)
    text = read_text(out)
    assert "# recovered\t1" in text
    assert "# curated\t2" in text


def test_benchmark_refuses_uncurated_truth_file(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import BenchmarkError, evaluate, read_truth

    truth_path = tmp_path / "truth.tsv"
    truth_path.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "xrRNA1\txrRNA\tflavivirus_3utr\tTBD\tTBD\tnot curated yet\n",
        encoding="utf-8",
    )
    truth = read_truth(truth_path)
    with pytest.raises(BenchmarkError, match="uncurated coordinates"):
        evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=False)

    results = evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=True)
    assert results[0]["recovered"] == "uncurated"


def test_shipped_flavivirus_truth_file_is_schema_valid() -> None:
    from rnaconsnake.tools.benchmark import TRUTH_COLUMNS, read_truth

    truth = read_truth(Path("resources/benchmark/flavivirus_elements.tsv"))
    assert {"xrRNA", "DB", "sHP", "3SL"} <= {element.element_class for element in truth}
    assert TRUTH_COLUMNS[0] == "element_id"
    # Coordinates are alignment specific and intentionally not shipped curated.
    assert all(not element.curated for element in truth)


def test_structural_domains_descend_through_long_range_pairs() -> None:
    """Long-range helices bracket several elements rather than defining one."""
    from rnaconsnake.tools.benchmark_scaffold import structural_domains

    # One long-range pair enclosing two element-scale hairpins.
    left = "(((((....)))))"  # 14 nt, 5 bp
    right = "((((......))))"  # 14 nt, 5 bp
    structure = "(" + left + "....." + right + ")"

    wide = structural_domains(structure, max_width=200, min_width=5, min_pairs=3)
    assert len(wide) == 1, "with a generous width the whole thing is one domain"

    narrow = structural_domains(structure, max_width=20, min_width=5, min_pairs=3)
    assert [(d.start, d.end) for d in narrow] == [(2, 15), (21, 34)]
    assert [d.n_pairs for d in narrow] == [5, 4]


def test_structural_domains_skip_trivial_helices() -> None:
    from rnaconsnake.tools.benchmark_scaffold import structural_domains

    assert structural_domains("((...))", max_width=50, min_width=5, min_pairs=3) == []
    assert structural_domains(".........", max_width=50) == []


def test_read_ss_cons_requires_a_reference_structure(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark_scaffold import read_ss_cons

    bare = tmp_path / "bare.stk"
    bare.write_text("# STOCKHOLM 1.0\nseqA ACGU\nseqB ACGA\n//\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no '#=GC SS_cons' line"):
        read_ss_cons(bare)

    annotated = tmp_path / "annotated.stk"
    annotated.write_text("# STOCKHOLM 1.0\nseqA ACGU\nseqB ACGA\n#=GC SS_cons <<>>\n//\n", encoding="utf-8")
    structure, n_seq = read_ss_cons(annotated)
    assert structure == "<<>>"
    assert n_seq == 2


def test_scaffold_marks_missed_domains_and_unexplained_loci() -> None:
    from rnaconsnake.tools.benchmark_scaffold import Domain, scaffold_rows

    domains = [Domain(start=10, end=60, n_pairs=12), Domain(start=200, end=260, n_pairs=14)]
    loci = [
        {"locus_id": "L1", "locus_start": "5", "locus_end": "70", "wbn": "w1", "rnazprob": "0.99"},
        {"locus_id": "L2", "locus_start": "400", "locus_end": "460", "wbn": "w2", "rnazprob": "0.97"},
        # Low-scoring and unmatched: not worth a curator's attention.
        {"locus_id": "L3", "locus_start": "600", "locus_end": "660", "wbn": "w3", "rnazprob": "0.1"},
    ]
    rows = scaffold_rows(domains, loci, "demo")
    by_id = {row["element_id"]: row for row in rows}

    assert "NOT REPORTED" not in by_id["TBD_01"]["notes"]
    assert "NOT REPORTED" in by_id["TBD_02"]["notes"], "a domain with no locus must be flagged"
    assert "TBD_extra_L2" in by_id, "a strong locus outside the reference needs a decision"
    assert "TBD_extra_L3" not in by_id
    # Every row is coordinate-complete but label-free.
    assert all(row["element_id"].startswith("TBD") for row in rows)
    assert all(row["start"] and row["end"] for row in rows)


def test_scaffold_does_not_double_report_sibling_loci() -> None:
    """The same element is reported once per window length."""
    from rnaconsnake.tools.benchmark_scaffold import Domain, scaffold_rows

    domains = [Domain(start=10, end=60, n_pairs=12)]
    loci = [
        {"locus_id": "len100_0001", "locus_start": "8", "locus_end": "62", "wbn": "a", "rnazprob": "0.99"},
        {"locus_id": "len200_0001", "locus_start": "5", "locus_end": "70", "wbn": "b", "rnazprob": "0.98"},
    ]
    rows = scaffold_rows(domains, loci, "demo")
    assert len(rows) == 1
    assert not any("extra" in row["element_id"] for row in rows)


def test_scaffold_output_is_readable_by_the_benchmark(tmp_path: Path) -> None:
    """The scaffold must be schema-compatible, and refused until curated."""
    from rnaconsnake.tools.benchmark import BenchmarkError, evaluate, read_truth
    from rnaconsnake.tools.benchmark_scaffold import Domain, scaffold_rows, write_scaffold

    path = tmp_path / "scaffold.tsv"
    write_scaffold(scaffold_rows([Domain(10, 60, 12)], [], "demo"), path, "demo")

    truth = read_truth(path)
    assert len(truth) == 1
    # Coordinates are filled in, so it is the *labels* that block it, and the
    # benchmark still refuses it until a curator has been through.
    assert truth[0].curated is True
    assert truth[0].element_id.startswith("TBD")
    with pytest.raises(BenchmarkError, match="placeholder"):
        evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=False)


def test_shipped_jevg_truth_file_is_fully_curated() -> None:
    """The JEV-group truth file is curated repository data, not a scaffold."""
    from rnaconsnake.tools.benchmark import BenchmarkError, evaluate, read_truth

    truth = read_truth(Path("resources/benchmark/jevg_3utr_elements.tsv"), "jevg_3utr")
    assert len(truth) == 9
    assert all(element.curated for element in truth)

    names = {element.element_id for element in truth}
    assert {"xrRNA1", "xrRNA2", "DB1", "DB2", "sHP", "3SL"} <= names
    # HP1 and HP2 were not previously annotated; they are part of the finding.
    assert {"HP1", "HP2"} <= names

    # Spans are ordered, non-overlapping and within the alignment.
    spans = sorted((e.start, e.end) for e in truth)
    assert all(a[1] < b[0] for a, b in itertools.pairwise(spans))
    assert spans[0][0] >= 1
    assert spans[-1][1] <= 711

    # No placeholders survive, so the benchmark accepts it.
    try:
        evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=False)
    except BenchmarkError as error:  # pragma: no cover - guards a regression
        raise AssertionError(f"curated truth file was refused: {error}") from error


def test_flavivirus_truth_file_remains_an_uncurated_stub() -> None:
    """The broad MBFV file is still a stub; it must not be mistaken for curated."""
    from rnaconsnake.tools.benchmark import read_truth

    truth = read_truth(Path("resources/benchmark/flavivirus_elements.tsv"))
    assert all(not element.curated for element in truth)


def test_null_baseline_summarises_chance_recovery(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import null_baseline

    truth = _two_elements(tmp_path)
    arms = [
        [_locus("n1", 90, 210)],  # covers e1 only
        [_locus("n1", 90, 210), _locus("n2", 390, 510)],  # covers both
        [],  # covers neither
    ]
    baseline = null_baseline(truth, arms, min_overlap_fraction=0.5)
    assert baseline["arms"] == 3
    assert baseline["min"] == 0
    assert baseline["max"] == 2
    assert baseline["mean"] == pytest.approx(1.0)
    assert null_baseline(truth, [], 0.5) == {}


def test_recovery_report_warns_when_null_matches_the_real_arm(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, write_recovery

    truth = _two_elements(tmp_path)
    results = evaluate(truth, [_locus("r1", 90, 210)], 0.5, allow_uncurated=True)
    output = tmp_path / "recovery.tsv"
    # Null arms recover more than the real arm did.
    write_recovery(results, output, {"arms": 5, "mean": 2.0, "min": 2, "max": 2})
    text = read_text(output)
    assert "# null_baseline_recovered\tmean 2.0" in text
    assert "recovery_margin_over_null\t-1.0" in text
    assert "not evidence of detection" in text


def test_recovery_report_warns_when_the_overlap_test_is_vacuous(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, write_recovery

    truth = _two_elements(tmp_path)
    results = evaluate(truth, [_locus("r1", 90, 210), _locus("r2", 390, 510)], 0.5, allow_uncurated=True)
    output = tmp_path / "recovery.tsv"
    write_recovery(results, output, {"arms": 10, "mean": 1.8, "min": 1, "max": 2})
    text = read_text(output)
    assert "recovery_margin_over_null\t+0.2" in text
    assert "close to vacuous" in text


def test_recovery_report_is_quiet_when_the_margin_is_real(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, write_recovery

    truth = _two_elements(tmp_path)
    results = evaluate(truth, [_locus("r1", 90, 210), _locus("r2", 390, 510)], 0.5, allow_uncurated=True)
    output = tmp_path / "recovery.tsv"
    write_recovery(results, output, {"arms": 10, "mean": 0.4, "min": 0, "max": 1})
    text = read_text(output)
    assert "recovery_margin_over_null\t+1.6" in text
    assert "WARNING" not in text


def test_one_null_loci_group_is_one_arm(tmp_path: Path) -> None:
    """An arm's window lengths pool into a single baseline sample.

    Counting each window length as its own arm would dilute the mean with
    partial views of the same arm, and understate the baseline.
    """
    truth_path = tmp_path / "truth.tsv"
    truth_path.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "e1\txrRNA\taln\t100\t199\tfirst\n"
        "e2\tDB\taln\t400\t499\tsecond\n",
        encoding="utf-8",
    )
    qvalues = tmp_path / "qvalues.tsv"
    qvalues.write_text(
        "locus_id\twlen\tstart\tend\trnazprob\talifoldzscore\nlen100_0001\t100\t90\t210\t0.99\t-3.0\n",
        encoding="utf-8",
    )
    header = "locus_id,locus_start,locus_end,n_windows,members,wbn,rnazprob,alifoldzscore\n"
    len100 = tmp_path / "len100.nr.csv"
    len100.write_text(header + "len100_0001,90,210,1,a,a,0.99,-3.0\n", encoding="utf-8")
    len200 = tmp_path / "len200.nr.csv"
    len200.write_text(header + "len200_0001,390,510,1,b,b,0.99,-3.0\n", encoding="utf-8")

    output = tmp_path / "recovery.tsv"
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.benchmark",
            "--truth",
            str(truth_path),
            "--qvalues",
            str(qvalues),
            "--output",
            str(output),
            "--alignment",
            "aln",
            "--null-loci",
            str(len100),
            str(len200),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    text = read_text(output)
    # One arm, recovering both elements -- not two arms recovering one each.
    assert "over 1 null arms" in text
    assert "# null_baseline_recovered\tmean 2.0" in text


def test_benchmark_rule_passes_every_null_arm_as_one_group() -> None:
    """The recovery table is uninterpretable without the null baseline, so the
    rule that builds it has to supply the arms' locus tables itself."""
    text = read_text(Path("snakefile"))
    rule = text.split("rule benchmark_recovery:", 1)[1].split("\nrule ", 1)[0]
    assert "null_loci=" in rule
    assert 'cmd += ["--null-loci", *params.null_loci_by_arm[arm]]' in rule
    assert "for arm in sorted(params.null_loci_by_arm)" in rule


def test_reciprocal_overlap_exposes_an_oversized_locus(tmp_path: Path) -> None:
    """A locus far larger than the element scores 1.0 on overlap_fraction.

    That is how a screen reporting one huge locus can look perfectly sensitive.
    """
    from rnaconsnake.tools.benchmark import RECOVERY_COLUMNS, evaluate

    assert "reciprocal_overlap" in RECOVERY_COLUMNS
    truth = _two_elements(tmp_path)
    # One locus spanning the whole alignment contains both elements entirely.
    results = evaluate(truth, [_locus("huge", 1, 1000)], 0.5, allow_uncurated=True)
    for row in results:
        assert row["recovered"] == "yes"
        assert row["overlap_fraction"] == "1.0000"
        # ...but the element accounts for only a tenth of the locus.
        assert float(row["reciprocal_overlap"]) == pytest.approx(0.1, abs=1e-3)


def test_every_curated_truth_file_is_packaged() -> None:
    """`benchmark_truth` names a file in resources/benchmark by name, and an
    installed run resolves it inside the package. A truth file left out of the
    build hook exists in the repository and nowhere else."""
    import importlib.util
    import types

    # setup.py is read for its copy list alone, so setuptools is stubbed rather
    # than required: the test env does not build the package.
    setuptools = types.ModuleType("setuptools")
    setuptools.setup = lambda **kwargs: None
    command = types.ModuleType("setuptools.command")
    build_py = types.ModuleType("setuptools.command.build_py")
    build_py.build_py = type("build_py", (), {"run": lambda self: None})
    stubs = {
        "setuptools": setuptools,
        "setuptools.command": command,
        "setuptools.command.build_py": build_py,
    }
    saved = {name: sys.modules.get(name) for name in stubs}
    sys.modules.update(stubs)
    try:
        spec = importlib.util.spec_from_file_location("_rnaconsnake_setup", Path("setup.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous

    packaged = {source.name for source in module.WORKFLOW_SOURCES}
    for truth in sorted(Path("resources/benchmark").glob("*.tsv")):
        assert truth.name in packaged, f"{truth.name} is not copied into the package"


def test_shipped_denvg_truth_file_is_the_held_out_test_set() -> None:
    """DENVG was curated after the configuration was locked on JEVG.

    Four elements only: the scaffold domains that turned out not to be xrRNAs
    are excluded rather than labelled uncertainly, because an unbelieved row
    would inflate recovery.
    """
    from rnaconsnake.tools.benchmark import read_truth

    truth = read_truth(Path("resources/benchmark/denvg_3utr_elements.tsv"), "denvg_3utr")
    assert len(truth) == 4
    assert {e.element_id for e in truth} == {"DB1", "DB2", "sHP", "3SL"}
    assert all(e.curated for e in truth)

    spans = sorted((e.start, e.end) for e in truth)
    assert all(a[1] < b[0] for a, b in itertools.pairwise(spans))
    assert spans[-1][1] <= 488

    text = read_text(Path("resources/benchmark/denvg_3utr_elements.tsv"))
    # The exclusions must stay documented, not silently dropped.
    assert "EXCLUDED" in text
    assert "not xrRNAs" in text
    assert "HELD-OUT" in text


def test_scaffold_flags_domains_that_need_attention(capsys) -> None:
    """Wide domains bracket several elements; narrow ones are easy to miss.

    Lowering max_width was tested and does not improve which elements are
    captured -- a wide domain already contains them -- so the scaffold warns
    rather than changing the decomposition.
    """
    from rnaconsnake.tools.benchmark_scaffold import DEFAULT_MIN_WIDTH, Domain, scaffold_rows

    # A confirmed 14 nt element was dropped by the previous min-width of 15.
    assert DEFAULT_MIN_WIDTH <= 14

    scaffold_rows([Domain(1, 131, 30), Domain(200, 214, 4)], [], "demo")
    # scaffold_rows itself is silent; the warnings live in main(). Assert the
    # thresholds those warnings use are the ones documented.
    import inspect

    from rnaconsnake.tools import benchmark_scaffold

    source = inspect.getsource(benchmark_scaffold.main)
    assert "exceed 100 nt" in source
    assert "under 20 nt" in source


def test_benchmark_scaffold_cli_drafts_a_truth_file(tmp_path: Path) -> None:
    alignment = tmp_path / "aln.stk"
    alignment.write_text(
        "# STOCKHOLM 1.0\n"
        "#=GF ID scaffold\n"
        "a GGGGAAAACCCCAAAAAAAAGGGGGAAAAACCCCC\n"
        "b GGGGAAAACCCCAAAAAAAAGGGGGAAAAACCCCC\n"
        "#=GC SS_cons ((((....))))........(((((.....)))))\n"
        "//\n",
        encoding="utf-8",
    )
    output = tmp_path / "scaffold.tsv"

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.benchmark_scaffold",
            "--alignment",
            str(alignment),
            "--alignment-id",
            "scaffold",
            "--output",
            str(output),
            "--min-width",
            "8",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 0, result.stderr
    assert "consensus domains" in result.stdout
    # The labels are placeholders on purpose; saying so is part of the contract.
    assert "Labels are placeholders" in result.stdout
    text = read_text(output)
    assert "element_id\telement_class\talignment\tstart\tend\tnotes" in text
    rows = [line for line in text.splitlines() if line and not line.startswith(("#", "element_id"))]
    assert len(rows) == 2, text
    for row in rows:
        assert row.split("\t")[2] == "scaffold"

    # ...and the benchmark refuses the draft until a curator has been through:
    # the coordinates are filled in, so it is the labels that block it.
    from rnaconsnake.tools.benchmark import BenchmarkError, evaluate, read_truth

    truth = read_truth(output, "scaffold")
    assert all(element.curated for element in truth)
    with pytest.raises(BenchmarkError, match="placeholder"):
        evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=False)
