"""The null-model arm: replicate generation, seeding, empirical
FDR / q-values, and the threshold sweep over them."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from tests.helpers import (
    PYTHON,
    TOY_ALIGNMENT,
    _candidate_rows,
    _sweep_arm_inputs,
    _sweep_base,
    _write_summary_csv,
    read_text,
    subprocess_env,
    write_fake_randomize_aln,
)


def test_null_settings_default_to_disabled_and_legacy_paths() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    disabled = NullSettings.from_config({})
    assert disabled.method == "none"
    assert disabled.enabled is False
    assert disabled.arms() == []

    enabled = NullSettings.from_config({"null": {"method": "sissiz", "replicates": 3}})
    assert enabled.enabled is True
    assert enabled.arms() == ["real", "null_000", "null_001", "null_002"]

    # replicates: 0 must be as inert as method: none
    assert NullSettings.from_config({"null": {"method": "sissiz", "replicates": 0}}).enabled is False


def test_null_settings_accepts_unquoted_yaml_null_key() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    # An unquoted "null:" key in YAML parses as the null scalar, so the section
    # arrives under None. It must not be silently ignored.
    settings = NullSettings.from_config({None: {"method": "sissiz", "replicates": 2}})
    assert settings.method == "sissiz"
    assert settings.replicates == 2


def test_null_settings_rejects_unknown_method() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    with pytest.raises(ValueError, match=r"Unknown null\.method"):
        NullSettings.from_config({"null": {"method": "dinucleotide-magic"}})


def test_arm_seed_is_deterministic_arm_specific_and_none_for_real() -> None:
    from rnaconsnake.workflow_helpers import arm_seed

    assert arm_seed("real", 20261101) is None
    assert arm_seed("null_000", 20261101) == arm_seed("null_000", 20261101)
    assert arm_seed("null_000", 20261101) != arm_seed("null_001", 20261101)
    assert arm_seed("null_000", 20261101) != arm_seed("null_000", 20261102)


def test_null_model_make_arm_copies_real_alignment_verbatim(tmp_path: Path) -> None:
    from rnaconsnake.tools.null_model import make_arm_alignment

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")
    output = tmp_path / "arms" / "real" / "input_alignment.stk"
    make_arm_alignment("real", source, output)

    # A copy, never a symlink: downstream tools rewrite alignments in place.
    assert not output.is_symlink()
    assert read_text(output) == TOY_ALIGNMENT


def test_null_model_pool_roundtrip_with_fake_backend(tmp_path: Path) -> None:
    from rnaconsnake.tools.alignment_io import read_stockholm_alignment
    from rnaconsnake.tools.null_model import make_arm_alignment, simulate_pool

    if shutil.which("perl") is None:  # pragma: no cover - perl is present on CI
        pytest.skip("perl is required for the rnazRandomizeAln backend")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    randomize = write_fake_randomize_aln(bin_dir)

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")

    metadata = simulate_pool(
        source_path=source,
        output_path=tmp_path / "pool.stk",
        metadata_path=tmp_path / "pool.json",
        method="rnazRandomizeAln",
        replicates=2,
        seed=20261101,
        sissiz_command=["SISSIz"],
        randomize_command=[str(randomize)],
        workdir=tmp_path / "work",
    )
    assert metadata["seeded"] is True
    assert metadata["determinism"] == "seeded"
    assert len(metadata["replicate_diagnostics"]) == 2

    original = read_stockholm_alignment(source)
    for arm in ["null_000", "null_001"]:
        out = tmp_path / arm / "input_alignment.stk"
        make_arm_alignment(arm, source, out, pool_path=tmp_path / "pool.stk")
        replicate = read_stockholm_alignment(out)
        assert replicate.order == original.order
        assert replicate.length == original.length
        # Column shuffling preserves per-sequence gap counts exactly.
        for name in original.order:
            assert replicate.seqs[name].count("-") == original.seqs[name].count("-")

    first = read_text(tmp_path / "null_000" / "input_alignment.stk")
    second = read_text(tmp_path / "null_001" / "input_alignment.stk")
    assert first != second


def test_null_model_pool_is_reproducible_for_the_same_seed(tmp_path: Path) -> None:
    from rnaconsnake.tools.null_model import simulate_pool

    if shutil.which("perl") is None:  # pragma: no cover
        pytest.skip("perl is required for the rnazRandomizeAln backend")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    randomize = write_fake_randomize_aln(bin_dir)
    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")

    def run(tag: str, seed: int) -> str:
        simulate_pool(
            source_path=source,
            output_path=tmp_path / f"pool_{tag}.stk",
            metadata_path=tmp_path / f"pool_{tag}.json",
            method="rnazRandomizeAln",
            replicates=2,
            seed=seed,
            sissiz_command=["SISSIz"],
            randomize_command=[str(randomize)],
            workdir=tmp_path / f"work_{tag}",
        )
        return read_text(tmp_path / f"pool_{tag}.stk")

    assert run("a", 20261101) == run("b", 20261101)
    assert run("a", 20261101) != run("c", 20261102)


def test_null_model_rejects_single_sequence_alignment(tmp_path: Path) -> None:
    from rnaconsnake.tools.null_model import NullModelError, simulate_pool

    source = tmp_path / "input.stk"
    source.write_text("# STOCKHOLM 1.0\nseqA ACGU\n//\n", encoding="utf-8")
    with pytest.raises(NullModelError, match="at least 2 sequences"):
        simulate_pool(
            source_path=source,
            output_path=tmp_path / "pool.stk",
            metadata_path=tmp_path / "pool.json",
            method="rnazRandomizeAln",
            replicates=1,
            seed=1,
            sissiz_command=["SISSIz"],
            randomize_command=["rnazRandomizeAln.pl"],
        )


def test_empirical_fdr_and_qvalue_envelope_are_monotone() -> None:
    from rnaconsnake.tools.calibration import empirical_fdr, qvalue_envelope

    real = [0.95, 0.90, 0.80, 0.60]
    nulls = [[0.60, 0.55], [0.62, 0.50]]
    fdr = empirical_fdr(real, nulls)
    envelope = qvalue_envelope(fdr)

    assert set(fdr) == {0.60, 0.80, 0.90, 0.95}
    # Nothing in the null arms reaches 0.80, so the high-score end is clean.
    assert fdr[0.95] == 0.0
    assert fdr[0.80] == 0.0
    assert fdr[0.60] == pytest.approx(1.0 / 4)

    values = [envelope[threshold] for threshold in sorted(envelope)]
    assert values == sorted(values, reverse=True) or all(
        values[i] >= values[i + 1] for i in range(len(values) - 1)
    )
    assert all(0.0 <= value <= 1.0 for value in values)


def test_calibration_summary_records_every_clustering_parameter(tmp_path: Path) -> None:
    """Clustering decides how many loci each arm reports, and the q-values are
    counts over loci. A summary missing those parameters cannot be reproduced
    from -- and the export manifest copies this block verbatim."""
    from dataclasses import fields

    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "2"))},
    }
    for index in range(2):
        arm = f"null_{index:03d}"
        arm_inputs[arm] = {100: _write_summary_csv(tmp_path / f"{arm}.csv", _candidate_rows(0.20, 0.5, "0"))}

    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "rnazRandomizeAln", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )

    recorded = summary["thresholds"]
    for field in fields(Thresholds):
        assert field.name in recorded, f"{field.name} is not recorded in summary.json"
    assert recorded["max_container_width"] == 120
    assert recorded["container_min_coverage"] == pytest.approx(0.8)
    assert recorded["representative_rule"] == "widest"


def test_calibrate_writes_funnel_qvalues_and_summary(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "2"))},
    }
    for index in range(3):
        arm = f"null_{index:03d}"
        arm_inputs[arm] = {100: _write_summary_csv(tmp_path / f"{arm}.csv", _candidate_rows(0.20, 0.5, "0"))}

    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "rnazRandomizeAln", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )

    assert summary["counting_unit"] == "merged_loci"
    assert summary["fdr_conditional_on_stage_one"] is True
    assert summary["cascade_fdr"] == 0.0
    assert summary["counts"]["real_loci"] == 4
    assert summary["q_resolution"] == pytest.approx(1 / 3)
    assert summary["warnings"] == []

    funnel = read_text(tmp_path / "calibration" / "funnel.tsv")
    assert "# fdr_conditional_on_stage_one\ttrue" in funnel
    assert "counts are on merged loci" in funnel
    # One row per filter stage per arm, plus the null_mean aggregate.
    for arm in ["real", "null_000", "null_001", "null_002", "null_mean"]:
        for stage in ["windows", "loci", "rnaz", "alifoldz", "rscape", "cascade"]:
            assert f"{arm}\t100\t{stage}\t" in funnel

    qvalues = read_text(tmp_path / "calibration" / "qvalues.tsv")
    assert qvalues.count("\tyes\t") == 4
    assert "q_cascade" in qvalues

    dists = read_text(tmp_path / "calibration" / "score_dists.tsv")
    assert "null_000\t100\t" in dists
    assert "real\t" not in dists.split("\n", 1)[1]


def test_calibrate_is_reproducible_for_identical_inputs(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "2"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.30, -0.5, "0"))},
        "null_001": {100: _write_summary_csv(tmp_path / "n1.csv", _candidate_rows(0.95, -3.0, "1"))},
    }
    thresholds = Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2)
    for tag in ["a", "b"]:
        calibrate(
            arm_inputs=arm_inputs,
            thresholds=thresholds,
            null_metadata={"method": "sissiz", "seed": 20261101, "warnings": []},
            output_dir=tmp_path / tag,
            two_stage=True,
        )
    assert read_text(tmp_path / "a" / "qvalues.tsv") == read_text(tmp_path / "b" / "qvalues.tsv")
    assert read_text(tmp_path / "a" / "funnel.tsv") == read_text(tmp_path / "b" / "funnel.tsv")


def test_calibrate_warns_when_collapse_ratio_diverges_between_arms(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    # Real arm: 4 windows collapsing into 1 locus. Null arm: 4 windows, 4 loci.
    real_rows = [
        {
            "wbn": f"RC_100_0001_aln_{1 + index * 10}_{100 - index * 10}",
            "rnazprob": "0.95",
            "alifoldzscore": "-3.0",
            "rscape_covary_count": "1",
            "nrseq": "6",
            "alilen": "100",
        }
        for index in range(4)
    ]
    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", real_rows)},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.10, 1.0, "0"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=False,
    )
    assert summary["warnings"]
    assert "collapse ratio" in summary["warnings"][0]
    assert "# WARNING\t" in read_text(tmp_path / "calibration" / "funnel.tsv")


def test_calibrate_rejects_stage_one_looser_than_reported_threshold(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "1"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.1, 1.0, "0"))},
    }
    with pytest.raises(ValueError, match="stage1_rnaz_prob"):
        calibrate(
            arm_inputs=arm_inputs,
            thresholds=Thresholds(0.9, -2.0, 1, 0.95, 1, 0.2),
            null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
            output_dir=tmp_path / "calibration",
            two_stage=True,
        )


def test_calibrate_treats_stage_one_skips_as_failing_alifoldz(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    skipped = _candidate_rows(0.97, -4.0, "1")
    for row in skipped:
        row["alifoldzscore"] = "NA"
    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", skipped)},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.1, 1.0, "0"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )
    # A missing AlifoldZ score must never be read as a good (very negative) one.
    assert summary["counts"]["real_cascade_survivors"] == 0


def test_calibrate_drops_rscape_from_cascade_when_it_never_ran(tmp_path: Path) -> None:
    """do_rscape: false makes every count NA; that must not zero the headline."""
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "NA"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.10, 1.0, "NA"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )
    assert summary["rscape_evaluated"] is False
    assert summary["cascade_filters"] == ["rnaz", "alifoldz"]
    assert summary["counts"]["real_cascade_survivors"] == 4
    assert summary["cascade_fdr"] == 0.0
    assert any("R-scape produced no covariation counts" in w for w in summary["warnings"])
    assert "# rscape_in_cascade\tfalse" in read_text(tmp_path / "calibration" / "funnel.tsv")


def test_calibrate_keeps_rscape_in_cascade_when_it_ran(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "0"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.10, 1.0, "0"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )
    assert summary["rscape_evaluated"] is True
    assert summary["cascade_filters"] == ["rnaz", "alifoldz", "rscape"]
    # Zero covarying pairs is a real, reported result, not a missing value.
    assert summary["counts"]["real_cascade_survivors"] == 0
    assert "# rscape_in_cascade\ttrue" in read_text(tmp_path / "calibration" / "funnel.tsv")


def test_sissiz_backend_is_reported_as_unseeded() -> None:
    """SISSIz draws its seed from the clock, so null.seed cannot reproduce it.

    Two invocations inside the same second agree, which makes this easy to
    mis-measure; the metadata must not claim reproducibility it cannot deliver.
    """
    import inspect

    from rnaconsnake.tools import null_model

    source = inspect.getsource(null_model._simulate_sissiz)
    assert '"determinism": "unseeded"' in source
    assert '"seeded": False' in source
    assert "tool_is_deterministic" not in source


def test_unseeded_backend_warns_that_the_seed_does_not_reproduce_the_pool(
    tmp_path: Path,
) -> None:
    from rnaconsnake.tools.null_model import simulate_pool

    if shutil.which("perl") is None:  # pragma: no cover
        pytest.skip("perl is required for the rnazRandomizeAln backend")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    randomize = write_fake_randomize_aln(bin_dir)
    # Make the fake report itself as a non-Perl command so it takes the
    # unseeded path, exercising the warning.
    plain = bin_dir / "randomize_plain"
    plain.write_text(randomize.read_text(encoding="utf-8"), encoding="utf-8")
    plain.chmod(0o755)

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")
    metadata = simulate_pool(
        source_path=source,
        output_path=tmp_path / "pool.stk",
        metadata_path=tmp_path / "pool.json",
        method="rnazRandomizeAln",
        replicates=1,
        seed=1,
        sissiz_command=["SISSIz"],
        randomize_command=[str(plain)],
        workdir=tmp_path / "work",
    )
    assert metadata["seeded"] is False
    assert any("not seedable" in warning for warning in metadata["warnings"])


def test_adopt_pool_reuses_a_pinned_pool_and_validates_it(tmp_path: Path) -> None:
    """Pinning is how a SISSIz calibration is made reproducible."""
    from rnaconsnake.tools.null_model import NullModelError, adopt_pool

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")

    pool = tmp_path / "pinned.stk"
    pool.write_text(
        "\n".join(
            [
                "# STOCKHOLM 1.0",
                "#=GF ID null_000",
                "seqA GGCUAGCUAGCUAACGUAGCUAGCUAGGCAUCGAUCGAUCG",
                "seqB GGCUAGCUAG---ACGUAGCUAGCUAGGCAUCGAUCGAUCG",
                "seqC GGCUAGCUAGCUAACGUAGCUCGCUAGGCAUCGAUCG---G",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )

    metadata = adopt_pool(
        pool_path=pool,
        source_path=source,
        output_path=tmp_path / "pool.stk",
        metadata_path=tmp_path / "pool.json",
        replicates=1,
        method="sissiz",
    )
    assert metadata["determinism"] == "pinned"
    assert metadata["pinned_from"] == str(pool)
    assert read_text(tmp_path / "pool.stk") == read_text(pool)

    with pytest.raises(NullModelError, match="need 2"):
        adopt_pool(
            pool_path=pool,
            source_path=source,
            output_path=tmp_path / "pool2.stk",
            metadata_path=tmp_path / "pool2.json",
            replicates=2,
            method="sissiz",
        )


def test_adopt_pool_rejects_a_pool_from_a_different_alignment(tmp_path: Path) -> None:
    """A pinned pool from another alignment would silently miscalibrate."""
    from rnaconsnake.tools.null_model import NullModelError, adopt_pool

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")
    wrong = tmp_path / "wrong.stk"
    wrong.write_text(
        "# STOCKHOLM 1.0\n#=GF ID null_000\nseqA ACGU\nseqB ACGA\nseqC ACGC\n//\n",
        encoding="utf-8",
    )
    with pytest.raises(NullModelError, match="length"):
        adopt_pool(
            pool_path=wrong,
            source_path=source,
            output_path=tmp_path / "pool.stk",
            metadata_path=tmp_path / "pool.json",
            replicates=1,
            method="sissiz",
        )


def test_perl_seeded_command_wraps_scripts_and_falls_back(tmp_path: Path) -> None:
    from rnaconsnake.workflow_helpers import perl_seed_env, perl_seeded_command

    script = tmp_path / "helper.pl"
    script.write_text("print 1;\n", encoding="utf-8")
    command, seeded = perl_seeded_command([str(script)], 42, ["-f"])
    assert seeded is True
    assert command[0] == "perl"
    assert command[3] == "42"
    assert command[-1] == "-f"

    # Non-Perl commands run unchanged, and say so rather than pretending.
    fallback, seeded = perl_seeded_command(["/bin/echo"], 42, ["hi"])
    assert seeded is False
    assert fallback == ["/bin/echo", "hi"]

    with pytest.raises(ValueError, match="Empty command"):
        perl_seeded_command([], 1)

    # Hash order must be pinned too; srand alone is not enough.
    env = perl_seed_env({})
    assert env["PERL_HASH_SEED"] == "0"
    assert env["PERL_PERTURB_KEYS"] == "0"


def test_derived_seed_is_deterministic_and_per_candidate() -> None:
    from rnaconsnake.workflow_helpers import derived_seed

    assert derived_seed(7, "cand_a") == derived_seed(7, "cand_a")
    assert derived_seed(7, "cand_a") != derived_seed(7, "cand_b")
    assert derived_seed(7, "cand_a") != derived_seed(8, "cand_a")


def test_snakefile_seeds_alifoldz_per_candidate() -> None:
    """alifoldz.pl shuffles internally and has no seed option.

    Its z-scores feed the q-values, so an unseeded AlifoldZ makes the whole
    calibration irreproducible no matter how the null pool was generated.
    """
    text = read_text(Path("snakefile"))
    assert 'ALIFOLDZ_SEED = config.get("alifoldz_seed")' in text
    assert "perl_seeded_command(" in text
    assert "derived_seed(int(ALIFOLDZ_SEED), wildcards.file)" in text


def test_config_seeds_alifoldz_by_default() -> None:
    import yaml as _yaml

    payload = _yaml.safe_load(read_text(Path("config.yaml")))
    assert payload.get("alifoldz_seed") is not None, "a calibrated run must be reproducible by default"


def test_null_settings_accepts_a_pinned_pool() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    settings = NullSettings.from_config({"null": {"method": "sissiz", "pool_file": "/tmp/pool.stk"}})
    assert settings.pool_file == "/tmp/pool.stk"
    # Snakemake stringifies nested --config values, so "None" must not be a path.
    for sentinel in ["None", "", "null", "~"]:
        assert (
            NullSettings.from_config({"null": {"method": "sissiz", "pool_file": sentinel}}).pool_file is None
        )


def test_threshold_sweep_honours_the_representative_rule(tmp_path: Path) -> None:
    """The sweep clusters with the rule it was given. It used to drop
    `base.representative_rule` on the floor and always use `widest`, so a sweep
    of a run configured otherwise described a different set of scores."""
    from rnaconsnake.tools.calibration import Thresholds
    from rnaconsnake.tools.threshold_sweep import sweep

    # One locus: a wide, weak window containing a narrow, strong one. "widest"
    # represents it by the weak parent; "best_rnaz" by the strong fragment.
    rows = [
        {
            "wbn": "RC_100_0001_aln_10_120",
            "rnazprob": "0.10",
            "alifoldzscore": "-0.1",
            "rscape_covary_count": "NA",
            "nrseq": "4",
            "alilen": "111",
            "maxcovarval": "0",
            "maxcovarcount": "0",
            "consensus_mfe": "-5",
            "sci": "0.5",
            "alifold_consstruc": "",
        },
        {
            "wbn": "RC_100_0001_aln_20_60",
            "rnazprob": "0.99",
            "alifoldzscore": "-5.0",
            "rscape_covary_count": "NA",
            "nrseq": "4",
            "alilen": "41",
            "maxcovarval": "0",
            "maxcovarcount": "0",
            "consensus_mfe": "-9",
            "sci": "0.9",
            "alifold_consstruc": "",
        },
    ]

    def survivors(rule: str) -> int:
        base = Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2, "containment", 0.9, 0, 0.8, rule)
        arm_inputs = {
            "real": {100: _write_summary_csv(tmp_path / f"real_{rule}.csv", rows)},
            "null_000": {100: _write_summary_csv(tmp_path / f"null_{rule}.csv", [])},
        }
        (point,) = sweep(arm_inputs, [0.9], [-2.0], base)
        return point.real_survivors

    assert survivors("widest") == 0
    assert survivors("best_rnaz") == 1


def test_threshold_sweep_reports_survivors_and_fdr(tmp_path: Path) -> None:
    from rnaconsnake.tools.threshold_sweep import sweep

    strong = _candidate_rows(0.97, -4.0, "NA")
    weak = _candidate_rows(0.20, 0.5, "NA")
    points = sweep(
        _sweep_arm_inputs(tmp_path, strong, [weak, weak]),
        rnaz_grid=[0.5, 0.9],
        alifoldz_grid=[-1.0, -3.0],
        base=_sweep_base(),
    )
    assert len(points) == 4
    by_key = {(p.rnaz_prob, p.alifoldz): p for p in points}

    # The real arm clears every combination; the null arm clears none.
    for point in points:
        assert point.real_survivors == 4
        assert point.null_mean == 0.0
        assert point.fdr == 0.0
    assert by_key[(0.9, -3.0)].null_sd == 0.0


def test_threshold_sweep_tightening_reduces_null_survivors(tmp_path: Path) -> None:
    from rnaconsnake.tools.threshold_sweep import sweep

    real = _candidate_rows(0.97, -4.0, "NA")
    # Null loci that clear a loose threshold but not a strict one.
    null = _candidate_rows(0.92, -1.2, "NA")
    points = sweep(
        _sweep_arm_inputs(tmp_path, real, [null, null]),
        rnaz_grid=[0.9],
        alifoldz_grid=[-1.0, -3.0],
        base=_sweep_base(),
    )
    loose = next(p for p in points if p.alifoldz == -1.0)
    strict = next(p for p in points if p.alifoldz == -3.0)
    assert loose.null_mean > strict.null_mean
    assert loose.fdr > strict.fdr


def test_threshold_sweep_counts_reference_recovery(tmp_path: Path) -> None:
    """Recovery is sensitivity only, and must be reported alongside FDR."""
    from rnaconsnake.tools.threshold_sweep import sweep

    real = _candidate_rows(0.97, -4.0, "NA")  # spans 1-100, 501-600, 1001-1100, 1501-1600
    null = _candidate_rows(0.10, 1.0, "NA")
    spans = [(10, 90), (510, 590), (9000, 9100)]  # third is never covered
    points = sweep(
        _sweep_arm_inputs(tmp_path, real, [null]),
        rnaz_grid=[0.9],
        alifoldz_grid=[-2.0],
        base=_sweep_base(),
        reference_spans=spans,
    )
    assert points[0].recovered == 2
    assert points[0].reference_total == 3


def test_threshold_sweep_stage_one_gate_follows_the_reported_threshold(tmp_path: Path) -> None:
    """stage1 must never exceed the reported RNAz cutoff, or the sweep would
    select candidates whose AlifoldZ was never computed."""
    from rnaconsnake.tools.threshold_sweep import sweep

    real = _candidate_rows(0.6, -4.0, "NA")
    points = sweep(
        _sweep_arm_inputs(tmp_path, real, [_candidate_rows(0.1, 1.0, "NA")]),
        rnaz_grid=[0.2],
        alifoldz_grid=[-2.0],
        base=_sweep_base(),  # stage1 default is 0.5, above the 0.2 grid point
    )
    # No exception, and the loose threshold is honoured.
    assert points[0].real_survivors == 4


def test_write_sweep_marks_undefined_fdr(tmp_path: Path) -> None:
    from rnaconsnake.tools.threshold_sweep import SweepPoint, write_sweep

    output = tmp_path / "sweep.tsv"
    write_sweep(
        [SweepPoint(0.9, -2.0, 0, 0.0, 0.0, None, None, None)],
        output,
    )
    text = read_text(output)
    # No real survivors means no FDR, and that must not read as zero.
    assert "\tNA\tNA\tNA" in text
    assert "sensitivity only" in text


def test_threshold_sweep_cli_writes_a_grid(tmp_path: Path) -> None:
    real = _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "2"))
    null = _write_summary_csv(tmp_path / "null.csv", _candidate_rows(0.20, 0.5, "0"))
    reference = tmp_path / "truth.tsv"
    reference.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\ne1\txrRNA\taln\t1\t100\tfirst\n",
        encoding="utf-8",
    )
    output = tmp_path / "sweep.tsv"

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.threshold_sweep",
            "--arm-input",
            f"real:100:{real}",
            "--arm-input",
            f"null_000:100:{null}",
            "--output",
            str(output),
            "--rnaz-grid",
            "0.5,0.9",
            "--alifoldz-grid=-1,-2",
            "--reference",
            str(reference),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 0, result.stderr
    assert "4 threshold combinations" in result.stdout
    lines = read_text(output).splitlines()
    header = next(line for line in lines if not line.startswith("#"))
    assert "rnaz_prob" in header
    assert "alifoldz" in header
    # One row per grid point, plus the recovery column the reference adds.
    assert len([line for line in lines if line.startswith("0.")]) == 4
    assert "recovered" in header


def test_threshold_sweep_reads_reference_spans_and_ignores_labels(tmp_path: Path) -> None:
    from rnaconsnake.tools.threshold_sweep import read_reference_spans

    path = tmp_path / "ref.tsv"
    path.write_text(
        "# a comment line\n"
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "e1\txrRNA\taln\t10\t20\tfine\n"
        "e2\tDB\taln\tTBD\tTBD\tuncurated\n",
        encoding="utf-8",
    )
    # The uncurated row has no coordinates, so it contributes no span.
    assert read_reference_spans(path) == [(10, 20)]


def test_configuration_is_locked_to_the_recorded_values() -> None:
    """The screen configuration was fixed on the JEV group before the dengue
    group was evaluated. Changing these values silently would forfeit the only
    out-of-sample number the study has."""
    import yaml as _yaml

    payload = _yaml.safe_load(read_text(Path("config.yaml")))
    assert payload["dereplicate"]["representative"] == "widest"
    assert payload["calibration"]["rnaz_prob_threshold"] == 0.8
    assert payload["calibration"]["alifoldz_threshold"] == -1.5
