"""Repository-level invariants: the scope boundary, packaging,
documentation, and the container build context."""

from __future__ import annotations

import subprocess
from pathlib import Path

from tests.helpers import (
    ROOT,
    WEB_FORBIDDEN_PATHS,
    WEB_FORBIDDEN_PATTERNS,
    _tracked_files,
    read_text,
)


def test_docs_reference_supported_input_formats() -> None:
    for doc in [
        Path("README.md"),
        Path("docs/usage.md"),
        Path("docs/pipeline_summary.md"),
    ]:
        text = read_text(doc)
        assert ".stk" in text
        if doc.name == "usage.md":
            assert "{stk,aln}" in text
        else:
            assert ".aln" in text


def test_dereplication_doc_matches_the_implementation() -> None:
    """The spec in docs/dereplication.md must not drift from the code."""
    from rnaconsnake.tools.dereplicate import METHODS, NR_COLUMNS
    from rnaconsnake.workflow_helpers import CalibrationSettings

    doc = read_text(Path("docs/dereplication.md"))
    defaults = CalibrationSettings.from_config({})

    # Every method is documented, and the documented default is the real one.
    for method in METHODS:
        assert f"`{method}`" in doc, f"method {method} is undocumented"
    assert f"`{defaults.dereplicate_method}` is the default" in doc

    # Documented column order matches what is written.
    assert "locus_id  locus_start  locus_end  n_windows  members" in doc
    assert NR_COLUMNS[:5] == ["locus_id", "locus_start", "locus_end", "n_windows", "members"]

    # Documented default thresholds match the shipped ones.
    assert f"default {defaults.pair_containment}" in doc
    assert f"default {defaults.locus_min_overlap}" in doc

    # The doc claims the step is unconditional, so the target must be among the
    # default targets rule all requests.
    snakefile = read_text(Path("snakefile"))
    targets = snakefile.split("def pipeline_targets():")[1].split("\nrule ")[0]
    assert "RNAConSnake.nr.csv" in targets
    assert "pipeline_targets()" in snakefile.split("rule all:")[1].split("\nrule ")[0]

    # And that the standalone invocation it prints actually exists.
    assert "python -m rnaconsnake.tools.dereplicate" in doc
    for flag in ["--input", "--output", "--method", "--label"]:
        assert flag in doc


def test_dereplication_doc_is_linked_from_the_user_docs() -> None:
    for path, needle in [
        (Path("README.md"), "docs/dereplication.md"),
        (Path("docs/usage.md"), "dereplication.md"),
        (Path("docs/pipeline_summary.md"), "dereplication.md"),
    ]:
        assert needle in read_text(path), f"{path} does not link the de-replication spec"


def test_no_web_assets_are_tracked() -> None:
    offenders = [path for path in _tracked_files() if path.suffix.lower() in WEB_FORBIDDEN_PATHS]
    assert not offenders, "web assets are out of scope for RNAConSnake: " + ", ".join(
        str(path) for path in offenders
    )


def test_no_rendering_code_or_dependency_is_tracked() -> None:
    import re as _re

    import tests.helpers

    pattern = _re.compile("|".join(WEB_FORBIDDEN_PATTERNS), _re.IGNORECASE)
    # The guard names the forbidden things, so it cannot police the file that
    # spells them out. That is the module defining the patterns, not
    # necessarily the one running the test -- derived from the module so the
    # exclusion follows the literals if they move again, and stays exactly two
    # files wide. Everything else in the repository is still checked.
    exempt = {Path(__file__).resolve(), Path(tests.helpers.__file__).resolve()}
    offenders: list[str] = []
    for path in _tracked_files():
        full = ROOT / path
        if not full.is_file():
            continue
        if full.resolve() in exempt:
            continue
        try:
            text = full.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for match in pattern.finditer(text):
            offenders.append(f"{path}: {match.group(0)!r}")
    assert not offenders, "presentation-layer code leaked into RNAConSnake:\n" + "\n".join(offenders)


def test_the_declared_python_floor_is_one_the_dependencies_support() -> None:
    """A floor below what snakemake supports is a promise nothing can keep:
    `pip install` fails outright on such an interpreter, and CI would test a
    version the package can never run on."""
    import re as _re
    from importlib.metadata import metadata

    def floor(spec: str) -> tuple[int, ...]:
        match = _re.search(r">=\s*(\d+)\.(\d+)", spec or "")
        assert match, f"no lower bound in {spec!r}"
        return (int(match.group(1)), int(match.group(2)))

    declared = floor(
        _re.search(r'requires-python\s*=\s*"([^"]+)"', read_text(Path("pyproject.toml"))).group(1)
    )
    required = floor(metadata("snakemake").get("Requires-Python", ""))
    assert declared >= required, (
        f"pyproject claims Python {declared[0]}.{declared[1]}, but snakemake needs "
        f"{required[0]}.{required[1]}"
    )

    # ...and CI must not test below it either.
    workflow = read_text(Path(".github/workflows/ci.yml"))
    versions = _re.search(r"python-version: \[([^\]]+)\]", workflow).group(1)
    tested = [tuple(int(part) for part in v.strip().strip('"').split(".")) for v in versions.split(",")]
    assert min(tested) >= declared, f"CI tests {min(tested)}, below the declared floor {declared}"


def test_packaging_declares_no_web_dependencies() -> None:
    text = read_text(Path("pyproject.toml")).lower()
    for forbidden in ["jinja", "flask", "fastapi", "markdown", "starlette"]:
        assert forbidden not in text, f"{forbidden} is a presentation concern, not a pipeline one"


def test_workflow_emits_no_html() -> None:
    """No browser-style HTML reporting in this project."""
    snakefile = read_text(Path("snakefile"))
    assert ".html" not in snakefile
    for path in _tracked_files():
        assert not str(path).endswith(".html")


def test_readme_documents_the_tools_conda_cannot_supply() -> None:
    """Three tools are on no package index; a user who does not know that is
    stuck before the first run."""
    readme = read_text(Path("README.md"))
    container = read_text(Path("container/README.md"))

    assert "https://github.com/mtw/SISSIz" in readme
    assert "https://github.com/mtw/SISSIz" in container
    for tool in ["SISSIz", "alifoldz.pl"]:
        assert tool in readme, f"{tool} is not mentioned in the README"
    # The --use-conda limitation is a stated submission blocker; it must not
    # silently disappear from the docs.
    assert "--use-conda" in readme

    # prepare-context.sh must name the sources in its failure path, since that
    # is where a user actually hits the problem.
    prepare = read_text(Path("container/prepare-context.sh"))
    assert "github.com/mtw/SISSIz" in prepare
    assert "RNAz source tarball" in prepare


def test_readme_covers_the_features_added_this_cycle() -> None:
    readme = read_text(Path("README.md"))
    for topic in [
        "null-model calibration",  # the calibration arm
        "De-replicating",  # locus de-replication
        "screenability",  # blind-region reporting
        "Positive control",  # benchmark
        "container",  # containerised toolchain
    ]:
        assert topic.lower() in readme.lower(), f"README does not mention {topic}"


def test_readme_links_resolve() -> None:
    """Links must point at files the repository actually ships.

    Checking only that a path exists locally is not enough: a file present in
    the working tree but untracked (or gitignored) resolves here and 404s for
    anyone who clones.
    """
    import re as _re

    readme = read_text(Path("README.md"))
    targets = [t for t in _re.findall(r"\[`[^`]+`\]\(([^)]+)\)", readme) if not t.startswith("http")]
    for target in targets:
        assert Path(target).exists(), f"README links to a missing path: {target}"

    # An unstaged file is merely not added yet; an *ignored* one will never be
    # in a clone, so the link is permanently broken for everyone else.
    ignored = subprocess.run(
        ["git", "check-ignore", *targets],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    ).stdout.split()
    assert not ignored, "README links to gitignored paths, which would 404 in a clone: " + ", ".join(ignored)


def test_public_container_dir_has_no_lab_tooling() -> None:
    """Host-specific job-distribution scripts are lab infrastructure and do not
    belong in a public repository."""
    names = {p.name for p in Path("container").iterdir() if p.is_file()}
    # Third-party patches live in their own directory and are scanned below.
    assert (Path("container") / "patches").is_dir()
    assert names == {
        "Dockerfile",
        "environment.container.yaml",
        "prepare-context.sh",
        "README.md",
    }
    for path in Path("container").rglob("*"):
        if path.is_file():
            text = path.read_text(encoding="utf-8", errors="ignore")
            for leak in ["/Users/mtw", "venus", "mercury"]:
                assert leak not in text, f"{path.name} leaks '{leak}'"


def test_no_text_file_discloses_a_non_public_downstream_project() -> None:
    """The public repository must not describe what is built on top of it.

    Scope statements ("this project does not render HTML") are fine; naming a
    separate closed-source project, or characterising its business status, is
    not.
    """
    import re as _re

    disclosive = _re.compile(
        r"proprietary|closed[- ]source|commercial(ly)? licen|sibling project"
        r"|rnaconsnake[_-]portal",
        _re.IGNORECASE,
    )
    offenders: list[str] = []
    for path in _tracked_files():
        full = ROOT / path
        if not full.is_file():
            continue
        # This test names the forbidden words, so it cannot police itself.
        if full.resolve() == Path(__file__).resolve():
            continue
        try:
            text = full.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for match in disclosive.finditer(text):
            line = text[: match.start()].count("\n") + 1
            offenders.append(f"{path}:{line}: {match.group(0)!r}")
    assert not offenders, "public files disclose a non-public downstream project:\n" + "\n".join(offenders)
