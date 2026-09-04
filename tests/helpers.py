"""Shared fixtures, fake external tools, and builders for the test suite.

Split out of the single test module so each test file can state what it
uses. Nothing here asserts; it only builds inputs and runs tools.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

from rnaconsnake import cli

ROOT = Path(__file__).resolve().parents[1]


FIXTURES = ROOT / "tests" / "fixtures"


PYTHON = sys.executable


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _parse_cli_args(argv: list[str]) -> argparse.Namespace:
    """RNAcs's parsed arguments for one command line, without running anything."""
    original = sys.argv
    sys.argv = [cli.PROGRAM_NAME, *argv]
    try:
        return cli.parse_args()
    finally:
        sys.argv = original


def subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    pythonpath = str(ROOT / "src")
    env["PYTHONPATH"] = pythonpath + (f":{env['PYTHONPATH']}" if env.get("PYTHONPATH") else "")
    return env


def write_fake_rnalalifold(bin_dir: Path) -> None:
    script = bin_dir / "RNALalifold"
    script.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import sys
            from pathlib import Path

            args = sys.argv[1:]
            if "--version" in args:
                import RNA

                print("RNALalifold " + RNA.__version__)
                raise SystemExit(0)
            prefix = "RC"
            wlen = "100"
            input_format = ""
            for i, arg in enumerate(args):
                if arg == "--id-prefix" and i + 1 < len(args):
                    prefix = args[i + 1]
                if arg == "-L" and i + 1 < len(args):
                    wlen = args[i + 1]
                if arg == "-f" and i + 1 < len(args):
                    input_format = args[i + 1]

            out_path = Path.cwd() / f"{prefix}_0001.stk"
            out_path.write_text(
                "\\n".join(
                    [
                        "# STOCKHOLM 1.0",
                        f"#=GF ID {prefix}_0001_aln_1_12",
                        "seqA ACGUACGU----",
                        "seqB ACGUACGU----",
                        "#=GC SS_cons <<<<....>>>>",
                        "//",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            sys.stdout.write(f"fake RNALalifold completed for window {wlen} format {input_format}\\n")
            """
        ),
        encoding="utf-8",
    )
    script.chmod(0o755)


def write_fake_postprocess_tools(bin_dir: Path) -> None:
    tools = {
        "esl-reformat": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[-1])
print("CLUSTAL W")
print()
for line in inp.read_text(encoding="utf-8").splitlines():
    if not line or line.startswith("#") or line == "//":
        continue
    name, seq = line.split(None, 1)
    print(f"{name} {seq}")
print("***")
""",
        "RNAz": """\
#!/usr/bin/env python3
print("Mean z-score: -3.21")
print("Mean MFE: -14.10")
print("Structure conservation index: 0.58")
print("SVM RNA-class probability: 0.95")
""",
        "alifoldz.pl": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

# The real alifoldz.pl shells out to RNAalifold without --noPS, so an
# undeclared "alirna.ps" lands in whatever directory it was run from. Emit one
# too, so the workflow's confinement of that stray file is actually exercised.
(Path.cwd() / "alirna.ps").write_text("%!PS stray from alifoldz\\n", encoding="utf-8")
_ = sys.stdin.read()
print("#           Input: 3 sequences of 41 columns")
print("fake alifoldz header")
print("-3.21")
""",
        "RNAalifold": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    # The real binary prints and exits; it does not fold, and must not write.
    # The version is the module's: RNAcs requires one ViennaRNA build.
    import RNA
    print("RNAalifold " + RNA.__version__)
    raise SystemExit(0)
prefix = "fake"
for i, arg in enumerate(args):
    if arg == "--id-prefix" and i + 1 < len(args):
        prefix = args[i + 1]

cwd = Path.cwd()
(cwd / "RNAalifold_results.stk").write_text(
    "\\n".join(
        [
            "# STOCKHOLM 1.0",
            f"#=GF ID {prefix}",
            "seqA ACGUACGU----",
            "seqB ACGUACGU----",
            "#=GC SS_cons <<<<....>>>>",
            "//",
            "",
        ]
    ),
    encoding="utf-8",
)
(cwd / f"{prefix}_0001_ali.out").write_text("2 a b c d 12\\n", encoding="utf-8")
(cwd / f"{prefix}_0001_dp.ps").write_text("%!PS\\n", encoding="utf-8")
(cwd / f"{prefix}_0001_aln.ps").write_text("%!PS\\n", encoding="utf-8")
(cwd / f"{prefix}_0001_ss.ps").write_text("%!PS\\n", encoding="utf-8")
(cwd / "alirna.ps").write_text("%!PS stray\\n", encoding="utf-8")
print("fake RNAalifold run")
""",
        "ps2eps": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[1])
out = inp.with_suffix(".eps")
out.write_text("%!EPS\\n", encoding="utf-8")
""",
        "epstopdf": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[1])
out = inp.with_suffix(".pdf")
out.write_text("%PDF-FAKE\\n", encoding="utf-8")
""",
        "magick": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

out = Path(sys.argv[-1])
out.write_text("PNG-FAKE\\n", encoding="utf-8")
""",
        "R-scape": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[-1])
(Path.cwd() / f"{inp.stem}.power").write_text("# BPAIRS observed to covary 1\\n", encoding="utf-8")
""",
    }

    for name, content in tools.items():
        script = bin_dir / name
        script.write_text(textwrap.dedent(content), encoding="utf-8")
        script.chmod(0o755)


def _write_run_tree(root: Path, cleaned_sequence: str) -> None:
    """One window length's deterministic outputs, under whatever root is given."""
    split_text = "# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA ACGU\n//\n"
    (root / "Lalifold" / "len_150" / "split").mkdir(parents=True)
    (root / "generated_files" / "stk" / "len_150").mkdir(parents=True)
    (root / "Lalifold" / "len_150" / "RC_150_0001.stk").write_text(split_text, encoding="utf-8")
    for manifest in [
        root / "Lalifold" / "len_150" / "split" / "manifest.txt",
        root / "generated_files" / "stk" / "len_150" / "manifest.txt",
    ]:
        manifest.write_text("RC_150_0001_aln_1_10.stk\n", encoding="utf-8")
    (root / "Lalifold" / "len_150" / "split" / "RC_150_0001_aln_1_10.stk").write_text(
        split_text, encoding="utf-8"
    )
    (root / "generated_files" / "stk" / "len_150" / "RC_150_0001_aln_1_10.stk").write_text(
        f"# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA {cleaned_sequence}\n//\n",
        encoding="utf-8",
    )


def _run_verify(left: Path, right: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [PYTHON, "-m", "rnaconsnake.tools.verify_run_consistency", str(left), str(right)],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )


_REFOLD_RECORD = re.compile(r"^> \S+\n[ACGUN]+\n[().]+ \( *-?\d+\.\d\d\)$", re.MULTILINE)


FAKE_RANDOMIZE_ALN = """\
#!/usr/bin/env perl
# Fake rnazRandomizeAln.pl for CI: permutes alignment columns, which preserves
# per-sequence gap counts exactly like the real column shuffler does.
use strict;
use warnings;

my @names;
my %seq;
while (my $line = <>) {
    chomp $line;
    next if $line =~ /^CLUSTAL/i;
    next if $line =~ /^\\s*$/;
    my ($name, $chunk) = split /\\s+/, $line, 2;
    next unless defined $chunk;
    push @names, $name unless exists $seq{$name};
    $seq{$name} .= $chunk;
}
die "no sequences\\n" unless @names;
my $len = length($seq{$names[0]});
my @order = (0 .. $len - 1);
for (my $i = $#order; $i > 0; $i--) {
    my $j = int(rand($i + 1));
    @order[$i, $j] = @order[$j, $i];
}
print "CLUSTAL W (1.81) multiple sequence alignment\\n\\n\\n";
for my $name (@names) {
    my $out = join '', map { substr($seq{$name}, $_, 1) } @order;
    print "$name $out\\n";
}
"""


def write_fake_randomize_aln(bin_dir: Path) -> Path:
    script = bin_dir / "rnazRandomizeAln.pl"
    script.write_text(FAKE_RANDOMIZE_ALN, encoding="utf-8")
    script.chmod(0o755)
    return script


TOY_ALIGNMENT = "\n".join(
    [
        "# STOCKHOLM 1.0",
        "#=GF ID toy",
        "seqA GGCUAGCUAGCUAACGUAGCUAGCUAGGCAUCGAUCGAUCG",
        "seqB GGCUAGCUAG---ACGUAGCUAGCUAGGCAUCGAUCGAUCG",
        "seqC GGCUAGCUAGCUAACGUAGCUCGCUAGGCAUCGAUCG---G",
        "//",
        "",
    ]
)


def _window_rows(spans: list[tuple[int, int]], structures: dict[str, str] | None = None):
    rows = []
    for start, end in spans:
        name = f"RC_100_0001_aln_{start}_{end}"
        rows.append(
            {
                "wbn": name,
                "rnazprob": "0.5",
                "alifoldzscore": "-1.0",
                "alifold_consstruc": (structures or {}).get(name, ""),
            }
        )
    return rows


def _write_summary_csv(path: Path, rows: list[dict[str, str]]) -> str:
    from rnaconsnake.workflow_helpers import SUMMARY_FIELDS

    path.parent.mkdir(parents=True, exist_ok=True)
    import csv as _csv

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in SUMMARY_FIELDS})
    return str(path)


def _candidate_rows(prob: float, alifoldz: float, rscape: str, count: int = 4) -> list[dict[str, str]]:
    rows = []
    for index in range(count):
        start = 1 + index * 500
        rows.append(
            {
                "wbn": f"RC_100_0001_aln_{start}_{start + 99}",
                "rnazprob": f"{prob:.4f}",
                "alifoldzscore": f"{alifoldz:.4f}",
                "rscape_covary_count": rscape,
                "rscape_avg_confidence": "0.65" if rscape != "NA" else "",
                "rscape_mutual_info": "0.20" if rscape != "NA" else "",
                "nrseq": "6",
                "alilen": "100",
            }
        )
    return rows


BENCHMARK_TRUTH = "\n".join(
    [
        "element_id\telement_class\talignment\tstart\tend\tnotes",
        "xrRNA1\txrRNA\tflavivirus_3utr\t100\t199\tprimary target",
        "xrRNA2\txrRNA\tflavivirus_3utr\t900\t999\tsecond copy",
        "",
    ]
)


def _fake_viennarna_tools(tmp_path: Path, version: str) -> dict[str, str]:
    """Stand-ins for every ViennaRNA binary the check looks at, so the machine's
    own install does not decide the outcome."""
    tools = {}
    for key in cli.VIENNARNA_BINARIES:
        script = tmp_path / key
        script.write_text(f"#!/usr/bin/env python3\nprint('{key} {version}')\n", encoding="utf-8")
        script.chmod(0o755)
        tools[key] = str(script)
    return tools


LOWERCASE_CLUSTAL = "\n".join(
    [
        "CLUSTAL 2.1 multiple sequence alignment",
        "",
        "NMV_NC_032088.1 aggcacagaacgccg",
        "KOKV_NC_009029. cggcacagaacgccg",
        "ZIKV_NC_012532. aggcacagatcgccg",
        "                *********:*****",
        "",
    ]
)


ALIFOLDZ_EMPTY_REPORT = """\
###################################################################
# alifoldz.pl
#
#           Input: 0 sequences of 0 columns
#   Sample Number: 100
###################################################################

  From      To    Strand    Native MFE    Mean MFE     STDV      Z
 ------------------------------------------------------------------

9999
"""


ALIFOLDZ_GOOD_REPORT = """\
###################################################################
#           Input: 20 sequences of 52 columns
###################################################################

  From      To    Strand    Native MFE    Mean MFE     STDV      Z
 ------------------------------------------------------------------
     1      52       +         -4.64       -1.77       1.61    -1.8
-1.8
"""


def _run_extract_alifoldz(tmp_path: Path, report: str):
    source = tmp_path / "a.alifoldz.txt"
    source.write_text(report, encoding="utf-8")
    output = tmp_path / "a.alifoldz.json"
    return (
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.legacy_postprocess",
                "extract-alifoldz",
                "--input",
                str(source),
                "--output",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
            env=subprocess_env(),
        ),
        output,
    )


def _render_markdown(tmp_path: Path, nr_rows, full_rows, method="containment"):
    import csv as _csv

    from rnaconsnake.tools.dereplicate import NR_COLUMNS
    from rnaconsnake.workflow_helpers import SUMMARY_FIELDS

    nr = tmp_path / "RNAConSnake.nr.csv"
    with open(nr, "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=NR_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in nr_rows:
            writer.writerow({field: row.get(field, "") for field in NR_COLUMNS})

    full = tmp_path / "RNAConSnake.log.csv"
    with open(full, "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in full_rows:
            writer.writerow({field: row.get(field, "") for field in SUMMARY_FIELDS})

    out = tmp_path / "RNAConSnake.md"
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "render-markdown",
            "--label",
            "len_200",
            "--nr",
            str(nr),
            "--full",
            str(full),
            "--output",
            str(out),
            "--method",
            method,
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    return read_text(out)


def _write_nr_table(run_dir: Path, wlen: int, rows: list[dict[str, str]], method="containment"):
    import csv as _csv

    from rnaconsnake.tools.dereplicate import NR_COLUMNS

    target = run_dir / "generated_files" / "summary" / f"len_{wlen}"
    target.mkdir(parents=True, exist_ok=True)
    with open(target / "RNAConSnake.nr.csv", "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=NR_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in NR_COLUMNS})
    (target / "RNAConSnake.nr.json").write_text(
        json.dumps({"method": method, "n_loci": len(rows)}) + "\n", encoding="utf-8"
    )


WEB_FORBIDDEN_PATTERNS = [
    r"\bjinja2?\b",
    r"\bdatatables\b",
    r"\btabler\b",
    r"<!doctype html",
    r"<html[ >]",
    r"render_.*_pages",
    r"site_builder",
]


WEB_FORBIDDEN_PATHS = (".html", ".jinja", ".j2", ".css", ".js")


def _tracked_files() -> list[Path]:
    result = subprocess.run(["git", "ls-files"], check=True, capture_output=True, text=True, cwd=str(ROOT))
    return [Path(line) for line in result.stdout.splitlines() if line]


def _sweep_arm_inputs(tmp_path: Path, real_rows, null_rows_per_arm):
    arm_inputs = {"real": {100: _write_summary_csv(tmp_path / "real.csv", real_rows)}}
    for index, rows in enumerate(null_rows_per_arm):
        arm = f"null_{index:03d}"
        arm_inputs[arm] = {100: _write_summary_csv(tmp_path / f"{arm}.csv", rows)}
    return arm_inputs


def _sweep_base():
    from rnaconsnake.tools.calibration import Thresholds

    return Thresholds(0.9, -2.0, 1, 0.5, 0.1, 0.5, 1, 0.2, "containment", 0.9, 120, 0.8)


def _locus(locus_id, start, end, prob="0.99"):
    return {
        "locus_id": locus_id,
        "start": str(start),
        "end": str(end),
        "rnazprob": prob,
        "alifoldzscore": "-3.0",
        "q_cascade": "",
        "cascade_pass": "",
    }


def _two_elements(tmp_path: Path):
    from rnaconsnake.tools.benchmark import read_truth

    path = tmp_path / "truth.tsv"
    path.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "e1\txrRNA\taln\t100\t199\tfirst\n"
        "e2\tDB\taln\t400\t499\tsecond\n",
        encoding="utf-8",
    )
    return read_truth(path, "aln")


def _envelope_alignment():
    from rnaconsnake.tools.alignment_io import Alignment

    # Two near-identical sequences and two divergent ones, so subsets span a
    # range of mean pairwise identity.
    seqs = {
        "a": "ACGUACGUACGUACGUACGU",
        "b": "ACGUACGUACGUACGUACGA",
        "c": "AGGUAGGUAGGUAGGUAGGU",
        "d": "UCAUUCAUUCAUUCAUUCAU",
    }
    return Alignment(order=list(seqs), seqs=seqs)


def _refold_alignment(rows: dict[str, str]):
    from rnaconsnake.tools.alignment_io import Alignment

    return Alignment(order=list(rows), seqs=dict(rows))


FOLD_REGION_ALIGNMENT = "\n".join(
    [
        "# STOCKHOLM 1.0",
        "#=GF ID region_source",
        "seqA GGGCUAGCUAGGCAUCGAUCGGCUAGCUAGCCGAUCG",
        "seqB GGGCUAGCUAGGCAUCGAUC-GCUAGCUAGCCGAUCG",
        "seqC GGGCUAGCAAGGCAUCGAUCGGCUAGCUAGCCGAUCG",
        "//",
        "",
    ]
)


def _fold_region_tools(bin_dir: Path) -> list[str]:
    """Point every external tool at the fakes, by path rather than via PATH."""
    return [
        "--rnaalifold",
        str(bin_dir / "RNAalifold"),
        "--rnaz",
        str(bin_dir / "RNAz"),
        "--alifoldz",
        str(bin_dir / "alifoldz.pl"),
        "--eslreformat",
        str(bin_dir / "esl-reformat"),
        "--ps2eps",
        str(bin_dir / "ps2eps"),
        "--epstopdf",
        str(bin_dir / "epstopdf"),
    ]
