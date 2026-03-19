from __future__ import annotations

from pathlib import Path
from shutil import copy2

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


ROOT = Path(__file__).resolve().parent
PACKAGE_WORKFLOW_DIR = Path("rnaconsnake") / "workflow"
WORKFLOW_SOURCES = {
    ROOT / "snakefile": PACKAGE_WORKFLOW_DIR / "snakefile",
    ROOT / "config.yaml": PACKAGE_WORKFLOW_DIR / "config.yaml",
}


class build_py(_build_py):
    def run(self):
        super().run()
        for src, relative_dest in WORKFLOW_SOURCES.items():
            dest = Path(self.build_lib) / relative_dest
            dest.parent.mkdir(parents=True, exist_ok=True)
            copy2(src, dest)


setup(cmdclass={"build_py": build_py})
