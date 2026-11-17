"""Find the CalculiX executable and make it runnable.

On Windows the conda build of ccx needs DLLs from ``Library/mingw-w64/bin``
(libarpack, gfortran runtime). When the env's python is called without
activating the env, neither folder is on PATH, so ``shutil.which("ccx")``
fails or ccx exits with code 127. ``setup_ccx()`` fixes PATH for this process.

Order of lookup: the ``CCX`` environment variable, then ccx next to the
running interpreter's conda env, then PATH.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _env_bins() -> list[Path]:
    root = Path(sys.prefix)
    return [root / "Library" / "bin", root / "Library" / "mingw-w64" / "bin"]


def setup_ccx() -> str:
    """Return the path of ccx after putting its DLL folders on PATH.

    Raises:
        FileNotFoundError: if no ccx executable is found.
    """
    bins = [str(p) for p in _env_bins() if p.is_dir()]
    path = os.environ.get("PATH", "").split(os.pathsep)
    missing = [b for b in bins if b not in path]
    if missing:   # idempotent: repeated calls must not grow PATH (32767-char limit on Windows)
        os.environ["PATH"] = os.pathsep.join(missing + path)
    cand = os.environ.get("CCX") or shutil.which("ccx")
    if not cand:
        raise FileNotFoundError("ccx not found: set CCX or install calculix in this env")
    return cand


def ccx_version(ccx: str | None = None) -> str:
    """Return the version line printed by ``ccx -v``."""
    ccx = ccx or setup_ccx()
    out = subprocess.run([ccx, "-v"], capture_output=True, text=True, timeout=60)
    txt = out.stdout + out.stderr
    return next((ln.strip() for ln in txt.splitlines() if "ersion" in ln), txt.strip()[:80])


def run_ccx(job: str, cwd: Path, threads: int = 4, timeout: int = 3600) -> subprocess.CompletedProcess:
    """Run ``ccx -i job`` in ``cwd`` (``job`` without the .inp suffix)."""
    ccx = setup_ccx()
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    return subprocess.run([ccx, "-i", job], cwd=cwd, capture_output=True, text=True,
                          timeout=timeout, env=env)
