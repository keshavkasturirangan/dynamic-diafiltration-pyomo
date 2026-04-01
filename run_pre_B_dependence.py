#!/usr/bin/env python3
"""Compatibility wrapper for relocated legacy script."""

from pathlib import Path
import runpy

_TARGET = Path(__file__).resolve().parent / "legacy" / "scripts" / "run_pre_B_dependence.py"
runpy.run_path(str(_TARGET), run_name="__main__")
