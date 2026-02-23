#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compatibility wrapper for legacy runner module name.

Use `unified_codebase_runfile.py` for all new development.
"""

try:
    from .unified_codebase_runfile import *  # noqa: F401,F403
except ImportError:
    from unified_codebase_runfile import *  # noqa: F401,F403
