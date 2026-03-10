#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compatibility wrapper for legacy module name.

Use `unified_codebase_library.py` for all new development.
"""

try:
    from .unified_codebase_library import *  # noqa: F401,F403
except ImportError:
    from unified_codebase_library import *  # noqa: F401,F403
