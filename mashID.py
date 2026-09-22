#!/usr/bin/env python3
"""Backward-compatible entry point: ``python mashID.py ...``.

After ``pip install .`` the ``mashID`` command is available directly.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mashid.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
