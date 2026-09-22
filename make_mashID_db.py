#!/usr/bin/env python3
"""Backward-compatible entry point: ``python make_mashID_db.py ...``.

After ``pip install .`` the ``make_mashID_db`` command is available directly.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mashid.makedb import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
