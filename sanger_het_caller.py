#!/usr/bin/env python
"""Entry point: run as `python sanger_het_caller.py <file.ab1> [options]`."""

import sys

from het_caller.cli import main

if __name__ == "__main__":
    sys.exit(main())
