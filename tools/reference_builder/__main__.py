"""Entry point for `python -m reference_builder`."""
import sys
from .cli import main

if __name__ == "__main__":
    sys.exit(main())