"""mashID: identify organisms from genome assemblies or raw reads using Mash."""

__version__ = "0.2.7"
__author__ = "Marc-Olivier Duceppe"


class MashIDError(Exception):
    """A user-facing error: bad input, missing tool, or a failed external command."""
