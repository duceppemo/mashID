"""``mashID_download_db``: fetch pre-built databases."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from mashid import MashIDError, __version__
from mashid.databases import ENV_DB_DIR, REGISTRY, db_dir, download, list_rows
from mashid.pipeline import format_table

log = logging.getLogger("mashid.download")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mashID_download_db",
        description="Download pre-built mashID databases. Without arguments, lists the available databases.",
    )
    parser.add_argument("names", nargs="*", metavar="NAME", help=f"Database(s) to download: {', '.join(REGISTRY)}")
    parser.add_argument("--all", action="store_true", help="Download every database in the registry.")
    parser.add_argument("--dir", metavar="DIR", type=Path, default=None,
                        help=f"Where to store databases. Default: ${ENV_DB_DIR} if set, else {db_dir()}")
    parser.add_argument("--force", action="store_true", help="Re-download even if already installed.")
    parser.add_argument("--no-verify", action="store_true", help="Skip MD5 verification.")
    parser.add_argument("-v", "--version", action="version", version=f"mashID_download_db {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
    names = list(REGISTRY) if args.all else args.names
    if not names:
        print(format_table(["Name", "Installed", "Size", "Description"], list_rows(args.dir)))
        print(f"\nDatabase directory: {args.dir or db_dir()}")
        print("Download with: mashID_download_db <NAME>   then use: mashID -d <NAME> ...")
        return 0
    try:
        for name in names:
            download(name, args.dir, force=args.force, verify=not args.no_verify)
    except MashIDError as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        log.error("Interrupted")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
