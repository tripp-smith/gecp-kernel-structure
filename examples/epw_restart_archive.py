"""Export a supported model-plus-lossless EPW archive to native restart text."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "archive",
        type=Path,
        help="EPWRestartArchive JSON manifest; payload NPZ is loaded automatically",
    )
    parser.add_argument("output", type=Path, help="new native anisotropic file")
    args = parser.parse_args()
    if args.archive.suffix == ".npz":
        parser.error("pass the archive JSON manifest, not its NPZ payload")
    from kernelgecp.integrations.epw_restart_archive import EPWRestartArchive

    manifest = args.output.with_name(args.output.name + ".json")
    if args.output.exists() or manifest.exists():
        raise FileExistsError("native output or context manifest already exists")
    archive = EPWRestartArchive.load(args.archive)
    context = archive.write_native_restart(args.output)
    manifest.write_text(json.dumps(context, indent=2, allow_nan=False) + "\n")
    print(f"Native restart: {args.output}; bound context: {manifest}")
    print("Uses explicit lossless corrections; this is not a pure fitted restart.")
    print("Validate the unmodified native EPW restart before making a workflow claim.")


if __name__ == "__main__":
    main()
