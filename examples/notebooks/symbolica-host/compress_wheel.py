"""Compress the final Pyodide wheel with Zstandard after its ABI-tag rewrite."""

import argparse
import copy
from pathlib import Path
from zipfile import ZIP_ZSTANDARD, ZipFile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    output = args.destination / args.wheel.name
    temporary = output.with_suffix(".whl.tmp")

    with ZipFile(args.wheel) as source, ZipFile(temporary, "w") as compressed:
        compressed.comment = source.comment
        for original in source.infolist():
            entry = copy.copy(original)
            entry.compress_type = ZIP_ZSTANDARD
            entry.compress_level = 22
            # File bytes, including RECORD and its hashes, remain unchanged.
            compressed.writestr(entry, source.read(original))

    with ZipFile(temporary) as compressed:
        bad_entry = compressed.testzip()
        if bad_entry is not None:
            raise ValueError(f"Wheel verification failed for {bad_entry}")
    temporary.replace(output)
    print(f"Zstandard level 22 wheel: {output} ({output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
