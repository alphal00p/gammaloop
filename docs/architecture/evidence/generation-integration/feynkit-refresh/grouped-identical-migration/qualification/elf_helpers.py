"""Unchanged pure ELF parser and predicate from qualified runtime common.py."""
import struct
from pathlib import Path

def require(ok, *detail):
    if not ok:
        raise RuntimeError(": ".join(map(str, detail)))

def elf(path):
    """Read ELF metadata without executing loader/readelf/ldd or searching host libraries."""
    path = Path(path).resolve()
    data = path.read_bytes()
    if data[:4] != b"\x7fELF":
        return {"path": str(path), "non_elf": True, "candidates": []}
    require(data[4:6] == b"\x02\x01", "only measured ELF64 LE tools supported")
    offset = struct.unpack_from("<Q", data, 32)[0]
    size, count = struct.unpack_from("<HH", data, 54)
    segments = [
        struct.unpack_from("<IIQQQQQQ", data, offset + i * size) for i in range(count)
    ]
    dynamic, interpreter = [], None
    for typ, _, off, _, _, length, _, _ in segments:
        if typ == 3:
            interpreter = data[off : off + length].rstrip(b"\0").decode()
        if typ == 2:
            dynamic = [
                struct.unpack_from("<qQ", data, j) for j in range(off, off + length, 16)
            ]
    strings = next((v for t, v in dynamic if t == 5), None)
    if strings is None:
        return {
            "path": str(path),
            "interpreter": interpreter,
            "needed": [],
            "runpath": [],
            "candidates": [],
        }
    base = next(
        off + strings - va
        for typ, _, off, va, _, length, _, _ in segments
        if typ == 1 and va <= strings < va + length
    )
    text = lambda off: data[base + off : data.index(b"\0", base + off)].decode()
    paths = [
        p.replace("$ORIGIN", str(path.parent))
        for tag, off in dynamic
        if tag in (15, 29)
        for p in text(off).split(":")
    ]
    return {
        "path": str(path),
        "interpreter": interpreter,
        "needed": [text(off) for tag, off in dynamic if tag == 1],
        "runpath": paths,
        "candidates": [p for p in paths if not p.startswith("/nix/store/")],
    }
