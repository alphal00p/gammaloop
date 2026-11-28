#!/usr/bin/env python3
"""Static server for the gamma5 documentation site.

Serves the directory `site/` next to this file on the loopback interface only. A reverse proxy in front
(`tailscale serve`, a Cloudflare tunnel) makes it reachable from elsewhere. Text resources are sent gzip-compressed
when the client accepts it; there is no directory listing and no write access.

    python3 serve.py [--port 21550] [--host 127.0.0.1] [--root site]
"""

from __future__ import annotations

import argparse
import gzip
import mimetypes
import posixpath
import sys
import urllib.parse
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

COMPRESSIBLE = (".html", ".js", ".mjs", ".css", ".json", ".svg", ".txt", ".md", ".toml", ".dot", ".typ", ".bib", ".py")
EXTRA_TYPES = {
    ".mjs": "text/javascript",
    ".js": "text/javascript",
    ".json": "application/json",
    ".svg": "image/svg+xml",
    ".woff2": "font/woff2",
    ".typ": "text/plain",
    ".toml": "text/plain",
    ".dot": "text/plain",
    ".bib": "text/plain",
    ".md": "text/plain",
    ".py": "text/plain",
}
# (path, mtime, size) -> compressed bytes
_gzip_cache: dict[tuple[str, float, int], bytes] = {}


class Handler(BaseHTTPRequestHandler):
    server_version = "gamma5docs/1"
    root: Path = Path(".")

    def log_message(self, fmt: str, *args) -> None:  # one line per request on stderr
        sys.stderr.write("%s %s\n" % (self.log_date_time_string(), fmt % args))

    def _resolve(self) -> Path | None:
        path = urllib.parse.unquote(urllib.parse.urlsplit(self.path).path)
        path = posixpath.normpath(path)
        parts = [p for p in path.split("/") if p and p not in (".", "..")]
        if any(p.startswith(".") for p in parts):
            return None
        target = self.root.joinpath(*parts)
        if target.is_dir():
            target = target / "index.html"
        try:
            resolved = target.resolve()
            resolved.relative_to(self.root)
        except (OSError, ValueError):
            return None
        return resolved if resolved.is_file() else None

    def _send(self, head_only: bool) -> None:
        target = self._resolve()
        if target is None:
            body = b"not found\n"
            self.send_response(HTTPStatus.NOT_FOUND)
            self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if not head_only:
                self.wfile.write(body)
            return
        suffix = target.suffix.lower()
        ctype = EXTRA_TYPES.get(suffix) or mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype in ("application/json", "image/svg+xml"):
            ctype += "; charset=utf-8"
        stat = target.stat()
        data = target.read_bytes()
        encoding = None
        if suffix in COMPRESSIBLE and "gzip" in self.headers.get("Accept-Encoding", "") and len(data) > 1024:
            key = (str(target), stat.st_mtime, stat.st_size)
            packed = _gzip_cache.get(key)
            if packed is None:
                packed = gzip.compress(data, compresslevel=6)
                _gzip_cache[key] = packed
            data, encoding = packed, "gzip"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        if encoding:
            self.send_header("Content-Encoding", encoding)
        self.send_header("Vary", "Accept-Encoding")
        # The site is rebuilt in place: let browsers revalidate rather than keep stale data.
        self.send_header("Cache-Control", "no-cache" if suffix in (".html", ".json", ".js", ".mjs", ".css") else "max-age=3600")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        if not head_only:
            self.wfile.write(data)

    def do_GET(self) -> None:  # noqa: N802 (http.server naming)
        self._send(head_only=False)

    def do_HEAD(self) -> None:  # noqa: N802
        self._send(head_only=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--port", type=int, default=21550)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--root", default=str(Path(__file__).resolve().parent / "site"))
    args = parser.parse_args()
    Handler.root = Path(args.root).resolve()
    if not Handler.root.is_dir():
        sys.exit(f"site directory not found: {Handler.root}")
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    server.daemon_threads = True
    print(f"serving {Handler.root} on http://{args.host}:{args.port}/", file=sys.stderr, flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
