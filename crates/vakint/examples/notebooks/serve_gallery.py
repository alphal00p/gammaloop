"""Serve the four interactive apps behind a tailnet-only reverse proxy."""

import hashlib
import os
from pathlib import Path

import marimo
import uvicorn
from shared_gallery import Gallery
from starlette.applications import Starlette
from starlette.responses import FileResponse, PlainTextResponse, RedirectResponse
from starlette.routing import Mount, Route


def evidence(request):
    try:
        resolved = Gallery.evidence_file(
            request.path_params["name"], request.path_params.get("artifact_index")
        )
        if resolved is None:
            return PlainTextResponse("Evidence not found", status_code=404)
        path, expected = resolved
        with path.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if actual != expected:
            return PlainTextResponse("Evidence checksum mismatch", status_code=409)
        return FileResponse(
            path, filename=path.name, headers={"X-Content-SHA256": actual}
        )
    except (OSError, ValueError, KeyError):
        return PlainTextResponse("Evidence unavailable", status_code=404)


def create_app():
    directory = Path(__file__).resolve().parent
    builder = marimo.create_asgi_app(include_code=False, show_tracebacks=False)
    for route, source in [
        ("/graphs", "01_graphs_and_numerators.py"),
        ("/inspect", "02_constants_and_masters.py"),
        ("/five-loop", "03_five_loop_uv.py"),
        ("/closure", "04_pinches_and_cache.py"),
    ]:
        builder = builder.with_app(path=route, root=str(directory / source))
    return Starlette(
        routes=[
            Route("/", lambda _: RedirectResponse("graphs/")),
            Route(
                "/assets/vakint-math.woff2",
                lambda _: FileResponse(
                    directory.parents[2] / "spynso3/typst/STIXTwoMath-Regular.woff2",
                    media_type="font/woff2",
                ),
            ),
            Route("/evidence/{name}", evidence),
            Route("/evidence/{name}/{artifact_index:int}", evidence),
            Mount("/", builder.build()),
        ]
    )


if __name__ == "__main__":
    uvicorn.run(
        create_app(),
        host="127.0.0.1",
        port=int(os.environ.get("VAKINT_GALLERY_PORT", "8417")),
    )
