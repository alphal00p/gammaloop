"""Private background worker for the explicitly requested gallery benchmark."""

import json
import os
import subprocess
import sys
from pathlib import Path
from time import monotonic

from shared_gallery import Gallery


def main():
    work = Path(sys.argv[1]).resolve()
    lock_fd = int(sys.argv[2])
    started = monotonic()

    def status(state, **fields):
        temporary = work / "status.tmp"
        temporary.write_text(
            json.dumps(
                {"state": state, "elapsed_seconds": monotonic() - started, **fields}
            )
        )
        temporary.replace(work / "status.json")

    def receipt():
        path = work / "benchmark.json"
        if not path.exists():
            return None
        if path.stat().st_size > Gallery.display_receipt_bytes:
            return {
                "artifact": str(path),
                "bytes": path.stat().st_size,
                "display": "Full receipt exceeds the browser display limit; scientific data remains in the artifact",
            }
        return Gallery.display_receipt(json.loads(path.read_text()))

    def log_tail():
        with (work / "run.log").open("rb") as stream:
            stream.seek(max(0, os.fstat(stream.fileno()).st_size - 5000))
            return stream.read().decode(errors="replace")

    try:
        command = json.loads((work / "request.json").read_text())["command"]
        status("running")
        with (work / "run.log").open("w") as log:
            result = subprocess.run(
                command,
                check=False,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        if result.returncode:
            status(
                "failed",
                exit_code=result.returncode,
                message=log_tail(),
                receipt=receipt(),
            )
        else:
            status("complete", receipt=receipt())
    except (OSError, ValueError) as error:
        status("failed", message=str(error))
    finally:
        os.close(lock_fd)


if __name__ == "__main__":
    main()
