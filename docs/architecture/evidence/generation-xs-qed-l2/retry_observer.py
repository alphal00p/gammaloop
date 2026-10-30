"""Historical observer API adaptation only; original failed build retained."""

import build
import common as M

build.compile_one("historical", "oracle_xs_qed_l2_capture", "example")
for label, argv in (("rustc", ["rustc", "-vV"]), ("cargo", ["cargo", "--version"])):
    proc, _ = M.run(argv, M.OUT / "build/historical" / label)
    M.require(proc.returncode == 0, "toolchain", label)
