"""Installed native hosts find resource wheels without environment overrides."""

import os
import subprocess
import sys
import unittest


@unittest.skipIf(sys.platform == "emscripten", "native installed-wheel discovery")
class InstalledDataTests(unittest.TestCase):
    def test_font_discovery(self):
        # Separate processes prevent the verified-byte cache hiding discovery bugs.
        for layout in ("installed", "temporary"):
            with self.subTest(layout=layout):
                environment = os.environ.copy()
                environment.pop("SYMBOLICA_DATA_DIRECTORY", None)
                subprocess.run(
                    [sys.executable, "-c", SCRIPT, layout],
                    env=environment,
                    check=True,
                )


SCRIPT = r'''
import shutil
import sys
import tempfile
from pathlib import Path

from symbolica import data
from symbolica.community import hepkit

name = "render/fonts/STIXTwoMath-Regular.woff2"
record = next(asset for asset in data.manifest()["assets"] if asset["name"] == name)
assert not record["loaded"], record
with tempfile.TemporaryDirectory() as temporary:
    if sys.argv[1] == "temporary":
        package = Path(temporary) / "symbolica"
        source = Path(data.__file__).parent / "_data/fonts" / name
        target = package / "_data/fonts" / name
        target.parent.mkdir(parents=True)
        shutil.copyfile(source, target)
        data.__file__ = str(package / "data.py")
        data.configure()
    view = hepkit.Model.phi3().process(["phi"], ["phi"]).render()
    html = view.to_html()
    request = view.measurement_request()
    assert request["labels"], request
    assert "@font-face" in request["font_css"]
    assert "<svg" in html
    assert next(asset for asset in data.manifest()["assets"] if asset["name"] == name)["loaded"]
    print(f"{sys.argv[1]}: labelled to_html() and measurement_request() discovered {record['size']} font bytes")
'''
