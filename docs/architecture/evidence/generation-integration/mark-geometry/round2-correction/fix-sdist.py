from pathlib import Path
p=Path('pyproject.toml');s=p.read_text().replace('include = [\n','include = [\n  { path = "vendor/kurbo/**", format = "sdist" },\n',1);p.write_text(s)
p=Path('tests/packaging/test_drawing_notices.py');s=p.read_text();anchor='    def test_api_embedded_payload(self):';idx=s.index(anchor);part=s[:idx];a=part.rfind('\n    def test_gammaloop_sdist_payload');assert a>=0;last=part[a:];needle='            with tarfile.open(archives[0]) as archive:\n                prefix = archive.getnames()[0].split("/")[0] + "/"\n';assert needle in last;last=last.replace(needle,needle+'''                # Cargo parses a root path patch even when the selected API
                # feature set only embeds prebuilt drawing WASM. Preserve its
                # complete published source, metadata and both license texts.
                for source in (ROOT / "vendor/kurbo").rglob("*"):
                    if source.is_file():
                        relative = source.relative_to(ROOT).as_posix()
                        with self.subTest(vendored_source=relative):
                            self.assertEqual(
                                archive.extractfile(prefix + relative).read(),
                                source.read_bytes(),
                            )
''',1);s=part[:a]+last+s[idx:];p.write_text(s)
