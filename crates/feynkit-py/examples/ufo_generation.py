from pathlib import Path

from symbolica.community.feynkit import Generator, Process, UfoLoader


loaded = UfoLoader().load(Path("/opt/ufo-models/sm"))
result = Generator(loaded.model).generate(
    Process.amplitude(["e-", "e+"], ["mu-", "mu+"])
)

for diagram in result.diagrams:
    print(diagram.to_dot())
