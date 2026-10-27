from pathlib import Path

from symbolica.community.hepkit import UfoLoader


loaded = UfoLoader().load(Path("/opt/ufo-models/sm"))
result = loaded.model.process(["e-", "e+"], ["mu-", "mu+"]).generate_diagrams()

for diagram in result.diagrams:
    print(diagram.to_dot())
