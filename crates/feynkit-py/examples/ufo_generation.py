from pathlib import Path

from symbolica.community.feynkit import UfoLoader


loaded = UfoLoader().load(Path("/opt/ufo-models/sm"))
result = loaded.model.generate_diagrams(
    incoming=["e-", "e+"],
    outgoing=["mu-", "mu+"],
)

for diagram in result.diagrams:
    print(diagram.to_dot())
