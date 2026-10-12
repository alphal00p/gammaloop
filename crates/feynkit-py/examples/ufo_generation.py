from pathlib import Path

from symbolica.community.hepkit import Loops, Naming, NumeratorOptions, UfoLoader


loaded = UfoLoader().load(Path("/opt/ufo-models/sm"))
request = loaded.model.process(["e-", "e+"], ["mu-", "mu+"]).amplitude(
    loops=Loops.total(0),
)
result = request.generate().build(NumeratorOptions(naming=Naming.named("D")))

for diagram in result.diagrams:
    print(diagram.to_dot())
