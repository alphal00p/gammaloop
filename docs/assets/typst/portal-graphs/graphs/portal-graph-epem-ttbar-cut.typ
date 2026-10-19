#import "../figure.typ": render

// Historical e+ e- -> t tbar illustration with explicit external legs.
// Preserve the drawing input from ae45d231fa54a22775923b3dcd993fba414e8dbb;
// canonical runtime fixtures now use sewn initial-state edges. This view does
// not apply the final-state physical cut.
#render(
  read("../fixtures/epem-ttbar-external-legs.dot"),
  cross-section-mode: true,
  momentum-arrows: true,
)
