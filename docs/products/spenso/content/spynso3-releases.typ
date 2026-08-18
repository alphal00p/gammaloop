#import "../../shared.typ": source-link

#let spynso3-releases = [
= Spynso3 release history

The Python adapter is versioned independently from Spenso. Its canonical
history remains the #source-link("crates/spynso3/CHANGELOG.typ", label:
"spynso3 changelog source"). This rendered copy integrates it with the manual,
navigation, search, and product PDF.

#[
  #show heading.where(level: 1): it => none
  #include "/crates/spynso3/CHANGELOG.typ"
]
]
