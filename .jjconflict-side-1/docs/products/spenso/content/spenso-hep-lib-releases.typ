#import "../../shared.typ": source-link

#let spenso-hep-lib-releases = [
= Spenso HEP library versions

The HEP helper library is versioned independently from Spenso. Its canonical
version history remains the #source-link("crates/spenso-hep-lib/CHANGELOG.typ", label:
"spenso-hep-lib changelog source"). This rendered copy integrates it with the
manual, navigation, search, and product PDF.

#[
  #show heading.where(level: 1): it => none
  #include "/crates/spenso-hep-lib/CHANGELOG.typ"
]
]
