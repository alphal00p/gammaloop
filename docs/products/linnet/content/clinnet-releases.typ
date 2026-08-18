#import "../../shared.typ": source-link

#let clinnet-releases = [
= Clinnet release history

Clinnet is versioned independently from the Linnet library. Its canonical
history remains the
#source-link("crates/clinnet/CHANGELOG.typ", label: "Clinnet changelog source").
This rendered copy makes that history part of the manual, navigation, search,
and product PDF.

#[
  #show heading.where(level: 1): it => none
  #include "/crates/clinnet/CHANGELOG.typ"
]
]
