#import "../../shared.typ": source-link

#let linnet-releases = [
= Linnet versions

This page renders the canonical
#source-link("crates/linnet/CHANGELOG.typ", label: "Linnet changelog source") as part of the
Linnet website, navigation, and search.

#[
  #show heading.where(level: 1): it => none
  #include "/crates/linnet/CHANGELOG.typ"
]
]
