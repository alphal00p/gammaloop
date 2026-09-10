#import "../../shared.typ": source-link

#let idenso-releases = [
= Idenso versions

This page renders the canonical
#source-link("crates/idenso/CHANGELOG.typ", label: "Idenso changelog source") as part of the
Idenso website, navigation, and search.

#[
  #show heading.where(level: 1): it => none
  #include "/crates/idenso/CHANGELOG.typ"
]
]
