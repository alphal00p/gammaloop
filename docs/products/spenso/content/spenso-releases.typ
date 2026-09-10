#import "../../shared.typ": source-link

#let spenso-releases = [
= Spenso versions

This page renders the canonical
#source-link("crates/spenso/CHANGELOG.typ", label: "Spenso changelog source") as part of the
Spenso website, navigation, and search.

#[
  #show heading.where(level: 1): it => none
  #include "/crates/spenso/CHANGELOG.typ"
]
]
