#import "../../shared.typ": source-link

#let idenso-releases = [
= Idenso versions

Release notes: #source-link("crates/idenso/CHANGELOG.typ", label: "Idenso changelog source").

#[
  #show heading.where(level: 1): it => none
  #show "All notable changes to this project will be documented in this file.": none
  #include "/crates/idenso/CHANGELOG.typ"
]
]
