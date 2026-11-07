#import "../../shared.typ": source-link

#let linnet-releases = [
= Linnet versions

Release notes: #source-link("crates/linnet/CHANGELOG.typ", label: "Linnet changelog source").

#[
  #show heading.where(level: 1): it => none
  #show "All notable changes to this project will be documented in this file.": none
  #include "/crates/linnet/CHANGELOG.typ"
]
]
