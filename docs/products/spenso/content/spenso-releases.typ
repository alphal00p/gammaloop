#import "../../shared.typ": source-link, release-note

#let spenso-releases = [
= Spenso versions

#release-note([
  This documentation covers Spenso 0.6.0. The latest recorded release below is 0.5.6;
  the canonical changelog has no 0.6.0 or 0.5.4 entry. Consult the source revision
  when assessing changes outside these recorded notes.
])

Release notes: #source-link("crates/spenso/CHANGELOG.typ", label: "Spenso changelog source").

#[
  #show heading.where(level: 1): it => none
  #show "All notable changes to this project will be documented in this file.": none
  #include "/crates/spenso/CHANGELOG.typ"
]
]
