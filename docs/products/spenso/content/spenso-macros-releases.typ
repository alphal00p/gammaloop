#import "../../shared.typ": source-link

#let spenso-macros-releases = [
= Spenso macros release history

The procedural macros are versioned independently from Spenso. Their canonical
history remains the #source-link("crates/spenso-macros/CHANGELOG.typ", label:
"spenso-macros changelog source"). This rendered copy integrates it with the
manual, navigation, search, and product PDF.

#[
  #show heading.where(level: 1): it => none
  #include "/crates/spenso-macros/CHANGELOG.typ"
]
]
