#import "shared.typ": *
#import "components.typ": *

#let profile-links(person) = [
  #link(person.url)[Professional profile #arrow()]
  #link(person.github)[GitHub #arrow()]
  #if "inspire_recid" in person {
    link("https://inspirehep.net/authors/" + str(person.inspire_recid))[
      INSPIRE HEP #arrow()
    ]
  }
  #if "orcid" in person {
    link("https://orcid.org/" + person.orcid)[ORCID #arrow()]
  }
]

#let person-entry(person) = person-card(person)[
  #kicker[Collaboration]

  == #person.name

  #person.role

  #navigation(person.name + " profiles")[#profile-links(person)]
]

#let people = [
  #page-hero[
    #kicker[People]

    = People building GammaLoop

    Researchers and collaborators developing GammaLoop, Local Unitarity methods, and the scientific software that supports them.
  ]
  #people-grid[
    #for person in portal.people { person-entry(person) }
  ]
]

#let people-page = portal-document(
  "People",
  "Researchers and collaborators building GammaLoop.",
  people,
  active: "people",
)
