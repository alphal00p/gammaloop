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
  == #person.name

  #person.role

  #navigation(person.name + " profiles")[#profile-links(person)]
]

#let people = [
  #page-hero(class: "portal-page-hero-compact")[
    = People
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
