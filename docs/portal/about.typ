#import "shared.typ": *
#import "components.typ": *

#let about = [
  #about-hero[
    #kicker[About the collaboration]

    = Research on Local Unitarity

    #portal.summary
  ]
  #about-origin[
    #about-origin-copy[
      #kicker[Why αLoop]

      == Precision is another path to discovery.

      Precision predictions test the Standard Model and help us search for indirect evidence of new physics through small departures from its expectations.

      Our collaboration develops theoretical methods and algorithms for the automated computation of cross-sections in quantum field theory at arbitrary perturbative orders. _Local Unitarity_ (LU) puts real and virtual contributions into a shared numerical representation, so infrared singularities cancel locally before integration.

      #navigation("Learn about Local Unitarity")[
        #button("https://arxiv.org/abs/2110.15662", [Read the introduction #arrow()], primary: true)
        #button("../publications/", [Explore publications #arrow(symbol: "→")])
      ]
    ]
    #about-equation(
      "../assets/about-double-triangle-light.svg",
      "../assets/about-double-triangle-dark.svg",
      "../assets/about-local-unitarity-equation-light.svg",
      "../assets/about-local-unitarity-equation-dark.svg",
      [Real and virtual contributions share one numerical representation.],
    )
  ]
  #about-pillars[
    #section-header[
      #kicker[From method to software]

      #anchored-heading(2, "about-pillars-title")[From Local Unitarity to software]
    ]
    #pillar-grid[#for pillar in portal.pillar { about-pillar(pillar) }]
  ]
  #about-affiliations[
    #section-header[
      #kicker[Affiliations]

      #anchored-heading(2, "about-affiliations-title")[Research across institutions.]

      αLoop connects collider-physics research, mathematical structures, and open scientific-software development.
    ]
    #affiliation-grid[
      #for affiliation in portal.affiliation { about-affiliation(affiliation) }
    ]
  ]
  #funding-note(
    portal.funding_url,
    title-id: "about-funding-title",
    about: true,
  )[
      #kicker[Funding]

      #anchored-heading(2, "about-funding-title")[Publicly funded research]

      #portal.funding
  ]
  #about-next[
    #kicker[The collaboration]

    == Meet the people doing the work.

    #navigation("Explore the collaboration")[
      #button("../people/", [People #arrow(symbol: "→")], primary: true)
      #button("../talks/", [Talks #arrow(symbol: "→")])
    ]
  ]
]

#let about-page = portal-document(
  "About",
  "The αLoop collaboration develops Local Unitarity methods and open research software for precision collider physics.",
  about,
  active: "about",
)
