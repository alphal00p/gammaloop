#import "shared.typ": *
#import "components.typ": *

#let about = [
  #about-hero[
    #kicker[About the collaboration]

    = Precision through local cancellation.

    #portal.summary
  ]
  #about-origin[
    #about-origin-copy[
      #kicker[Why αLoop]

      == Precision is another path to discovery.

      The lack of obvious sign of new physics phenomenon in collider experiments is an opportunity to take a step back and reflect on the amazing theory we have discovered so far: the Standard Model. In particular, we must now strive to make ever more precise predictions so as to hunt for indirect evidence of new physics in small departure from expectations.

      For this reason, our collaboration is dedicated to theoretical and algorithmic research for the automated computation of cross-sections in Quantum Field Theories at arbitrary perturbative orders. In particular, we develop a new theoretical framework called _Local Unitarity_ (LU) which approaches this problem from an unorthodox way, particularly suited to numerical computations.

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

      #anchored-heading(2, "about-pillars-title")[One research programme, connected structures.]
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
