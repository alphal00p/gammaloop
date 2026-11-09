#import "shared.typ": *
#import "components.typ": *

#let topic-entry(topic) = topic-card(topic.id)[
  === #topic.title

  #topic.summary

  #topic-links(topic.title)[
    #link(page-link(topic.guide))[Read the guide #arrow()]
    #link(page-link(topic.examples))[Examples #arrow()]
  ]
  #topic-libraries[
    #for (index, id) in topic.products.enumerate() {
      let product = registry.product.find(product => product.id == id)
      if index > 0 { [ · ] }
      link(product-root(id))[#product.title]
    }
  ]
]

#let landing = [
  #landing-hero(landing-art(portal.graphs))[
    #kicker[The #el("span", class: "keep-alpha", attrs: (style: "text-transform: none"))[α]Loop collaboration]

    = #portal.title

    #hero-lede[#portal.summary]

    #portal.libraries_summary

    #hero-actions[
      #button("#examples", [Explore the examples #arrow(symbol: "→")], primary: true)
      #button("libraries/", [Browse the libraries #arrow(symbol: "→")])
    ]
  ]
  #topic-section[
    #topic-heading[
      #anchored-heading(2, "topics-title")[Work through a calculation]
    ]
    #topic-grid[
      #for topic in portal.topic { topic-entry(topic) }
    ]
  ]
  #region(class: "portal-section portal-research", id: "research", labelled-by: "research-title")[
    #section-header[
      #anchored-heading(2, "research-title")[Local Unitarity]
    ]
    #about-origin[
      #about-origin-copy[
        Our research develops a shared numerical representation for perturbative
        cross-section calculations, with local cancellation of infrared singularities.

        #navigation("Local Unitarity research")[
          #button("about/", [About the research #arrow(symbol: "→")])
          #button("publications/", [Publications #arrow(symbol: "→")])
        ]
      ]
      #about-equation(
        "assets/about-double-triangle-light.svg",
        "assets/about-double-triangle-dark.svg",
        "assets/about-local-unitarity-equation-light.svg",
        "assets/about-local-unitarity-equation-dark.svg",
        [Real and virtual contributions share one numerical representation.],
      )
    ]
  ]
  #examples-section[
    #topic-heading[
      #anchored-heading(2, "examples-title")[Worked examples]
    ]
    #for example in portal.example {
      example-row(page-link(example.target), example.runtime)[
        #group[
          === #example.title

          #example.summary
        ]
        #example-runtime[
          #if example.runtime == "browser" { [Browser example] } else { [Run locally] }
          #arrow()
        ]
      ]
    }
  ]
  #library-note[
    FeynKit is a collection of tools for perturbative quantum field theory.
    Use its crates independently, alongside Linnet and Spenso, or together in GammaLoop.
    #link("libraries/")[Browse the libraries and API references #arrow()]
  ]
  #tour-section[
    #topic-heading[
      #anchored-heading(2, "tour-title")[See GammaLoop run]

      An animated walk through process generation, numerical integration, and the Rust and Python APIs. The commands are real; the terminal output is illustrative.
    ]
    #tour-frame("Animated tour of GammaLoop", "assets/showcase-hard-boiled.mp3")
    #tour-credit[
      Music: #link("https://incompetech.com/music/royalty-free/index.html?isrc=USUAN1700076")[Hard Boiled]
      by Kevin MacLeod, licensed under
      #link("https://creativecommons.org/licenses/by/4.0/")[CC BY 4.0].
    ]
  ]
  #funding-note(portal.funding_url)[
    #anchored-heading(2, "funding-title")[Research funding]

    #portal.funding
  ]
]

#let index-page = portal-document(
  "Software for perturbative quantum field theory",
  portal.summary,
  landing,
  root: "",
  home: "#overview",
  landing: true,
)
