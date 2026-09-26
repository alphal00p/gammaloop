#import "shared.typ": *
#import "components.typ": *

#let ordinal(index) = if index < 9 { "0" + str(index + 1) } else { str(index + 1) }

#let task-link(task, index) = {
  let product = registry.product.find(product => product.id == task.product)
  let quickstart = product.pages.find(page => page.id == "quickstart")
  let href = "products/" + product.id + "/" + channel-route + "/" + quickstart.route
  task-choice(
    product.id,
    href,
    ordinal(index),
    task.label,
    [#product.title · #task.role],
    product.tagline,
  )
}

#let project-entry(product, index) = {
  let task = portal.task.find(task => task.product == product.id)
  let components = product.rust_components + product.python_components
  let packages = components.map(component => component.package).sorted().dedup()
  let quickstart = product.pages.find(page => page.id == "quickstart")
  let guide = product.pages.find(page => page.group == "Guides")
  let root = "products/" + product.id + "/" + channel-route + "/"
  project-card(product.id)[
    #project-meta(ordinal(index), task.role)

    === #link(root)[#product.title]

    #project-summary[#product.tagline]

    #package-list(product.title, packages)

    #project-links(product.title)[
      #project-overview-link(root)[Overview #arrow()]
      #link(root + quickstart.route)[Get started]
      #link(root + guide.route)[Guides]
      #link(root + "reference/")[Reference]
    ]
  ]
}

#let funding = funding-note(portal.funding_url)[
  #kicker[Funding]

  #anchored-heading(2, "funding-title")[Publicly funded research]

  #portal.funding
]

#let tour = tour-section[
  #portal-section-heading[
    #heading-copy[
      #kicker[Guided tour]

      #anchored-heading(2, "tour-title")[See GammaLoop run]
    ]

    An animated walk through the command line, process generation, Monte Carlo integration, the Rust and Python APIs, and the crates behind them, drawn live in the browser. The commands are real; the terminal output is illustrative.
  ]
  #tour-frame("Animated tour of GammaLoop", "assets/showcase-hard-boiled.mp3")
  #tour-credit[
    Music: #link("https://incompetech.com/music/royalty-free/index.html?isrc=USUAN1700076")[Hard Boiled]
    by Kevin MacLeod (incompetech.com), licensed under
    #link("https://creativecommons.org/licenses/by/4.0/")[Creative Commons Attribution 4.0].
  ]
]

#let landing = [
  #landing-hero(landing-art(portal.graphs)[
    Local cancellation. \
    Global precision.
  ])[
    #kicker[#portal.eyebrow]

    = #portal.title

    #hero-lede[#portal.summary]

    #hero-actions[
      #button("#tasks", [Choose a task #arrow(symbol: "↓")], primary: true)
      #button("products/gammaloop/" + channel-route + "/quickstart/", [
        Start with GammaLoop #arrow()
      ])
    ]
  ]
  #task-section[
    #task-heading[
      #heading-copy[
        #kicker[Choose by task · 01—05]

        #anchored-heading(2, "tasks-title")[What do you want to work on?]
      ]

      Start with the scientific object or operation you have in hand. Each route opens the maintained first workflow for the component that owns it.
    ]
    #task-grid[
      #for (index, task) in portal.task.enumerate() { task-link(task, index) }
    ]
  ]
  #projects-section[
    #portal-section-heading[
      #heading-copy[
        #kicker[Research software · 01—05]

        #anchored-heading(2, "projects-title")[Projects & crates]
      ]

      Five connected codebases spanning numerical cross-sections, graph algorithms, tensor networks, symbolic identities, and integral evaluation.
    ]
    #project-grid[
      #for (index, product) in registry.product.enumerate() { project-entry(product, index) }
    ]
  ]
  #tour
  #funding
]

#let index-page = portal-document(
  "Research software for collider physics",
  portal.summary,
  landing,
  root: "",
  home: "#overview",
  landing: true,
)
