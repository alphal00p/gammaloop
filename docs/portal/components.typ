// Presentation boundary for the collaboration portal. This is the only Typst
// source that knows the generated HTML structure and CSS selectors.

#let el(tag, class: none, attrs: (:), body) = {
  let attrs = if class == none { attrs } else { (class: class) + attrs }
  html.elem(tag, attrs: attrs, body)
}

#let void(tag, class: none, attrs: (:)) = {
  let attrs = if class == none { attrs } else { (class: class) + attrs }
  html.elem(tag, attrs: attrs)
}

#let arrow(symbol: "↗") = el("span", attrs: ("aria-hidden": "true"))[#symbol]
#let kicker(body) = el("p", class: "portal-kicker", body)
#let button(href, body, primary: false) = el(
  "a",
  class: "portal-button" + if primary { " portal-button-primary" } else { "" },
  attrs: (href: href),
  body,
)

// Page sources use these as Typst layout components. Low-level HTML stays here so
// headings, paragraphs, emphasis, links, and lists in the pages remain ordinary
// Typst content.
#let region(class: none, id: none, labelled-by: none, attrs: (:), body) = {
  if id != none { attrs.insert("id", id) }
  if labelled-by != none { attrs.insert("aria-labelledby", labelled-by) }
  el("section", class: class, attrs: attrs, body)
}

#let group(class: none, attrs: (:), body) = el("div", class: class, attrs: attrs, body)

#let panel(class: none, id: none, attrs: (:), body) = {
  if id != none { attrs.insert("id", id) }
  el("article", class: class, attrs: attrs, body)
}

#let page-hero(class: none, body) = {
  let classes = "portal-page-hero"
  if class != none { classes += " " + class }
  el("header", class: classes, body)
}

#let navigation(label, class: none, body) = el(
  "nav",
  class: class,
  attrs: ("aria-label": label),
  body,
)

#let styled-paragraph(class, body) = el("p", class: class, body)

#let anchored-heading(level, id, body) = el(
  "h" + str(level),
  attrs: (id: id),
  body,
)

#let copyable-code(id, value, label) = [
  #el("pre", attrs: (id: id))[#el("code")[#value]]
  #el("button", class: "portal-button", attrs: (
    type: "button",
    "data-copy-target": id,
  ))[#label]
]

#let disclosure(summary, body) = el("details")[
  #el("summary")[#summary]
  #body
]

// Semantic portal components. Authored pages use these names and never need to
// know the CSS selectors or incidental HTML wrappers behind them.
#let landing-hero(art, body) = region(
  class: "portal-hero portal-section",
  id: "overview",
)[
  #group(class: "portal-hero-copy")[#body]
  #art
]

#let hero-lede(body) = styled-paragraph("portal-lede", body)
#let hero-actions(body) = group(class: "portal-hero-actions", body)

#let portal-process-graph(graph) = el("span", class: "portal-process-graph")[
  #void("img", class: "portal-graph-theme portal-graph-theme-light", attrs: (
    src: "assets/graphs/portal-graph-" + graph + "-light.svg",
    alt: "",
  ))
  #void("img", class: "portal-graph-theme portal-graph-theme-dark", attrs: (
    src: "assets/graphs/portal-graph-" + graph + "-dark.svg",
    alt: "",
  ))
]

#let landing-art(graphs, body) = group(class: "portal-hero-art")[
  #group(class: "portal-wordmark", attrs: (
    "aria-label": "αLoop collaboration mark",
    role: "img",
  ))[]
  #group(class: "portal-graph-field", attrs: (
    role: "img",
    "aria-label": "A jumble of Feynman graphs rendered by GammaLoop from real process and test data",
  ))[
    #for graph in graphs { portal-process-graph(graph) }
  ]
  #el("p")[#body]
]

#let task-choice(product, href, order, label, product-line, tagline) = el(
  "a",
  class: "portal-task-link",
  attrs: ("data-product": product, href: href),
)[
  #el("span", class: "portal-task-order")[#order]
  #el("span", class: "portal-task-copy")[
    #el("strong")[#label]
    #el("span")[#product-line]
    #el("small")[#tagline]
  ]
  #el("span", class: "portal-task-arrow", attrs: ("aria-hidden": "true"))[→]
]

#let task-section(body) = region(
  class: "portal-section portal-task-chooser",
  id: "tasks",
  labelled-by: "tasks-title",
  body,
)
#let task-heading(body) = el("header", class: "portal-task-heading", body)
#let task-grid(body) = navigation("Documentation by task", class: "portal-task-grid", body)
#let heading-copy(body) = group(body)

#let projects-section(body) = region(
  class: "portal-section portal-projects",
  id: "projects",
  labelled-by: "projects-title",
  body,
)
#let portal-section-heading(body) = group(class: "portal-section-heading", body)
#let project-grid(body) = group(class: "portal-project-grid", body)
#let project-card(product, body) = panel(
  class: "portal-project-card",
  attrs: ("data-product": product),
  body,
)
#let project-meta(order, role) = group(class: "portal-project-meta")[
  #el("span")[#order]#el("span")[#role]
]
#let project-summary(body) = styled-paragraph("portal-project-summary", body)
#let package-list(title, packages) = el("div", class: "portal-packages", attrs: (
  "aria-label": title + " crates and modules",
))[
  #el("span", class: "portal-packages-label")[Crates & modules]
  #for package in packages { el("span", class: "portal-package")[#package] }
]
#let project-links(title, body) = navigation(
  title + " documentation",
  class: "portal-card-links",
  body,
)
#let project-overview-link(href, body) = el(
  "a",
  class: "portal-card-primary",
  attrs: (href: href),
  body,
)
#let project-cite-link(href, body) = el(
  "a",
  class: "portal-card-cite",
  attrs: (href: href),
  body,
)

#let funding-note(href, title-id: "funding-title", about: false, body) = el(
  "aside",
  class: "portal-funding" + if about { " about-funding" } else { "" },
  attrs: ("aria-labelledby": title-id),
)[
  #el("span", class: "portal-funding-mark", attrs: ("aria-hidden": "true"))[α]
  #group(class: "portal-funding-copy")[#body]
  #el("a", class: "portal-text-link", attrs: (href: href))[
    Funding record #arrow()
  ]
]

#let about-hero(body) = page-hero(class: "about-page-hero", body)
#let about-origin(body) = region(class: "about-origin", body)
#let about-origin-copy(body) = group(class: "about-origin-copy", body)
#let about-equation(graph-light, graph-dark, formula-light, formula-dark, note) = el(
  "aside",
  class: "about-equation",
  attrs: ("aria-label": "Schematic Local Unitarity cross-section"),
)[
  #group(class: "about-equation-illustration")[
    #void("img", class: "about-equation-graph portal-graph-theme-light", attrs: (
      src: graph-light,
      alt: "",
    ))
    #void("img", class: "about-equation-graph portal-graph-theme-dark", attrs: (
      src: graph-dark,
      alt: "",
    ))
  ]
  #group(class: "about-equation-formula", attrs: (
    role: "img",
    "aria-label": "The differential cross section is a sum over graphs of loop-momentum integrals and a sum over cuts of the Local Unitarity integrand, constrained by the observable.",
  ))[
    #void("img", class: "portal-graph-theme-light", attrs: (src: formula-light, alt: ""))
    #void("img", class: "portal-graph-theme-dark", attrs: (src: formula-dark, alt: ""))
  ]
  #el("small")[#note]
]
#let about-pillars(body) = region(
  class: "about-pillars",
  labelled-by: "about-pillars-title",
  body,
)
#let about-affiliations(body) = region(
  class: "about-affiliations",
  labelled-by: "about-affiliations-title",
  body,
)
#let section-header(body) = el("header", body)
#let pillar-grid(body) = group(body)
#let affiliation-grid(body) = group(body)
#let about-pillar(pillar) = panel(class: "about-pillar")[
  #kicker[#pillar.label]

  == #pillar.title

  #pillar.summary
]
#let about-affiliation(affiliation) = el("a", class: "about-affiliation", attrs: (
  href: affiliation.url,
))[
  #el("span")[#affiliation.url]
  #el("strong")[#affiliation.location]
  #el("small")[#affiliation.name]

  #affiliation.summary

  #el("b", attrs: ("aria-hidden": "true"))[↗]
]
#let about-next(body) = region(class: "about-next", body)

#let person-portrait(person) = if "portrait" in person {
  void("img", class: "people-card-portrait", attrs: (
    src: "../assets/people/" + person.portrait,
    alt: "",
    width: "720",
    height: "720",
    loading: "lazy",
  ))
} else {
  el("span", class: "people-card-initials", attrs: ("aria-hidden": "true"))[
    #person.initials
  ]
}
#let person-card(person, body) = panel(class: "people-card", id: person.id)[
  #person-portrait(person)
  #group(class: "people-card-copy")[#body]
]
#let people-grid(body) = region(class: "people-page-grid", body)

#let talks-hero(body) = page-hero(class: "talks-page-hero", body)
#let talks-provenance(body) = styled-paragraph("talks-provenance", body)
#let talk-card(id, datetime, date, year, body) = panel(class: "talk-card", id: id)[
  #group(class: "talk-card-date")[
    #el("time", attrs: (datetime: datetime))[#date]
    #el("span")[#year]
  ]
  #group(class: "talk-card-copy")[#body]
]
#let talk-timeline(body) = group(class: "talk-timeline", body)
#let talk-year(year, body) = region(
  class: "talk-year",
  labelled-by: "talk-year-" + year,
)[
  #anchored-heading(2, "talk-year-" + year)[#year]
  #group[#body]
]

#let publication-record(publication, body) = panel(class: "publication-card", attrs: (
  "data-publication": "",
  "data-title": publication.title,
  "data-people": publication.people.join(" "),
  "data-year": str(publication.year),
  "data-date": publication.date,
  "data-types": publication.types.join("|"),
  "data-citations": str(publication.citations),
), body)
#let publication-metadata(publication, venue) = group(class: "publication-card-meta")[
  #el("time", attrs: (datetime: publication.date))[#str(publication.year)]
  #if venue != none { el("span")[#venue] }
  #el("span")[#str(publication.citations) citations]
]
#let publication-authors(body) = styled-paragraph("publication-authors", body)
#let publication-provenance(body) = styled-paragraph("publication-provenance", body)
#let publication-filters(authors, years, kinds, kind-label, count) = el(
  "form",
  class: "publication-filters",
  attrs: ("data-publication-filters": ""),
)[
  #el("label")[Search #void("input", attrs: (
    type: "search",
    "data-publication-search": "",
    placeholder: "Title or author",
  ))]
  #el("label")[Author #el("select", attrs: ("data-publication-author": ""))[
    #el("option", attrs: (value: ""))[All authors]
    #for author in authors { el("option", attrs: (value: author.id))[#author.name] }
  ]]
  #el("label")[Year #el("select", attrs: ("data-publication-year": ""))[
    #el("option", attrs: (value: ""))[All years]
    #for year in years { el("option", attrs: (value: str(year)))[#str(year)] }
  ]]
  #el("label")[Type #el("select", attrs: ("data-publication-type": ""))[
    #el("option", attrs: (value: ""))[All types]
    #for kind in kinds { el("option", attrs: (value: kind))[#kind-label(kind)] }
  ]]
  #el("label")[Sort #el("select", attrs: ("data-publication-sort": ""))[
    #el("option", attrs: (value: "newest"))[Newest]
    #el("option", attrs: (value: "cited"))[Most cited]
  ]]
  #el("output", attrs: (
    "data-publication-count": "",
    "aria-live": "polite",
  ))[#str(count) publications]
]
#let publication-list(body) = region(
  class: "publication-list",
  attrs: ("data-publication-list": ""),
  body,
)
#let publication-no-script(body) = el("noscript")[
  #styled-paragraph("portal-page-note", body)
]

#let citation-card(id, body) = panel(class: "citation-card", id: id, body)
#let citation-status(body) = styled-paragraph("citation-status", body)
#let citation-grid(body) = region(class: "citation-grid", body)

#let nav-attrs(href, item, active) = if item == active {
  (href: href, "aria-current": "page")
} else {
  (href: href)
}

#let portal-header(root: "", home: "", active: none, search: false) = el(
  "header",
  class: "portal-header",
)[
  #el("a", class: "portal-brand", attrs: (href: home, "aria-label": "αLoop home"))[
    #el("span", class: "portal-brand-logo", attrs: ("aria-hidden": "true"))[]
    #el("span", class: "portal-brand-copy")[
      #el("strong")[αLoop]
      #el("small")[Local Unitarity research]
    ]
  ]
  #el("nav", class: "portal-nav", attrs: ("aria-label": "Primary"))[
    #el("a", attrs: nav-attrs(root + "about/", "about", active))[About]
    #el("a", attrs: (href: root + "#projects"))[Projects]
    #el("a", attrs: nav-attrs(root + "people/", "people", active))[People]
    #el("a", attrs: nav-attrs(root + "talks/", "talks", active))[Talks]
    #el("a", attrs: nav-attrs(root + "publications/", "publications", active))[Publications]
    #el("a", attrs: (href: root + "developers/"))[Developers]
  ]
  #el("div", class: "portal-header-actions")[
    #if search {
      el("button", class: "portal-search-button", attrs: (
        type: "button",
        "data-search-open": "",
      ))[
        Search #el("span", class: "header-button-label")[⌘K]
      ]
    }
    #el("a", class: "portal-source-link", attrs: (
      href: "https://github.com/alphal00p/gammaloop",
    ))[GitHub #arrow()]
    #el("button", class: "portal-theme-button", attrs: (
      type: "button",
      "data-theme-toggle": "",
      "aria-label": "Toggle color theme",
    ))[#el("span", attrs: ("aria-hidden": "true"))[◐]]
  ]
]

#let portal-footer(root: "") = el("footer", class: "portal-footer")[
  #el("div")[
    #el("span", class: "portal-footer-mark", attrs: ("aria-hidden": "true"))[]
    #el("p")[#el("strong")[αLoop] #el("br")[] Local Unitarity research software]
  ]
  #el("nav", attrs: ("aria-label": "Footer"))[
    #el("a", attrs: (href: root + "about/"))[About]
    #el("a", attrs: (href: root + "#projects"))[Projects]
    #el("a", attrs: (href: root + "people/"))[People]
    #el("a", attrs: (href: root + "talks/"))[Talks]
    #el("a", attrs: (href: root + "publications/"))[Publications]
    #el("a", attrs: (href: root + "citations/"))[Cite]
    #el("a", attrs: (href: root + "developers/"))[Developers]
    #el("a", attrs: (href: "https://github.com/alphal00p/gammaloop"))[Source]
  ]
  #el("p")[Physics, algorithms, and software #el("br")[] developed in the open.]
]

#let search-dialog = el("dialog", class: "search-dialog", attrs: (
  "data-search-dialog": "",
))[
  #el("form", class: "search-form", attrs: (method: "dialog"))[
    #void("input", class: "search-input", attrs: (
      type: "search",
      "data-search-input": "",
      placeholder: "Search all projects and developer notes",
      "aria-label": "Search all documentation",
    ))
    #el("button", class: "header-button", attrs: (value: "close"))[Close]
  ]
  #el("ul", class: "search-results", attrs: (
    "data-search-results": "",
    "aria-live": "polite",
  ))[]
]

#let portal-document(
  title,
  description,
  body,
  root: "../",
  home: "../",
  active: none,
  landing: false,
) = el("html", attrs: (lang: "en"))[
  #el("head")[
    #void("meta", attrs: (charset: "utf-8"))
    #void("meta", attrs: (
      name: "viewport",
      content: "width=device-width,initial-scale=1",
    ))
    #void("meta", attrs: (name: "description", content: description))
    #void("meta", attrs: (name: "theme-color", content: "#f9f6f0"))
    #void("link", attrs: (
      rel: "icon",
      type: "image/svg+xml",
      href: root + "assets/local-unitarity-light.svg",
    ))
    #void("link", attrs: (
      rel: "icon",
      type: "image/svg+xml",
      href: root + "assets/local-unitarity-dark.svg",
      media: "(prefers-color-scheme: dark)",
    ))
    #el("title")[#if landing { [αLoop · #title] } else { [#title · αLoop] }]
    #void("link", attrs: (rel: "stylesheet", href: root + "assets/site.css"))
    #el("script", attrs: (defer: "", src: root + "assets/site.js"))[]
  ]
  #el(
    "body",
    class: if landing { "portal-body" } else { "portal-body portal-subpage-body" },
    attrs: if landing {
      ("data-search-index": "search-index.json", "data-search-root": "")
    } else {
      (:)
    },
  )[
    #el("a", class: "skip-link", attrs: (href: "#main-content"))[Skip to content]
    #portal-header(root: root, home: home, active: active, search: landing)
    #el(
      "main",
      class: if landing { "portal-main" } else { "portal-main portal-subpage-main" },
      attrs: (id: "main-content"),
    )[
      // Typst reserves one HTML heading level for document metadata. The portal
      // supplies its own <title>, so preserve the levels authors write (`=` is
      // h1, `==` is h2, and so on).
      #show heading: it => html.elem("h" + str(it.level), it.body)
      #body
    ]
    #portal-footer(root: root)
    #if landing { search-dialog }
  ]
]
