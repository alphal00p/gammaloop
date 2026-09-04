#import "shared.typ": *
#import "components.typ": *

#let product-version(product) = {
  let source = product.rust_components.first().version_source
  toml("/" + source).package.version
}

#let software-citation(product, version) = {
  let citation = product.citation
  let doi = citation.at("doi", default: none)
  let target = if doi == none { citation.repository } else { "https://doi.org/" + doi }
  (
    citation.creators.join(", ") + " (" + str(citation.year) + "). " + citation.title
      + " (Version " + version + ") [Computer software]. " + target
  )
}

#let software-bibtex(product, version) = {
  let citation = product.citation
  let doi = citation.at("doi", default: none)
  let target = if doi == none { citation.repository } else { "https://doi.org/" + doi }
  let doi-line = if doi == none { "" } else { "  doi = {" + doi + "},\n" }
  (
    "@software{" + product.id + str(citation.year) + ",\n"
      + "  author = {" + citation.creators.join(" and ") + "},\n"
      + "  title = {" + citation.title + "},\n"
      + "  version = {" + version + "},\n"
      + "  year = {" + str(citation.year) + "},\n"
      + doi-line + "  url = {" + target + "}\n}"
  )
}

#let citation-entry(product) = {
  let version = product-version(product)
  let doi = product.citation.at("doi", default: none)
  citation-card(product.id)[
    #kicker[Version #version]

    == #product.title

    #if doi == none {
      citation-status[No registered software DOI is currently configured; this citation uses the versioned source repository.]
    } else {
      citation-status[
        #link("https://doi.org/" + doi)[doi:#doi]
      ]
    }

    === Suggested citation

    #copyable-code(
      "citation-" + product.id,
      software-citation(product, version),
      [Copy citation],
    )

    #disclosure([BibTeX])[
      #copyable-code(
        "bibtex-" + product.id,
        software-bibtex(product, version),
        [Copy BibTeX],
      )
    ]
  ]
}

#let citations = [
  #page-hero[
    #kicker[Credit the software]

    = Cite αLoop projects

    Use the version-specific records below. Where a Zenodo DOI exists, it is the persistent citation target; otherwise the citation names the versioned source repository without inventing an identifier.
  ]
  #citation-grid[
    #for product in registry.product { citation-entry(product) }
  ]
]

#let citations-page = portal-document(
  "Citations",
  "Version-specific citations for αLoop research software.",
  citations,
)
