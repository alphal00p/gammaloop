#import "shared.typ": *
#import "components.typ": *

#let compact-authors(authors) = if authors.len() <= 8 {
  authors.join(", ")
} else {
  authors.slice(0, 6).join(", ") + ", et al."
}

#let publication-type-label(kind) = (
  article: "Article",
  "book chapter": "Book Chapter",
  "conference paper": "Conference Paper",
  proceedings: "Proceedings",
  report: "Report",
  thesis: "Thesis",
).at(kind, default: kind)

#let publication-entry(publication) = {
  let venue = publication.at("venue", default: none)
  let doi = publication.at("doi", default: none)
  let arxiv = publication.at("arxiv", default: none)
  publication-record(publication)[
    #publication-metadata(publication, venue)

    == #link(publication.url)[#publication.title]

    #publication-authors[#compact-authors(publication.authors)]

    #navigation("Citation links")[
      #link(publication.url)[INSPIRE]
      #if doi != none {
        link("https://doi.org/" + doi)[DOI #arrow()]
      }
      #if arxiv != none {
        link("https://arxiv.org/abs/" + arxiv)[arXiv #arrow()]
      }
      #link(publication.bibtex_url)[BibTeX]
    ]
  ]
}

#let years = publication-cache.publications.map(publication => publication.year).dedup()
#let publication-types = (
  publication-cache.publications
    .map(publication => publication.types)
    .flatten()
    .sorted()
    .dedup()
)

#let publications = [
  #page-hero[
    #kicker[Research output]

    = Publications

    Automatically assembled from stable INSPIRE HEP author identifiers, deduplicated across coauthors, and cached for a reproducible documentation build.

    #publication-provenance[
      Updated #publication-cache.updated · #link(publication-cache.api_url)[
        Open the live INSPIRE query
      ]
    ]
  ]
  #publication-filters(
    publication-cache.authors,
    years,
    publication-types,
    publication-type-label,
    publication-cache.publications.len(),
  )
  #publication-list[
    #for publication in publication-cache.publications { publication-entry(publication) }
  ]
  #publication-no-script[All records are shown above. Enable JavaScript only if you want to filter or sort them in the browser.]
]

#let publications-page = portal-document(
  "Publications",
  "Filterable publications by verified αLoop contributors, sourced from INSPIRE HEP.",
  publications,
  active: "publications",
)
