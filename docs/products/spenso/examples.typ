#import "../shared.typ": source-link

#let live-example(title, summary, notebook) = [
= #title

#summary

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": notebook,
    "aria-label": title,
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      Open the #source-link("examples/notebooks/" + notebook + ".py", label: "Marimo source") to run this example locally.
    ]
  ]
} else {
  [Run the #source-link("examples/notebooks/" + notebook + ".py", label: "Marimo source") locally.]
}
]

