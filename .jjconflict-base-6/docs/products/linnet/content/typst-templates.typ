#import "/crates/linnest/typst/docs/manual.typ": domain-style-concepts
#import "../../shared.typ": source-link

#let typst-templates = [
= Domain styles and templates

Linnest provides generic graph styling and drawing callbacks. Application domains
compose those primitives in ordinary Typst code or provide a custom V1 template;
there is no built-in physics module.

GammaLoop owns its particle-line, momentum, amplitude, and cross-section template.
Those policies are available when rendering through GammaLoop's template bundle,
not by importing them from Linnest. See GammaLoop's
#source-link(
  "assets/embedded/drawing/templates/physics-edge-style.typ",
  label: "physics drawing template",
).

== Concepts

#domain-style-concepts
]
