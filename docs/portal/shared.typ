// Shared portal data and route selection. Presentation lives in components.typ.

#let portal = toml("/docs/portal.toml")
#let registry = toml("/docs/products/registry.toml")
#let talk-catalog = toml("/docs/talks.toml")
#let publication-cache = json("/docs/data/publications.json")

#let docs-channel = sys.inputs.at("channel", default: "latest")
#let docs-snapshot-tag = sys.inputs.at("snapshot-tag", default: "")
#let channel-route = if docs-channel == "snapshot" {
  "snapshots/" + docs-snapshot-tag
} else {
  "latest"
}
