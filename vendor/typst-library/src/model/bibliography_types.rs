// Modified for optional notebook compiler features; see ../typst.typ.
use crate::foundations::{Content, elem};
use crate::introspection::Location;

/// The rendered parts for a bibliography.
pub struct RenderedBibliography {
    /// Lists all entries in the bibliography, with optional prefix.
    pub entries: Vec<RenderedEntry>,
    /// Whether the bibliography should have hanging indent applied.
    pub hanging_indent: bool,
}

/// The rendered parts for a bibliography entry.
pub struct RenderedEntry {
    /// An optional prefix. This is exposed separately because this will go into
    /// its own column for grid-based styles.
    pub prefix: Option<Content>,
    /// The main content of the rendered bibliography entry.
    pub body: Content,
    /// A location that should be attached to the rendered entry in some way.
    /// Citations will link there.
    pub backlink: Location,
}

/// Translation of `font-weight="light"` in CSL.
///
/// We translate `font-weight: "bold"` to `<strong>` since it's likely that the
/// CSL spec just talks about bold because it has no notion of semantic
/// elements. The benefits of a strict reading of the spec are also rather
/// questionable, while using semantic elements makes the bibliography more
/// accessible, easier to style, and more portable across export targets.
#[elem]
pub struct CslLightElem {
    #[required]
    pub body: Content,
}

/// Translation of `display="indent"` in CSL.
///
/// A `display="block"` is simply translated to a Typst `BlockElem`. Similarly,
/// we could translate `display="indent"` to a `PadElem`, but (a) it does not
/// yet have support in HTML and (b) a `PadElem` described a fixed padding while
/// CSL leaves the amount of padding user-defined so it's not a perfect fit.
#[elem]
pub struct CslIndentElem {
    #[required]
    pub body: Content,
}
