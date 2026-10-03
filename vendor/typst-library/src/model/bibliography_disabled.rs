// Modified for optional notebook compiler features; see ../typst.typ.
//! The notebook compiler has no citation database or CSL processor.
//! Shared realization rules still mention these elements, but the global
//! library does not expose constructors and no rendered works can exist.
use super::RenderedBibliography;
use crate::diag::{SourceResult, bail};
use crate::engine::Engine;
use crate::foundations::{Content, Derived, Label, Packed, StyleChain, elem};
use crate::introspection::{Locatable, Location};
use crate::loading::DataSource;
use std::sync::Arc;
use typst_syntax::{Span, Spanned};

#[elem(Locatable)]
pub struct BibliographyElem {}

impl BibliographyElem {
    pub fn has(_engine: &mut Engine, _key: Label, _span: Span) -> bool {
        false
    }
}

impl Packed<BibliographyElem> {
    pub fn realize_title(&self, _styles: StyleChain) -> Option<Content> {
        None
    }
}

pub type CslSource = DataSource;

#[derive(Debug, Clone, PartialEq, Hash)]
pub enum CslStyle {}

impl CslStyle {
    pub fn load(
        _engine: &mut Engine,
        source: Spanned<CslSource>,
    ) -> SourceResult<Derived<CslSource, Self>> {
        bail!(source.span, "bibliographies are not enabled")
    }
}

pub enum Works {}

impl Works {
    pub fn generate(_engine: &mut Engine, span: Span) -> SourceResult<Arc<Self>> {
        bail!(span, "bibliographies are not enabled")
    }
    pub fn citation(&self, _loc: Location, _span: Span) -> SourceResult<Content> {
        match *self {}
    }
    pub fn bibliography(
        &self,
        _loc: Location,
        _span: Span,
    ) -> SourceResult<&RenderedBibliography> {
        match *self {}
    }
}
