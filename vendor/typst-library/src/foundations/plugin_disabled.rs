// Modified for optional notebook compiler features; see ../typst.typ.
//! Without a plugin loader, plugin functions cannot be constructed.
use crate::diag::StrResult;
use crate::foundations::Bytes;
use ecow::EcoString;

/// An unavailable WebAssembly plugin function.
#[derive(Debug, Clone, PartialEq, Hash)]
pub enum PluginFunc {}

impl PluginFunc {
    /// The name of the plugin function.
    pub fn name(&self) -> &EcoString {
        match *self {}
    }
    /// Call the plugin function.
    pub fn call(&self, _args: Vec<Bytes>) -> StrResult<Bytes> {
        match *self {}
    }
}
