// Modified for optional notebook compiler features; see ../typst.typ.
//! PDF images cannot be constructed when PDF loading is disabled.
/// An unavailable PDF image.
#[derive(Clone, Hash)]
pub enum PdfImage {}
impl PdfImage {
    /// The image width.
    pub fn width(&self) -> f32 {
        match *self {}
    }
    /// The image height.
    pub fn height(&self) -> f32 {
        match *self {}
    }
}
