use std::sync::LazyLock;

use idenso::representations::Bispinor;
use spenso::structure::representation::{LibraryRep, Minkowski};

pub struct GammaloopRepresentations {
    pub mink: LibraryRep,
    pub bis: LibraryRep,
}

pub static GR: LazyLock<GammaloopRepresentations> = LazyLock::new(|| GammaloopRepresentations {
    mink: Minkowski {}.into(),
    bis: Bispinor {}.into(),
});
