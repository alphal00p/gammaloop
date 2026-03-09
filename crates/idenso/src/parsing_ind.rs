use spenso::structure::slot::{AbsInd, DummyAind, ParseableAind, SlotError};
use std::sync::atomic::AtomicUsize;
use symbolica::{
    atom::{Atom, AtomView},
    function, symbol,
};

static DUMMYCOUNTER: AtomicUsize = AtomicUsize::new(0);
/// A type that represents the name of an index in a tensor.
#[derive(Debug, Copy, Clone, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Parsind {
    Dummy(usize),
}

impl AbsInd for Parsind {}

impl DummyAind for Parsind {
    fn new_dummy() -> Self {
        Parsind::Dummy(DUMMYCOUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed))
    }

    fn new_dummy_at(i: usize) -> Self {
        Parsind::Dummy(i)
    }
    fn is_dummy(&self) -> bool {
        true
    }
}

impl ParseableAind for Parsind {
    type Error = SlotError;
    fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
        view.try_into()
    }
    fn to_atom(&self) -> Atom {
        (*self).into()
    }
}

impl std::fmt::Display for Parsind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "parsind")
    }
}

impl From<usize> for Parsind {
    fn from(value: usize) -> Self {
        Parsind::Dummy(value)
    }
}

impl From<Parsind> for Atom {
    fn from(value: Parsind) -> Self {
        match value {
            Parsind::Dummy(i) => function!(symbol!("dummy"), i as i64),
        }
    }
}
impl From<Parsind> for symbolica::atom::AtomOrView<'_> {
    fn from(value: Parsind) -> Self {
        symbolica::atom::AtomOrView::Atom(Atom::from(value))
    }
}

impl TryFrom<AtomView<'_>> for Parsind {
    type Error = SlotError;

    fn try_from(_: AtomView<'_>) -> Result<Self, Self::Error> {
        Ok(Parsind::Dummy(0))
    }
}
