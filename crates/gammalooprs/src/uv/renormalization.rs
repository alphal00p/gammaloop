use std::{fmt, ops::Deref};

use symbolica::atom::Atom;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RenormalizationStats {
    pub kernel_hits: usize,
    pub forest_node_count: usize,
}

#[derive(Clone, Debug)]
pub struct RenormalizationPart {
    pub expression: Atom,
    pub stats: RenormalizationStats,
}

impl RenormalizationPart {
    pub(crate) fn new(expression: Atom, kernel_hits: usize, forest_node_count: usize) -> Self {
        Self {
            expression,
            stats: RenormalizationStats {
                kernel_hits,
                forest_node_count,
            },
        }
    }
}

impl Deref for RenormalizationPart {
    type Target = Atom;

    fn deref(&self) -> &Self::Target {
        &self.expression
    }
}

impl fmt::Display for RenormalizationPart {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.expression.fmt(f)
    }
}

impl PartialEq for RenormalizationPart {
    fn eq(&self, other: &Self) -> bool {
        self.expression == other.expression
    }
}

impl Eq for RenormalizationPart {}
