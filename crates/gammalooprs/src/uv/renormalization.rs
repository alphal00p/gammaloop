use std::{fmt, ops::Deref};

use color_eyre::Result;
use idenso::tensor::SymbolicTensor;
use idenso::{CookMode, CookSettings};
use spenso::structure::IndexedStructure;
use symbolica::atom::Atom;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RenormalizationStats {
    pub forest_node_count: usize,
}

#[derive(Clone, Debug)]
pub struct RenormalizationPart {
    pub(crate) expression: Atom,
    pub stats: RenormalizationStats,
    // The public expression uses uncooked labels; its producer retains the
    // corresponding cooked interface for subsequent tensor operations.
    structure: IndexedStructure,
}

impl RenormalizationPart {
    pub fn expression(&self) -> &Atom {
        &self.expression
    }

    pub fn into_expression(self) -> Atom {
        self.expression
    }

    pub(crate) fn legacy(expression: Atom, structure: IndexedStructure) -> Self {
        Self {
            expression,
            stats: RenormalizationStats::default(),
            structure,
        }
    }

    pub(crate) fn new(
        expression: Atom,
        structure: IndexedStructure,
        forest_node_count: usize,
    ) -> Self {
        Self {
            expression,
            stats: RenormalizationStats { forest_node_count },
            structure,
        }
    }

    /// Recover the producer's tensor boundary without inferring the assembled sum.
    pub fn cooked_tensor(&self) -> Result<SymbolicTensor<IndexedStructure>> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        Ok(SymbolicTensor::from_produced_indexed_parts(
            cooking.try_cook(self.expression.as_view())?,
            self.structure.clone(),
        )?)
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
        self.expression == other.expression && self.structure == other.structure
    }
}

impl Eq for RenormalizationPart {}
