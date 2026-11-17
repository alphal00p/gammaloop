//! Collector factor selection: hand the selected-tensor collector only the
//! factors it has to split.
//!
//! The collector (`coefficient_list_with`) observes and parses every node of
//! its input, including scalar factors (dots, couplings, scalar sums) that
//! only end up in the coefficients. [`HiddenScalars::hide`] replaces, in every
//! product outside function arguments, the product of the factors that carry
//! no slot and that the collector would not select by one opaque placeholder
//! `hidden_scalar(n)`. Factors the collector selects (a propagator `denom(...)`
//! on the block path, every Minkowski dot on the monomial path) and
//! tensor-carrying factors stay where they are, so the collector sees the same
//! selected leaves in the same products and returns the same rows.
//! [`HiddenScalars::restore`] substitutes the placeholders back into each row.
//! This is a substitution of a scalar by its value, so every row is exact.

use std::collections::HashMap;
use std::sync::LazyLock;

use idenso::tensor::SymbolicTensor;
use spenso::structure::partial::PartialStructure;
use spenso::structure::slot::{SlotMatch, SlotMatcher};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::symbol;
use symbolica::utils::resources;

use super::TensorReductionError;

static HIDDEN: LazyLock<Symbol> = LazyLock::new(|| symbol!("feynkit_tensor::hidden_scalar"));

/// Collector rows: the selected structure and its coefficient.
pub(super) type Rows = Vec<(
    SymbolicTensor<PartialStructure>,
    SymbolicTensor<PartialStructure>,
)>;

/// What a subtree holds, as far as hiding is concerned.
#[derive(Clone, Copy, Default)]
struct Flags {
    /// An explicit slot anywhere, or a function with a compact representation
    /// argument (an anonymous port).
    slot: bool,
    /// The selector accepts the node or a function outside function arguments.
    selected: bool,
}

impl Flags {
    fn or(self, other: Flags) -> Flags {
        Flags {
            slot: self.slot || other.slot,
            selected: self.selected || other.selected,
        }
    }
}

/// The hidden scalar products, by placeholder number.
pub(super) struct HiddenScalars {
    values: Vec<Atom>,
}

struct Hider<'a> {
    select: &'a mut dyn FnMut(AtomView<'_>) -> Result<bool, TensorReductionError>,
    matcher: SlotMatcher,
    index: HashMap<Atom, usize>,
    values: Vec<Atom>,
    /// The flags of every product factor, by [`node_key`].
    factors: HashMap<(usize, usize), Flags>,
}

impl Hider<'_> {
    /// A function factor: never descended into for hiding, scanned once.
    fn function(&mut self, value: AtomView<'_>) -> Result<Flags, TensorReductionError> {
        let AtomView::Fun(function) = value else {
            unreachable!()
        };
        let mut flags = Flags {
            slot: false,
            selected: (self.select)(value)?,
        };
        for argument in function.iter() {
            if self.matcher.compact_representation(argument).is_some() {
                flags.slot = true;
            }
        }
        if !flags.slot {
            let matcher = &mut self.matcher;
            let mut slot = false;
            value.visitor(&mut |node| {
                if matches!(node, AtomView::Fun(_))
                    && matches!(matcher.classify(node), SlotMatch::Explicit(_))
                {
                    slot = true;
                }
                !slot
            });
            flags.slot = slot;
        }
        Ok(flags)
    }

    /// The flags of `value`. Records the flags of every product factor so that
    /// [`Hider::rebuild`] decides which factors to hide without walking them.
    fn flags(&mut self, value: AtomView<'_>) -> Result<Flags, TensorReductionError> {
        resources::checkpoint()?;
        Ok(match value {
            AtomView::Num(_) => Flags::default(),
            AtomView::Var(_) => Flags {
                slot: false,
                selected: (self.select)(value)?,
            },
            AtomView::Fun(_) => self.function(value)?,
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let mut flags = self.flags(base)?.or(self.flags(exponent)?);
                flags.selected |= (self.select)(value)?;
                flags
            }
            AtomView::Add(sum) => {
                let mut flags = Flags::default();
                for term in sum.iter() {
                    flags = flags.or(self.flags(term)?);
                }
                flags.selected |= (self.select)(value)?;
                flags
            }
            AtomView::Mul(product) => {
                let mut flags = Flags::default();
                for factor in product.iter() {
                    let factor_flags = self.flags(factor)?;
                    self.factors.insert(node_key(factor), factor_flags);
                    flags = flags.or(factor_flags);
                }
                flags.selected |= (self.select)(value)?;
                flags
            }
        })
    }

    /// Rebuild `value` with hidden scalar products; `None` if unchanged.
    /// Needs the factor flags recorded by [`Hider::flags`] on `value`; a factor
    /// hidden whole is neither walked nor rebuilt.
    fn rebuild(&mut self, value: AtomView<'_>) -> Result<Option<Atom>, TensorReductionError> {
        resources::checkpoint()?;
        Ok(match value {
            AtomView::Num(_) | AtomView::Var(_) | AtomView::Fun(_) => None,
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                self.rebuild(base)?
                    .map(|base| base.pow(exponent.to_owned()))
            }
            AtomView::Add(sum) => {
                let mut summands = Vec::new();
                let mut changed = false;
                for term in sum.iter() {
                    let rebuilt = self.rebuild(term)?;
                    changed |= rebuilt.is_some();
                    summands.push(rebuilt.unwrap_or_else(|| term.to_owned()));
                }
                changed.then(|| Atom::add_many(summands.iter()))
            }
            AtomView::Mul(product) => {
                let mut kept = Vec::new();
                let mut hidden = Vec::new();
                let mut changed = false;
                for factor in product.iter() {
                    let flags = self.factors[&node_key(factor)];
                    if !flags.slot && !flags.selected && !matches!(factor, AtomView::Num(_)) {
                        hidden.push(factor);
                        continue;
                    }
                    let rebuilt = self.rebuild(factor)?;
                    changed |= rebuilt.is_some();
                    kept.push(rebuilt.unwrap_or_else(|| factor.to_owned()));
                }
                // One hidden variable alone gains nothing.
                if hidden.len() > 1
                    || hidden
                        .first()
                        .is_some_and(|f| !matches!(f, AtomView::Var(_)))
                {
                    let scalar = Atom::mul_many(hidden.iter().copied());
                    let next = self.values.len();
                    let number = *self.index.entry(scalar.clone()).or_insert(next);
                    if number == next {
                        self.values.push(scalar);
                    }
                    kept.push(FunctionBuilder::new(*HIDDEN).add_arg(number).finish());
                    changed = true;
                } else {
                    kept.extend(hidden.into_iter().map(|factor| factor.to_owned()));
                }
                changed.then(|| Atom::mul_many(kept.iter()))
            }
        })
    }
}

/// A node by its position in the walked expression: product factors borrow
/// their bytes from the expression, so distinct factors have distinct keys.
fn node_key(value: AtomView<'_>) -> (usize, usize) {
    let data = value.get_data();
    (data.as_ptr() as usize, data.len())
}

impl HiddenScalars {
    /// Hide the scalar factors of every product of `expression` that carry no
    /// slot and are not selected by `select`. `None` if nothing is hidden or
    /// the expression already uses the placeholder.
    pub(super) fn hide(
        expression: AtomView<'_>,
        select: &mut dyn FnMut(AtomView<'_>) -> Result<bool, TensorReductionError>,
    ) -> Result<Option<(Atom, HiddenScalars)>, TensorReductionError> {
        let mut hider = Hider {
            select,
            matcher: SlotMatcher::default(),
            index: HashMap::new(),
            values: Vec::new(),
            factors: HashMap::new(),
        };
        if expression.contains_symbol(*HIDDEN) {
            return Ok(None);
        }
        hider.flags(expression)?;
        let rebuilt = hider.rebuild(expression)?;
        if hider.values.is_empty() {
            return Ok(None);
        }
        Ok(rebuilt.map(|rebuilt| {
            (
                rebuilt,
                HiddenScalars {
                    values: hider.values,
                },
            )
        }))
    }

    /// Whether `value` is a placeholder.
    pub(super) fn is_placeholder(value: AtomView<'_>) -> bool {
        matches!(value, AtomView::Fun(function) if function.get_symbol() == *HIDDEN)
    }

    /// Substitute the hidden products back into an expression.
    pub(super) fn substitute(&self, expression: AtomView<'_>) -> Option<Atom> {
        let hidden = *HIDDEN;
        if !expression.contains_symbol(hidden) {
            return None;
        }
        Some(expression.replace_map(|value, _, out| {
            if let AtomView::Fun(function) = value
                && function.get_symbol() == hidden
                && function.get_nargs() == 1
                && let Ok(number) = i64::try_from(function.iter().next().unwrap())
                && let Some(scalar) = usize::try_from(number)
                    .ok()
                    .and_then(|number| self.values.get(number))
            {
                **out = scalar.clone();
            }
        }))
    }

    /// Substitute the hidden products back into every row. The pieces keep
    /// their interfaces: a placeholder and its value are both scalars.
    pub(super) fn restore(&self, rows: Rows) -> Result<Rows, TensorReductionError> {
        rows.into_iter()
            .map(|(selected, coefficient)| {
                resources::checkpoint()?;
                let restore = |piece: SymbolicTensor<PartialStructure>| match self
                    .substitute(piece.expression().as_view())
                {
                    Some(expression) => piece.with_algebra_result(expression),
                    None => Ok(piece),
                };
                Ok((restore(selected)?, restore(coefficient)?))
            })
            .collect()
    }

    /// Number of distinct hidden products.
    #[cfg(test)]
    pub(super) fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether `expression` uses every placeholder of the table.
    #[cfg(test)]
    pub(super) fn all_used_in(&self, expression: AtomView<'_>) -> bool {
        (0..self.values.len()).all(|number| {
            let placeholder = FunctionBuilder::new(*HIDDEN).add_arg(number).finish();
            let mut found = false;
            expression.visitor(&mut |node| {
                found |= node == placeholder.as_view();
                !found
            });
            found
        })
    }
}

/// Collect `source` with `select`, hiding unselected scalar factors from the
/// collector first (see the module documentation). Falls back to collecting
/// `source` itself if nothing can be hidden or the hidden expression does not
/// keep the source interface. The flag tells whether anything was hidden.
pub(super) fn collect_rows(
    source: &SymbolicTensor<PartialStructure>,
    hide: bool,
    select: &mut dyn FnMut(AtomView<'_>) -> Result<bool, TensorReductionError>,
) -> Result<(Rows, bool), TensorReductionError> {
    resources::checkpoint()?;
    if hide
        && let Some((hidden, scalars)) = HiddenScalars::hide(source.expression().as_view(), select)?
        && let Ok(hidden) = source.with_checked_expression(hidden)
    {
        let rows = collect(&hidden, &mut |value: AtomView<'_>| {
            if HiddenScalars::is_placeholder(value) {
                return Ok(false);
            }
            select(value)
        })?;
        return Ok((scalars.restore(rows)?, true));
    }
    Ok((collect(source, select)?, false))
}

fn collect(
    source: &SymbolicTensor<PartialStructure>,
    select: &mut dyn FnMut(AtomView<'_>) -> Result<bool, TensorReductionError>,
) -> Result<Rows, TensorReductionError> {
    let mut error = None;
    let rows = source
        .coefficient_list_with(|value| match select(value) {
            Ok(selected) => selected,
            Err(failure) => {
                error.get_or_insert(failure);
                false
            }
        })
        .map_err(|e| TensorReductionError::TensorAlgebra(e.to_string()))?;
    resources::checkpoint()?;
    if let Some(error) = error {
        return Err(error);
    }
    Ok(rows)
}
