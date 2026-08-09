//! Exact Young-projector terms for carriers with structural column antisymmetry.
//!
//! Under the pullback convention used by tensor arguments, `C_T R_T` produces
//! selectors `column.compose(row)`. Once every tableau column is represented by
//! a structural antisymmetry group, selectors that differ by the right action
//! of the column group represent the same carrier up to that action's sign.
//! Equal-height columns may also exchange as unsigned whole blocks because
//! their pointwise exchange lies in the row group. Occurrence-local reduction
//! therefore uses the signed carrier stabilizer `C ⋊ W` without
//! interchanging `C_T` and `R_T`, which generally do not commute.
//!
//! Expressions without a Young-containing Power, a custom normalizer on any
//! exposed LocalTensor head or Function in the complete root, or repeated
//! occurrences of one general-Young head take a factored path. Each declared
//! tensor is replaced by its exact reduced projector after occurrence-local
//! column-class reduction, while Products keep distinct projector sums
//! factored. The executed projected root is parsed with the canonical policy,
//! including tensor-symmetry and Power-grammar validation.
//!
//! A lone general-Young tensor whose index lines are all distinct and external
//! sends the declared projected sum directly through the ordinary graph driver,
//! then selects the shared numeric normal form. It needs no carrier to align
//! declared heads or dummy frames. Other eligible factored expressions emit
//! projector terms through one deterministic private `Linear` carrier. Each
//! declared head is an opaque carrier argument; its tableau columns are direct
//! slots or ordered slot bundles.
//! Strict private metadata promotes those bundles back to graph-owned Young
//! columns, so the ordinary whole-root graph fixed-point driver canonizes the
//! carrier network with signed internal-column order and unsigned exchange of
//! equal-height columns. Decoding restores the original heads, collects root
//! numeric content, and chooses a deterministic orientation for every primitive
//! Add factor. It does not apply the projector or run the graph again.
//! The rebuilt Atom is fully validated and canonically parsed into a
//! graph-canonical `CanonicalPolicyNet`.
//!
//! Only the successful lone-root direct route returns a terminal Atom to the
//! Atom-facing canonicalization entry point without reparsing it. Carrier,
//! composite, and staged routes retain `CanonicalPolicyNet`; the test-only
//! policy-returning API reparses the direct terminal Atom when a network/Atom
//! pair is specifically required.
//!
//! A carrier `ConvergenceCycle` retries as an exact post-projector composite
//! loop: one ordinary graph step is followed by the same factored projector,
//! consecutive projected Atoms certify stability, and the middle graph result
//! is returned. A graph-size limit instead falls back to the staged fold.
//!
//! The staged fold distributes only projector alternatives. Each intermediate
//! candidate is executed, reparsed, and canonized as a complete current graph;
//! unrelated sums and nonlinear boundaries remain factored. Root candidates
//! use the caller's dummy allocator and are aggregated without a final
//! whole-Sum graph call. Candidate graphs retain the ordinary signed fixed-point
//! and graph-budget contracts, while projector actions and logical live terms
//! have separate checked limits. A Young-containing Power enters this path
//! directly. So does a source with a custom normalizer on any exposed
//! LocalTensor head or Function in the complete Young-containing root, including
//! a sibling outside the Young subtree, or one with repeated occurrences of one
//! general Young head. The staged result is an exact public section on a fresh
//! call and avoids passing the graph-rebuilt root back through user normalizers
//! during carrier decode and canonical reparse.
//!
//! A retained even Power gives its nondistributable base one anonymous boundary
//! frame before its copies are materialized. Concrete boundary names and the
//! base's overall sign are gauge there: the even magnitude binds every boundary
//! and cancels a uniform sign. This prevents the ordinary complete-Power pass
//! from depending on the temporary frame chosen for the projected base.

#[cfg(test)]
use std::sync::{LazyLock, Mutex};
use std::{
    cell::{OnceCell, RefCell},
    collections::{BTreeMap, BTreeSet},
};

use itertools::Itertools;
use linnet::{
    half_edge::{NodeIndex, tree::SimpleTraversalTree},
    permutation::Permutation,
    tree::child_pointer::ParentChildStore,
};
use spenso::{
    contraction::Contract,
    network::{
        Network,
        graph::{NetworkLeaf, NetworkNode, NetworkOp},
        parsing::{AtomStructureExt, StrictTensorFilter},
        store::TensorScalarStore,
    },
    structure::{
        OrderedStructure, YoungTableau, YoungTableauClass,
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind, Slot},
    },
    tensor_symbol,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    domains::{integer::Integer, rational::Rational},
};
use thiserror::Error;

use super::{
    CanonicalPolicyNet, CanonicalizationError,
    driver::{CanonicalRequest, execute_atom},
    projection::ExternalLineMode,
    semantic::SemanticAtomKey,
};
use crate::tensor::{SymbolicNet, SymbolicTensor};

const MAX_YOUNG_PROJECTOR_ACTIONS: usize = 4_096;
const MAX_YOUNG_EXPANSION_TERMS: usize = 4_096;

/// Private graph carrier for an already-projected general Young tensor.
///
/// Its first, opaque argument stores the declared tensor head. The remaining
/// arguments are the tableau columns, represented by direct slots for
/// singletons and ordered `aind(...)` bundles otherwise. Strict carrier
/// metadata, rather than the public syntax, makes those bundles graph-owned
/// Young columns. Every declaration shares one private head, so differently
/// declared tensors compare through the opaque original-head color in the same
/// semantic order as their decoded public heads. `Linear` makes a column sign
/// lift to the complete tensor factor.
pub(super) struct YoungColumnCarrier;

impl YoungColumnCarrier {
    const SYMBOL_NAMESPACE_PREFIX: &str = "idenso::young_column_carrier_";
    const SYMBOL_NAME: &str = "idenso::young_column_carrier_v1";

    fn symbol() -> Symbol {
        tensor_symbol!(Self::SYMBOL_NAME; Linear)
    }

    /// Recover the declaration carried by one generated private head.
    ///
    /// The reserved namespace is strict: malformed or unsupported internal
    /// carriers are errors rather than ordinary tensors. The opaque original
    /// head remains the authoritative identity, tableau payload, and graph
    /// color discriminator; the fixed generated head makes its semantic order
    /// agree with the decoded declared-head order.
    pub(super) fn declaration(
        function: symbolica::atom::representation::FunView<'_>,
    ) -> Result<Option<(Symbol, YoungTableau)>, CanonicalizationError> {
        let carrier = function.get_symbol();
        let name = carrier.get_name();
        if !name.starts_with(Self::SYMBOL_NAMESPACE_PREFIX) {
            return Ok(None);
        }
        if carrier != Self::symbol() {
            return Err(CanonicalizationError::Projection(format!(
                "unsupported internal Young-column carrier `{name}`"
            )));
        }
        let Some(AtomView::Var(original)) = function.iter().next() else {
            return Err(CanonicalizationError::Projection(
                "internal Young-column carrier lost its declared tensor head".to_owned(),
            ));
        };
        let original = original.get_symbol();
        let tableau = YoungTableau::from_symbol(original)
            .map_err(|error| CanonicalizationError::InvalidYoungTableauMetadata {
                head: original,
                reason: error.to_string(),
            })?
            .filter(|tableau| tableau.class() == YoungTableauClass::General)
            .ok_or_else(|| {
                CanonicalizationError::Projection(
                    "internal Young-column carrier references a non-general tableau".to_owned(),
                )
            })?;
        let columns = tableau.columns().len();
        if function.get_nargs() != columns + 1 {
            return Err(CanonicalizationError::Projection(format!(
                "internal Young-column carrier for {original} has {} column arguments, expected {columns}",
                function.get_nargs().saturating_sub(1),
            )));
        }
        Ok(Some((original, tableau)))
    }

    fn decode_policy<Aind>(
        policy: CanonicalPolicyNet<Aind>,
    ) -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
    {
        CanonicalPolicyNet::parse(Self::decode::<Aind>(policy.into_atom())?)
    }

    #[cfg(test)]
    pub(super) fn tensor<Aind>(
        projected: SymbolicTensor<Aind>,
        tableau: &YoungTableau,
    ) -> Result<SymbolicTensor<Aind>, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
    {
        let columns = tableau.columns().collect::<Vec<_>>();
        Self::tensor_with_symbol(projected, tableau, &columns, Self::symbol())
    }

    /// Replace one normalized projected tensor by its ordered-column carrier.
    ///
    /// Occurrence-local `C ⋊ W` normalization has already chosen each column's
    /// order and lifted its parity into the projector coefficient. Keeping the
    /// payload ordered avoids a second Symbolica structural-group normalization;
    /// strict layout metadata restores the signed column sites in Graphica.
    fn tensor_with_symbol<Aind>(
        projected: SymbolicTensor<Aind>,
        tableau: &YoungTableau,
        columns: &[Vec<usize>],
        carrier_symbol: Symbol,
    ) -> Result<SymbolicTensor<Aind>, CanonicalizationError>
    where
        Aind: AbsInd + DummyAind + ParseableAind,
    {
        let expression = projected.expression.clone();
        let AtomView::Fun(function) = expression.as_view() else {
            return Err(CanonicalizationError::StructureMismatch { expression });
        };
        let head = function.get_symbol();
        let arguments = function
            .iter()
            .map(|argument| argument.to_owned())
            .collect::<Vec<_>>();
        if arguments.len() != tableau.rank() {
            return Err(CanonicalizationError::InvalidYoungTableauArity {
                head,
                expected: tableau.rank(),
                actual: arguments.len(),
                expression,
            });
        }

        let mut carrier = FunctionBuilder::new(carrier_symbol).add_arg(Atom::var(head));
        let mut carrier_slots = Vec::with_capacity(tableau.rank());
        for column in columns {
            if let [position] = column.as_slice() {
                carrier = carrier.add_arg(arguments[*position].clone());
                carrier_slots.push(
                    Slot::<LibraryRep, Aind>::try_from(arguments[*position].as_view()).map_err(
                        |_| CanonicalizationError::StructureMismatch {
                            expression: expression.clone(),
                        },
                    )?,
                );
                continue;
            }

            let mut bundle = FunctionBuilder::new(AIND_SYMBOLS.aind);
            for position in column {
                let argument = &arguments[*position];
                carrier_slots.push(
                    Slot::<LibraryRep, Aind>::try_from(argument.as_view()).map_err(|_| {
                        CanonicalizationError::StructureMismatch {
                            expression: expression.clone(),
                        }
                    })?,
                );
                bundle = bundle.add_arg(argument.clone());
            }
            carrier = carrier.add_arg(bundle.finish());
        }

        let expression = carrier.finish();
        if !matches!(expression.as_view(), AtomView::Fun(_)) {
            return Err(CanonicalizationError::StructureMismatch { expression });
        }
        Ok(SymbolicTensor {
            structure: OrderedStructure::new(carrier_slots).structure,
            is_metric: false,
            is_composite: false,
            expression,
        })
    }

    fn decode<Aind>(value: Atom) -> Result<Atom, CanonicalizationError>
    where
        Aind: AbsInd + ParseableAind,
    {
        Ok(Self::normalize_numeric_content(Self::decode_view::<Aind>(
            value.as_view(),
        )?))
    }

    fn decode_view<Aind>(value: AtomView<'_>) -> Result<Atom, CanonicalizationError>
    where
        Aind: AbsInd + ParseableAind,
    {
        match value {
            AtomView::Fun(function) => {
                if let Some((head, tableau)) = Self::declaration(function)? {
                    return Self::decode_carrier::<Aind>(value, function, head, &tableau);
                }
                if !value.contains_exposed_tensor_topology(StrictTensorFilter::Tagged) {
                    return Ok(value.to_owned());
                }
                let mut rebuilt = FunctionBuilder::new(function.get_symbol());
                for argument in function.iter() {
                    rebuilt = rebuilt.add_arg(Self::decode_view::<Aind>(argument)?);
                }
                Ok(rebuilt.finish())
            }
            AtomView::Add(sum) => sum.iter().try_fold(Atom::Zero, |result, term| {
                Ok(result + Self::decode_view::<Aind>(term)?)
            }),
            AtomView::Mul(product) => product.iter().try_fold(Atom::num(1), |result, factor| {
                Ok(result * Self::decode_view::<Aind>(factor)?)
            }),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Ok(Self::decode_view::<Aind>(base)?.pow(Self::decode_view::<Aind>(exponent)?))
            }
            AtomView::Num(_) | AtomView::Var(_) => Ok(value.to_owned()),
        }
    }

    /// Hoist the content of every projected sum to one scalar, then choose a
    /// process-independent orientation for each primitive Add factor.
    fn normalize_numeric_content(value: Atom) -> Atom {
        let collected = value.collect_num();
        let mut coefficient = Atom::num(1);
        let factors = match collected.as_view() {
            AtomView::Mul(product) => product.iter().collect_vec(),
            _ => vec![collected.as_view()],
        };
        let mut primitive = Vec::with_capacity(factors.len());
        for factor in factors {
            let AtomView::Add(terms) = factor else {
                if matches!(factor, AtomView::Num(_)) {
                    coefficient *= factor.to_owned();
                } else {
                    primitive.push(factor.to_owned());
                }
                continue;
            };
            let positive = factor.to_owned();
            let negative = terms
                .iter()
                .fold(Atom::Zero, |sum, term| sum - term.to_owned());
            if SemanticAtomKey::new(negative.as_view()) < SemanticAtomKey::new(positive.as_view()) {
                coefficient = -coefficient;
                primitive.push(negative);
            } else {
                primitive.push(positive);
            }
        }
        coefficient
            * primitive
                .into_iter()
                .fold(Atom::num(1), |product, factor| product * factor)
    }

    fn decode_carrier<Aind>(
        value: AtomView<'_>,
        function: symbolica::atom::representation::FunView<'_>,
        head: Symbol,
        tableau: &YoungTableau,
    ) -> Result<Atom, CanonicalizationError>
    where
        Aind: AbsInd + ParseableAind,
    {
        let expression = value.to_owned();
        let mut carrier_arguments = function.iter();
        let _original = carrier_arguments
            .next()
            .expect("a validated Young-column carrier stores its declared tensor head");
        let columns = tableau.columns().collect::<Vec<_>>();

        let mut arguments = vec![None; tableau.rank()];
        for (column, encoded) in columns.iter().zip(carrier_arguments) {
            if let [position] = column.as_slice() {
                let slot = Slot::<LibraryRep, Aind>::try_from(encoded).map_err(|_| {
                    CanonicalizationError::StructureMismatch {
                        expression: expression.clone(),
                    }
                })?;
                arguments[*position] = Some(slot.to_atom());
                continue;
            }

            let AtomView::Fun(bundle) = encoded else {
                return Err(CanonicalizationError::StructureMismatch {
                    expression: expression.clone(),
                });
            };
            if bundle.get_symbol() != AIND_SYMBOLS.aind || bundle.get_nargs() != column.len() {
                return Err(CanonicalizationError::StructureMismatch {
                    expression: expression.clone(),
                });
            }
            for (&position, member) in column.iter().zip(bundle.iter()) {
                let slot = Slot::<LibraryRep, Aind>::try_from(member).map_err(|_| {
                    CanonicalizationError::StructureMismatch {
                        expression: expression.clone(),
                    }
                })?;
                arguments[position] = Some(slot.to_atom());
            }
        }
        let decoded = arguments
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| CanonicalizationError::StructureMismatch {
                expression: expression.clone(),
            })?
            .into_iter()
            .fold(FunctionBuilder::new(head), FunctionBuilder::add_arg)
            .finish();
        Ok(decoded)
    }
}

#[cfg(test)]
static FACTORED_GRAPH_OUTCOMES: LazyLock<Mutex<BTreeMap<Symbol, (usize, usize)>>> =
    LazyLock::new(|| Mutex::new(BTreeMap::new()));

#[cfg(test)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct YoungRouteOutcomes {
    carrier: usize,
    composite: usize,
    staged: usize,
}

#[cfg(test)]
static YOUNG_ROUTE_OUTCOMES: LazyLock<Mutex<BTreeMap<Symbol, YoungRouteOutcomes>>> =
    LazyLock::new(|| Mutex::new(BTreeMap::new()));

#[cfg(test)]
fn reset_factored_graph_outcomes(head: Symbol) {
    FACTORED_GRAPH_OUTCOMES
        .lock()
        .expect("factored graph outcome lock is not poisoned")
        .insert(head, (0, 0));
}

#[cfg(test)]
fn factored_graph_outcomes(head: Symbol) -> (usize, usize) {
    FACTORED_GRAPH_OUTCOMES
        .lock()
        .expect("factored graph outcome lock is not poisoned")
        .get(&head)
        .copied()
        .unwrap_or_default()
}

#[cfg(test)]
fn reset_young_route_outcomes(head: Symbol) {
    YOUNG_ROUTE_OUTCOMES
        .lock()
        .expect("Young route outcome lock is not poisoned")
        .insert(head, YoungRouteOutcomes::default());
}

#[cfg(test)]
fn young_route_outcomes(head: Symbol) -> YoungRouteOutcomes {
    YOUNG_ROUTE_OUTCOMES
        .lock()
        .expect("Young route outcome lock is not poisoned")
        .get(&head)
        .copied()
        .unwrap_or_default()
}

#[cfg(test)]
#[derive(Clone, Copy)]
enum YoungRoute {
    Carrier,
    Composite,
    Staged,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum YoungProjectorQuantity {
    RowActionCount,
    ColumnActionCount,
    FullActionCount,
    HookProduct,
}

/// Checked sizes that callers can inspect before projector actions are allocated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct YoungProjectorPreflight {
    rank: usize,
    row_action_count: usize,
    column_action_count: usize,
    full_action_count: usize,
    hook_product: usize,
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub(super) enum YoungProjectorPlanError {
    #[error("Young-projector {quantity:?} overflows usize while multiplying {partial} by {factor}")]
    ArithmeticOverflow {
        quantity: YoungProjectorQuantity,
        partial: usize,
        factor: usize,
    },
    #[error("cannot reserve {requested} permutations for Young-projector {quantity:?}")]
    Allocation {
        quantity: YoungProjectorQuantity,
        requested: usize,
    },
}

impl YoungProjectorPreflight {
    pub(super) fn new(tableau: &YoungTableau) -> Result<Self, YoungProjectorPlanError> {
        let row_action_count = Self::subgroup_order(
            tableau.rows().map(<[usize]>::len),
            YoungProjectorQuantity::RowActionCount,
        )?;
        let column_action_count = Self::subgroup_order(
            tableau.columns().map(|column| column.len()),
            YoungProjectorQuantity::ColumnActionCount,
        )?;
        let full_action_count = row_action_count.checked_mul(column_action_count).ok_or(
            YoungProjectorPlanError::ArithmeticOverflow {
                quantity: YoungProjectorQuantity::FullActionCount,
                partial: row_action_count,
                factor: column_action_count,
            },
        )?;
        let hook_product = Self::checked_hook_product(tableau.shape())?;

        Ok(Self {
            rank: tableau.rank(),
            row_action_count,
            column_action_count,
            full_action_count,
            hook_product,
        })
    }

    pub(super) const fn rank(self) -> usize {
        self.rank
    }

    pub(super) const fn row_action_count(self) -> usize {
        self.row_action_count
    }

    pub(super) const fn column_action_count(self) -> usize {
        self.column_action_count
    }

    pub(super) const fn full_action_count(self) -> usize {
        self.full_action_count
    }

    pub(super) const fn hook_product(self) -> usize {
        self.hook_product
    }

    fn subgroup_order(
        group_sizes: impl IntoIterator<Item = usize>,
        quantity: YoungProjectorQuantity,
    ) -> Result<usize, YoungProjectorPlanError> {
        group_sizes.into_iter().try_fold(1usize, |order, size| {
            let factorial = (2..=size).try_fold(1usize, |partial, factor| {
                partial
                    .checked_mul(factor)
                    .ok_or(YoungProjectorPlanError::ArithmeticOverflow {
                        quantity,
                        partial,
                        factor,
                    })
            })?;
            order
                .checked_mul(factorial)
                .ok_or(YoungProjectorPlanError::ArithmeticOverflow {
                    quantity,
                    partial: order,
                    factor: factorial,
                })
        })
    }

    fn checked_hook_product(shape: &[usize]) -> Result<usize, YoungProjectorPlanError> {
        let mut product = 1usize;
        for (row, &row_length) in shape.iter().enumerate() {
            for column in 0..row_length {
                let below = shape[row + 1..]
                    .iter()
                    .filter(|&&length| length > column)
                    .count();
                let arm = row_length - column;
                let hook_length =
                    arm.checked_add(below)
                        .ok_or(YoungProjectorPlanError::ArithmeticOverflow {
                            quantity: YoungProjectorQuantity::HookProduct,
                            partial: arm,
                            factor: below,
                        })?;
                product = product.checked_mul(hook_length).ok_or(
                    YoungProjectorPlanError::ArithmeticOverflow {
                        quantity: YoungProjectorQuantity::HookProduct,
                        partial: product,
                        factor: hook_length,
                    },
                )?;
            }
        }
        Ok(product)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct YoungProjectorTerm {
    weight: Rational,
    selector: Vec<usize>,
}

impl YoungProjectorTerm {
    pub(super) fn weight(&self) -> &Rational {
        &self.weight
    }

    pub(super) fn selector(&self) -> &[usize] {
        &self.selector
    }
}

/// The normalized `C_T R_T / h_T` action after exact structural-column reduction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct YoungProjectorPlan {
    preflight: YoungProjectorPreflight,
    terms: Vec<YoungProjectorTerm>,
}

impl YoungProjectorPlan {
    #[cfg(test)]
    pub(super) fn new(tableau: &YoungTableau) -> Result<Self, YoungProjectorPlanError> {
        let preflight = YoungProjectorPreflight::new(tableau)?;
        Self::from_preflight(tableau, preflight)
    }

    fn from_preflight(
        tableau: &YoungTableau,
        preflight: YoungProjectorPreflight,
    ) -> Result<Self, YoungProjectorPlanError> {
        let rows = tableau.rows().map(<[usize]>::to_vec).collect::<Vec<_>>();
        let columns = tableau.columns().collect::<Vec<_>>();
        let row_group = Self::subgroup(
            preflight.rank(),
            &rows,
            preflight.row_action_count(),
            YoungProjectorQuantity::RowActionCount,
        )?;
        let column_group = Self::subgroup(
            preflight.rank(),
            &columns,
            preflight.column_action_count(),
            YoungProjectorQuantity::ColumnActionCount,
        )?;
        debug_assert_eq!(
            row_group.len() * column_group.len(),
            preflight.full_action_count()
        );
        let mut structural_columns = columns.clone();
        structural_columns
            .iter_mut()
            .for_each(|column| column.sort_unstable());
        let mut coefficients = BTreeMap::<Vec<usize>, Integer>::new();

        for row in &row_group {
            for column in &column_group {
                let action = column.compose(row);
                let (selector, structural_sign) =
                    Self::canonical_right_column_coset(&action, &structural_columns);
                let contribution = Integer::from((column.sign() * structural_sign) as isize);
                *coefficients.entry(selector).or_insert_with(|| 0.into()) += contribution;
            }
        }

        let denominator = Integer::from(preflight.hook_product());
        let terms = coefficients
            .into_iter()
            .filter(|(_, coefficient)| coefficient != &0)
            .map(|(selector, coefficient)| YoungProjectorTerm {
                weight: Rational::from((coefficient, denominator.clone())),
                selector,
            })
            .collect();

        Ok(Self { preflight, terms })
    }

    pub(super) const fn hook_product(&self) -> usize {
        self.preflight.hook_product()
    }

    pub(super) fn terms(&self) -> &[YoungProjectorTerm] {
        &self.terms
    }

    fn subgroup(
        rank: usize,
        groups: &[Vec<usize>],
        count: usize,
        quantity: YoungProjectorQuantity,
    ) -> Result<Vec<Permutation>, YoungProjectorPlanError> {
        let groups = groups
            .iter()
            .filter(|group| group.len() > 1)
            .map(Vec::as_slice)
            .collect::<Vec<_>>();
        let mut permutations = Vec::new();
        permutations
            .try_reserve_exact(count)
            .map_err(|_| YoungProjectorPlanError::Allocation {
                quantity,
                requested: count,
            })?;
        Self::extend_subgroup(
            &groups,
            0,
            &mut (0..rank).collect::<Vec<_>>(),
            &mut permutations,
        );
        Ok(permutations)
    }

    fn extend_subgroup(
        groups: &[&[usize]],
        group_index: usize,
        map: &mut [usize],
        permutations: &mut Vec<Permutation>,
    ) {
        let Some(group) = groups.get(group_index) else {
            permutations.push(Permutation::from_map(map.to_vec()));
            return;
        };

        for sources in group.iter().copied().permutations(group.len()) {
            for (&target, source) in group.iter().zip(sources) {
                map[target] = source;
            }
            Self::extend_subgroup(groups, group_index + 1, map, permutations);
        }
    }

    /// Lexicographically minimize a selector under the right column action.
    ///
    /// Right composition permutes selector entries only within the positions
    /// of each structural column. Sorting those entries therefore gives the
    /// same representative as enumerating the entire right column subgroup.
    fn canonical_right_column_coset(
        action: &Permutation,
        sorted_columns: &[Vec<usize>],
    ) -> (Vec<usize>, i8) {
        let selector = action.map();
        let mut representative = selector.to_vec();
        let mut structural_map = (0..selector.len()).collect::<Vec<_>>();

        for column in sorted_columns.iter().filter(|column| column.len() > 1) {
            let mut sources = column.clone();
            sources.sort_unstable_by_key(|&source| selector[source]);
            for (&target, source) in column.iter().zip(sources) {
                representative[target] = selector[source];
                structural_map[target] = source;
            }
        }

        (representative, Permutation::from_map(structural_map).sign())
    }
}

/// Signed normal form for the symmetries already carried by a projected head.
///
/// `C` permutes members inside each antisymmetric column and contributes its
/// parity. `W` exchanges whole columns of equal height without a sign: such an
/// exchange lies in the tableau row group, and therefore fixes `C_T R_T` on
/// the right. The resulting carrier stabilizer is `C ⋊ W` with character
/// `chi(c w) = sign(c)`, not the sign of the complete slot permutation.
struct YoungCarrierSymmetry {
    columns: Vec<Vec<usize>>,
    equal_height_blocks: Vec<Vec<usize>>,
}

impl YoungCarrierSymmetry {
    fn new(tableau: &YoungTableau) -> Self {
        let columns = tableau.columns().collect::<Vec<_>>();
        let mut blocks = BTreeMap::<usize, Vec<usize>>::new();
        for (column, positions) in columns.iter().enumerate() {
            blocks.entry(positions.len()).or_default().push(column);
        }
        Self {
            columns,
            equal_height_blocks: blocks
                .into_values()
                .filter(|block| block.len() > 1)
                .collect(),
        }
    }

    /// Return the source position assigned to every manifest target position.
    ///
    /// Keys are occurrence-local equality classes rather than concrete dummy
    /// names. Equal members inside one nontrivial column make the carrier zero.
    fn canonicalize<K: Ord>(&self, keys: &[K]) -> Option<(Vec<usize>, i8)> {
        let mut source_for_target = (0..keys.len()).collect::<Vec<_>>();
        let mut sign = 1i8;

        for column in self.columns.iter().filter(|column| column.len() > 1) {
            for (left, &source) in column.iter().enumerate() {
                for &right in &column[left + 1..] {
                    if keys[source] > keys[right] {
                        sign = -sign;
                    }
                }
            }
            let mut sources = column.clone();
            sources.sort_by(|left, right| {
                keys[*left].cmp(&keys[*right]).then_with(|| left.cmp(right))
            });
            if sources
                .windows(2)
                .any(|pair| keys[pair[0]] == keys[pair[1]])
            {
                return None;
            }
            for (&target, source) in column.iter().zip(sources) {
                source_for_target[target] = source;
            }
        }

        for block in &self.equal_height_blocks {
            let mut source_columns = block.clone();
            source_columns.sort_by(|left, right| {
                self.columns[*left]
                    .iter()
                    .map(|&position| &keys[source_for_target[position]])
                    .cmp(
                        self.columns[*right]
                            .iter()
                            .map(|&position| &keys[source_for_target[position]]),
                    )
                    .then_with(|| left.cmp(right))
            });
            let assignments = block
                .iter()
                .zip(source_columns)
                .flat_map(|(&target_column, source_column)| {
                    self.columns[target_column]
                        .iter()
                        .copied()
                        .zip(
                            self.columns[source_column]
                                .iter()
                                .map(|&position| source_for_target[position]),
                        )
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            for (target, source) in assignments {
                source_for_target[target] = source;
            }
        }

        Some((source_for_target, sign))
    }
}

#[derive(Clone)]
struct WeightedNetwork<Aind: AbsInd> {
    coefficient: Rational,
    network: SymbolicNet<Aind>,
    /// Declared-head form used only to preserve pre-carrier local collection.
    equivalence: Option<SymbolicNet<Aind>>,
}

impl<Aind: AbsInd> WeightedNetwork<Aind> {
    fn into_network(self) -> SymbolicNet<Aind> {
        if self.coefficient == 1 {
            self.network
        } else {
            Network::from_scalar(Atom::num(self.coefficient)) * self.network
        }
    }
}

#[derive(Clone)]
struct YoungCanonicalValue<Aind: AbsInd> {
    terms: Vec<WeightedNetwork<Aind>>,
    /// C-only projector multiplicity retained for expansion-budget accounting.
    logical_terms: usize,
    distributable: bool,
    contains_young: bool,
}

#[derive(Clone, Copy, Default)]
struct YoungScope {
    contains_young: bool,
    contains_young_power: bool,
    contains_sum: bool,
    repeated_young_head: bool,
    /// No exposed LocalTensor head or Function anywhere in the root may rerun
    /// its normalizer after graph slot assignment.
    contains_custom_normalizer: bool,
    distinct_young_tensor_root: bool,
}

#[derive(Clone, Copy)]
enum YoungTensorForm {
    Declared,
    Carrier,
}

struct YoungFactoredValue<Aind: AbsInd> {
    network: SymbolicNet<Aind>,
    linear_terms: Option<Vec<WeightedNetwork<Aind>>>,
    numeric_scalar: Option<Rational>,
    distributable: bool,
    contains_young: bool,
    logical_terms: usize,
}

#[derive(Clone, Copy)]
enum CanonicalStage {
    Intermediate,
    Complete,
}

impl<Aind: AbsInd> YoungCanonicalValue<Aind> {
    fn ordinary(network: SymbolicNet<Aind>) -> Self {
        Self {
            terms: vec![WeightedNetwork {
                coefficient: 1.into(),
                network,
                equivalence: None,
            }],
            logical_terms: 1,
            distributable: false,
            contains_young: false,
        }
    }

    fn into_network(self) -> SymbolicNet<Aind> {
        self.terms
            .into_iter()
            .map(WeightedNetwork::into_network)
            .reduce(|sum, term| sum + term)
            .unwrap_or_else(Network::zero)
    }

    fn into_factor_terms(self) -> Vec<WeightedNetwork<Aind>> {
        if self.distributable {
            self.terms
        } else {
            vec![WeightedNetwork {
                coefficient: 1.into(),
                network: self.into_network(),
                equivalence: None,
            }]
        }
    }

    fn boundary(network: SymbolicNet<Aind>, contains_young: bool) -> Self {
        let mut value = Self::ordinary(network);
        value.contains_young = contains_young;
        value
    }
}

impl<Aind: AbsInd> YoungFactoredValue<Aind> {
    fn ordinary(network: SymbolicNet<Aind>) -> Self {
        Self {
            network,
            linear_terms: None,
            numeric_scalar: None,
            distributable: false,
            contains_young: false,
            logical_terms: 1,
        }
    }

    fn boundary(network: SymbolicNet<Aind>, contains_young: bool) -> Self {
        Self {
            network,
            linear_terms: None,
            numeric_scalar: None,
            distributable: false,
            contains_young,
            logical_terms: 1,
        }
    }
}

struct YoungNetworkFold<'a, Aind: AbsInd> {
    network: &'a SymbolicNet<Aind>,
    tree: SimpleTraversalTree<ParentChildStore<()>>,
    source_indices: OnceCell<Vec<Atom>>,
    young_tableaux: RefCell<BTreeMap<Symbol, Option<YoungTableau>>>,
}

impl<'a, Aind> YoungNetworkFold<'a, Aind>
where
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
    SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
{
    fn new(network: &'a SymbolicNet<Aind>) -> Self {
        Self {
            network,
            tree: network.graph.expr_tree().cast(),
            source_indices: OnceCell::new(),
            young_tableaux: RefCell::new(BTreeMap::new()),
        }
    }

    fn source_indices(&self) -> &[Atom] {
        self.source_indices.get_or_init(|| {
            self.network
                .store
                .tensors
                .iter()
                .flat_map(|tensor| {
                    (&tensor.structure)
                        .into_iter()
                        .map(|slot| slot.aind.to_atom())
                })
                .collect()
        })
    }

    #[cfg(test)]
    fn record_factored_graph_outcome(&self, attempt: bool) {
        let heads = self
            .network
            .store
            .tensors
            .iter()
            .filter_map(|tensor| match tensor.expression.as_view() {
                AtomView::Fun(function) => Some(function.get_symbol()),
                _ => None,
            })
            .collect::<BTreeSet<_>>();
        let mut outcomes = FACTORED_GRAPH_OUTCOMES
            .lock()
            .expect("factored graph outcome lock is not poisoned");
        for head in heads {
            let (attempts, fallbacks) = outcomes.entry(head).or_default();
            if attempt {
                *attempts += 1;
            } else {
                *fallbacks += 1;
            }
        }
    }

    #[cfg(test)]
    fn record_young_route(&self, route: YoungRoute) {
        let heads = self
            .network
            .store
            .tensors
            .iter()
            .filter_map(|tensor| self.general_young_tableau(tensor).ok().flatten())
            .map(|(head, _)| head)
            .collect::<BTreeSet<_>>();
        let mut outcomes = YOUNG_ROUTE_OUTCOMES
            .lock()
            .expect("Young route outcome lock is not poisoned");
        for head in heads {
            let outcome = outcomes.entry(head).or_default();
            match route {
                YoungRoute::Carrier => outcome.carrier += 1,
                YoungRoute::Composite => outcome.composite += 1,
                YoungRoute::Staged => outcome.staged += 1,
            }
        }
    }

    fn run<F>(
        &self,
        request: &mut CanonicalRequest<'_, Aind, F>,
    ) -> Result<YoungCanonicalValue<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        let root = self.network.graph.graph.node_id(self.network.graph.head());
        self.node(root, request)
    }

    fn scope(&self) -> Result<YoungScope, CanonicalizationError> {
        let root = self.network.graph.graph.node_id(self.network.graph.head());
        let mut scope = self.node_scope(root, &mut BTreeSet::new())?;
        if let NetworkNode::Leaf(NetworkLeaf::LocalTensor(tensor)) = self.network.graph.graph[root]
        {
            let tensor = &self.network.store.tensors[tensor];
            scope.distinct_young_tensor_root = self.general_young_tableau(tensor)?.is_some()
                && !self
                    .network
                    .graph
                    .graph
                    .iter_crown(root)
                    .any(|hedge| self.network.graph.graph.is_self_loop(hedge));
        }
        Ok(scope)
    }

    fn node_scope(
        &self,
        node: NodeIndex,
        young_heads: &mut BTreeSet<Symbol>,
    ) -> Result<YoungScope, CanonicalizationError> {
        let mut scope = YoungScope::default();
        for child in self
            .tree
            .iter_children(node, self.network.graph.graph.as_ref())
        {
            let child = self.node_scope(child, young_heads)?;
            scope.contains_young |= child.contains_young;
            scope.contains_young_power |= child.contains_young_power;
            scope.contains_sum |= child.contains_sum;
            scope.repeated_young_head |= child.repeated_young_head;
            scope.contains_custom_normalizer |= child.contains_custom_normalizer;
        }

        match self.network.graph.graph[node].clone() {
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(tensor)) => {
                let tensor = &self.network.store.tensors[tensor];
                if let AtomView::Fun(function) = tensor.expression.as_view() {
                    scope.contains_custom_normalizer =
                        function.get_symbol().get_normalization_function().is_some();
                }
                if let Some((head, _)) = self.general_young_tableau(tensor)? {
                    scope.contains_young = true;
                    scope.repeated_young_head = !young_heads.insert(head);
                }
            }
            NetworkNode::Op(NetworkOp::Function(function)) => {
                scope.contains_custom_normalizer |= function.get_normalization_function().is_some();
            }
            NetworkNode::Op(NetworkOp::Power(_)) => {
                scope.contains_young_power |= scope.contains_young;
            }
            NetworkNode::Op(NetworkOp::Sum) => scope.contains_sum = true,
            NetworkNode::Leaf(_) | NetworkNode::Op(_) => {}
        }
        Ok(scope)
    }

    fn run_factored(
        &self,
        form: YoungTensorForm,
        collect_linear_terms: bool,
    ) -> Result<YoungFactoredValue<Aind>, CanonicalizationError> {
        let root = self.network.graph.graph.node_id(self.network.graph.head());
        self.factored_node(root, form, collect_linear_terms)
    }

    fn factored_policy(
        &self,
        form: YoungTensorForm,
        collect_linear_terms: bool,
    ) -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError> {
        let value = self.run_factored(form, collect_linear_terms)?;
        CanonicalPolicyNet::parse(execute_atom(value.network)?)
    }

    fn factored_node(
        &self,
        node: NodeIndex,
        form: YoungTensorForm,
        collect_linear_terms: bool,
    ) -> Result<YoungFactoredValue<Aind>, CanonicalizationError> {
        let children = self
            .tree
            .iter_children(node, self.network.graph.graph.as_ref())
            .collect::<Vec<_>>();
        match self.network.graph.graph[node].clone() {
            NetworkNode::Leaf(NetworkLeaf::Scalar(reference)) => {
                let scalar = self.network.store.get_scalar_ref(reference).clone();
                let numeric_scalar = Rational::try_from(scalar.as_view()).ok();
                let mut value = YoungFactoredValue::ordinary(Network::from_scalar(scalar));
                value.numeric_scalar = numeric_scalar;
                Ok(value)
            }
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(tensor)) => {
                let tensor = &self.network.store.tensors[tensor];
                let Some(projected) = self.tensor(tensor, form, collect_linear_terms)? else {
                    return Ok(YoungFactoredValue::ordinary(Network::from_tensor(
                        tensor.clone(),
                    )));
                };
                let logical_terms = projected.logical_terms;
                let (network, linear_terms) = if collect_linear_terms {
                    (projected.clone().into_network(), Some(projected.terms))
                } else {
                    (projected.into_network(), None)
                };
                Ok(YoungFactoredValue {
                    network,
                    linear_terms,
                    numeric_scalar: None,
                    distributable: true,
                    contains_young: true,
                    logical_terms,
                })
            }
            NetworkNode::Leaf(_) => Err(CanonicalizationError::UnsupportedLeaf),
            NetworkNode::Op(NetworkOp::Product) => {
                let mut accumulator = None;
                for child in children {
                    let value = self.factored_node(child, form, collect_linear_terms)?;
                    accumulator = Some(match accumulator {
                        Some(accumulator) => Self::factored_product(accumulator, value)?,
                        None => value,
                    });
                }
                Ok(accumulator.unwrap_or_else(|| YoungFactoredValue::ordinary(Network::one())))
            }
            NetworkNode::Op(NetworkOp::Sum) => {
                let mut ordinary = Vec::with_capacity(children.len());
                let mut combined = Vec::<WeightedNetwork<Aind>>::new();
                let mut contains_young = false;
                let mut requested_terms = 0usize;
                for child in children {
                    let mut value = self.factored_node(child, form, collect_linear_terms)?;
                    requested_terms = requested_terms.saturating_add(if value.distributable {
                        value.logical_terms
                    } else {
                        1
                    });
                    let previously_contained_young = contains_young;
                    contains_young |= value.contains_young;
                    if contains_young {
                        Self::check_live_terms(requested_terms)?;
                    }
                    if !contains_young {
                        ordinary.push(value);
                        continue;
                    }

                    if !previously_contained_young {
                        for value in ordinary.drain(..) {
                            Self::insert_local_term(
                                &mut combined,
                                WeightedNetwork {
                                    coefficient: 1.into(),
                                    network: value.network,
                                    equivalence: None,
                                },
                            );
                            Self::check_live_terms(combined.len())?;
                        }
                    }
                    let candidates = value.linear_terms.take().unwrap_or_else(|| {
                        vec![WeightedNetwork {
                            coefficient: 1.into(),
                            network: value.network,
                            equivalence: None,
                        }]
                    });
                    for candidate in candidates {
                        Self::insert_local_term(&mut combined, candidate);
                        Self::check_live_terms(combined.len())?;
                    }
                }
                if !contains_young {
                    let network = ordinary
                        .into_iter()
                        .map(|value| value.network)
                        .reduce(|sum, term| sum + term)
                        .unwrap_or_else(Network::zero);
                    return Ok(YoungFactoredValue::boundary(network, false));
                }
                let mut candidates = combined;
                let content = candidates
                    .iter()
                    .map(|term| term.coefficient.abs())
                    .reduce(|left, right| left.gcd(&right))
                    .unwrap_or_else(Rational::one);
                if content != 1 {
                    let inverse = content.inv();
                    candidates.iter_mut().for_each(|term| {
                        term.coefficient = &term.coefficient * &inverse;
                    });
                }
                let primitive = candidates
                    .into_iter()
                    .map(WeightedNetwork::into_network)
                    .reduce(|sum, term| sum + term)
                    .unwrap_or_else(Network::zero);
                let network = if content == 1 {
                    primitive
                } else {
                    Network::from_scalar(Atom::num(content)) * primitive
                };
                Ok(YoungFactoredValue::boundary(network, true))
            }
            NetworkNode::Op(NetworkOp::Neg) => {
                let mut value = self.only_factored_child(children, form, collect_linear_terms)?;
                value.network = -value.network;
                if let Some(terms) = &mut value.linear_terms {
                    terms.iter_mut().for_each(|term| {
                        term.coefficient = -term.coefficient.clone();
                    });
                }
                value.numeric_scalar = value.numeric_scalar.map(|coefficient| -coefficient);
                Ok(value)
            }
            NetworkNode::Op(NetworkOp::Function(function)) => {
                let value = self.only_factored_child(children, form, collect_linear_terms)?;
                Ok(YoungFactoredValue::boundary(
                    value.network.fun(function),
                    value.contains_young,
                ))
            }
            NetworkNode::Op(NetworkOp::Power(0)) => {
                Ok(YoungFactoredValue::ordinary(Network::one()))
            }
            NetworkNode::Op(NetworkOp::Power(exponent)) => {
                let value = self.only_factored_child(children, form, collect_linear_terms)?;
                debug_assert!(!value.contains_young);
                Ok(YoungFactoredValue::ordinary(value.network.pow(exponent)))
            }
        }
    }

    fn only_factored_child(
        &self,
        children: Vec<NodeIndex>,
        form: YoungTensorForm,
        collect_linear_terms: bool,
    ) -> Result<YoungFactoredValue<Aind>, CanonicalizationError> {
        let [child] = children.as_slice() else {
            return Err(CanonicalizationError::Projection(
                "unary factored Young operation does not have exactly one child".to_owned(),
            ));
        };
        self.factored_node(*child, form, collect_linear_terms)
    }

    fn factored_product(
        mut accumulator: YoungFactoredValue<Aind>,
        mut value: YoungFactoredValue<Aind>,
    ) -> Result<YoungFactoredValue<Aind>, CanonicalizationError> {
        let contains_young = accumulator.contains_young || value.contains_young;
        if accumulator.distributable || value.distributable {
            let left_terms = if accumulator.distributable {
                accumulator.logical_terms
            } else {
                1
            };
            let right_terms = if value.distributable {
                value.logical_terms
            } else {
                1
            };
            let requested_terms = left_terms.saturating_mul(right_terms);
            Self::check_live_terms(requested_terms)?;
            let linear_terms = match (accumulator.linear_terms.take(), value.linear_terms.take()) {
                (Some(mut terms), None) if !value.contains_young => {
                    if let Some(coefficient) = &value.numeric_scalar {
                        terms.iter_mut().for_each(|term| {
                            term.coefficient = &term.coefficient * coefficient;
                        });
                    } else {
                        terms.iter_mut().for_each(|term| {
                            if let Some(equivalence) = &mut term.equivalence {
                                *equivalence = equivalence.clone() * value.network.clone();
                            }
                            term.network = term.network.clone() * value.network.clone();
                        });
                    }
                    Some(terms)
                }
                (None, Some(mut terms)) if !accumulator.contains_young => {
                    if let Some(coefficient) = &accumulator.numeric_scalar {
                        terms.iter_mut().for_each(|term| {
                            term.coefficient = coefficient * &term.coefficient;
                        });
                    } else {
                        terms.iter_mut().for_each(|term| {
                            if let Some(equivalence) = &mut term.equivalence {
                                *equivalence = accumulator.network.clone() * equivalence.clone();
                            }
                            term.network = accumulator.network.clone() * term.network.clone();
                        });
                    }
                    Some(terms)
                }
                _ => None,
            };
            Ok(YoungFactoredValue {
                network: accumulator.network * value.network,
                linear_terms,
                numeric_scalar: None,
                distributable: true,
                contains_young,
                logical_terms: requested_terms,
            })
        } else {
            let numeric_scalar = if contains_young {
                None
            } else {
                accumulator
                    .numeric_scalar
                    .as_ref()
                    .zip(value.numeric_scalar.as_ref())
                    .map(|(left, right)| left * right)
            };
            Ok(YoungFactoredValue {
                network: accumulator.network * value.network,
                linear_terms: None,
                numeric_scalar,
                distributable: false,
                contains_young,
                logical_terms: 1,
            })
        }
    }

    fn insert_local_term(
        combined: &mut Vec<WeightedNetwork<Aind>>,
        candidate: WeightedNetwork<Aind>,
    ) {
        if candidate.coefficient == 0 {
            return;
        }
        if let Some(position) = combined.iter().position(|existing| {
            existing.equivalence.as_ref().unwrap_or(&existing.network)
                == candidate.equivalence.as_ref().unwrap_or(&candidate.network)
        }) {
            combined[position].coefficient += &candidate.coefficient;
            if combined[position].coefficient == 0 {
                combined.swap_remove(position);
            }
        } else {
            combined.push(candidate);
        }
    }

    fn check_live_terms(requested_terms: usize) -> Result<(), CanonicalizationError> {
        if requested_terms > MAX_YOUNG_EXPANSION_TERMS {
            Err(CanonicalizationError::YoungExpansionSizeLimit {
                requested_terms,
                term_limit: MAX_YOUNG_EXPANSION_TERMS,
            })
        } else {
            Ok(())
        }
    }

    fn node<F>(
        &self,
        node: NodeIndex,
        request: &mut CanonicalRequest<'_, Aind, F>,
    ) -> Result<YoungCanonicalValue<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        let children = self
            .tree
            .iter_children(node, self.network.graph.graph.as_ref())
            .collect::<Vec<_>>();
        match self.network.graph.graph[node].clone() {
            NetworkNode::Leaf(NetworkLeaf::Scalar(reference)) => Ok(YoungCanonicalValue::ordinary(
                Network::from_scalar(self.network.store.get_scalar_ref(reference).clone()),
            )),
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(tensor)) => {
                let tensor = &self.network.store.tensors[tensor];
                self.tensor(tensor, YoungTensorForm::Declared, false)
                    .map(|projected| {
                        projected.unwrap_or_else(|| {
                            YoungCanonicalValue::ordinary(Network::from_tensor(tensor.clone()))
                        })
                    })
            }
            NetworkNode::Leaf(_) => Err(CanonicalizationError::UnsupportedLeaf),
            NetworkNode::Op(NetworkOp::Product) => {
                let mut children = children.into_iter().peekable();
                let Some(first) = children.next() else {
                    return Ok(YoungCanonicalValue::ordinary(Network::one()));
                };
                let mut accumulator = self.node(first, request)?;
                while let Some(child) = children.next() {
                    let value = self.node(child, request)?;
                    accumulator =
                        self.product_pair(accumulator, value, children.peek().is_some(), request)?;
                }
                Ok(accumulator)
            }
            NetworkNode::Op(NetworkOp::Sum) => {
                let mut values = Vec::with_capacity(children.len());
                let mut contains_young = false;
                let mut requested_terms = 0usize;
                for child in children {
                    let value = self.node(child, request)?;
                    requested_terms = requested_terms.saturating_add(if value.distributable {
                        value.logical_terms
                    } else {
                        1
                    });
                    contains_young |= value.contains_young;
                    if contains_young {
                        Self::check_live_terms(requested_terms)?;
                    }
                    values.push(value);
                }
                if !contains_young {
                    let network = values
                        .into_iter()
                        .map(YoungCanonicalValue::into_network)
                        .reduce(|sum, term| sum + term)
                        .unwrap_or_else(Network::zero);
                    return Ok(YoungCanonicalValue::boundary(network, false));
                }
                let candidates = values
                    .into_iter()
                    .flat_map(|value| {
                        if value.distributable {
                            value.terms
                        } else {
                            vec![WeightedNetwork {
                                coefficient: 1.into(),
                                network: value.into_network(),
                                equivalence: None,
                            }]
                        }
                    })
                    .collect();
                Ok(YoungCanonicalValue {
                    terms: self.canonicalize_terms(
                        candidates,
                        request,
                        CanonicalStage::Intermediate,
                    )?,
                    logical_terms: 1,
                    distributable: false,
                    contains_young: true,
                })
            }
            NetworkNode::Op(NetworkOp::Neg) => {
                let mut value = self.only_child(children, request)?;
                value
                    .terms
                    .iter_mut()
                    .for_each(|term| term.coefficient = -term.coefficient.clone());
                Ok(value)
            }
            NetworkNode::Op(NetworkOp::Function(function)) => {
                let value = self.only_child(children, request)?;
                let contains_young = value.contains_young;
                let argument = if value.terms.is_empty() {
                    Network::from_scalar(Atom::Zero)
                } else {
                    value.into_network()
                };
                Ok(YoungCanonicalValue::boundary(
                    argument.fun(function),
                    contains_young,
                ))
            }
            NetworkNode::Op(NetworkOp::Power(0)) => {
                Ok(YoungCanonicalValue::ordinary(Network::one()))
            }
            NetworkNode::Op(NetworkOp::Power(exponent)) => {
                let mut value = self.only_child(children, request)?;
                let contains_young = value.contains_young;
                if !contains_young {
                    return Ok(YoungCanonicalValue::boundary(
                        value.into_network().pow(exponent),
                        false,
                    ));
                }

                let magnitude = usize::from(exponent.unsigned_abs());
                let requested_terms = if value.distributable {
                    (0..magnitude)
                        .fold(1usize, |terms, _| terms.saturating_mul(value.logical_terms))
                } else {
                    value.logical_terms.saturating_mul(magnitude)
                };
                Self::check_live_terms(requested_terms)?;
                value.terms =
                    self.canonicalize_terms(value.terms, request, CanonicalStage::Intermediate)?;
                if value.terms.is_empty() {
                    return Ok(YoungCanonicalValue::boundary(
                        Network::from_scalar(Atom::Zero).pow(exponent),
                        true,
                    ));
                }

                if !value.distributable && magnitude.is_multiple_of(2) {
                    let atom = execute_atom(value.into_network())?;
                    let policy = CanonicalPolicyNet::parse(atom)?;
                    let source_indices = self.source_indices();
                    let namespace = request.reserve_temporary_namespace(source_indices.len())?;
                    let policy = request.canonize_temporary_graph(
                        policy,
                        source_indices,
                        namespace,
                        ExternalLineMode::AnonymousEvenPower,
                    )?;
                    value = YoungCanonicalValue::boundary(policy.network().clone(), true);
                }
                let mut copies = Vec::with_capacity(magnitude);
                copies.push(value.clone());
                for _ in 1..magnitude {
                    copies.push(self.fresh_power_copy(&value, request)?);
                }
                let product = self.product(copies, request)?;
                if exponent < 0 {
                    if product.terms.is_empty() {
                        Ok(YoungCanonicalValue::boundary(
                            Network::from_scalar(Atom::Zero).pow(exponent),
                            true,
                        ))
                    } else {
                        Ok(YoungCanonicalValue::boundary(
                            product.into_network().pow(-1),
                            true,
                        ))
                    }
                } else {
                    Ok(product)
                }
            }
        }
    }

    /// Materialize one Power copy through the graph-aware temporary rebuilder.
    ///
    /// That existing path retains external slots, assigns internal lines by
    /// representation and ordinal, and rebuilds tensor payloads structurally.
    /// A fresh namespace per copy prevents later execution from merging copy
    /// scopes when it renders the product back into one Atom.
    fn fresh_power_copy<F>(
        &self,
        value: &YoungCanonicalValue<Aind>,
        request: &mut CanonicalRequest<'_, Aind, F>,
    ) -> Result<YoungCanonicalValue<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        Ok(YoungCanonicalValue {
            terms: self.canonicalize_terms(
                value.terms.clone(),
                request,
                CanonicalStage::Intermediate,
            )?,
            logical_terms: value.logical_terms,
            distributable: value.distributable,
            contains_young: value.contains_young,
        })
    }

    fn only_child<F>(
        &self,
        children: Vec<NodeIndex>,
        request: &mut CanonicalRequest<'_, Aind, F>,
    ) -> Result<YoungCanonicalValue<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        let [child] = children.as_slice() else {
            return Err(CanonicalizationError::Projection(
                "unary Young fold operation does not have exactly one child".to_owned(),
            ));
        };
        self.node(*child, request)
    }

    fn general_young_tableau(
        &self,
        tensor: &SymbolicTensor<Aind>,
    ) -> Result<Option<(Symbol, YoungTableau)>, CanonicalizationError> {
        let AtomView::Fun(function) = tensor.expression.as_view() else {
            return Ok(None);
        };
        let head = function.get_symbol();
        let tableau = if let Some(cached) = self.young_tableaux.borrow().get(&head).cloned() {
            cached
        } else {
            let tableau = YoungTableau::from_symbol(head).map_err(|error| {
                CanonicalizationError::InvalidYoungTableauMetadata {
                    head,
                    reason: error.to_string(),
                }
            })?;
            self.young_tableaux
                .borrow_mut()
                .insert(head, tableau.clone());
            tableau
        };
        let Some(tableau) = tableau else {
            return Ok(None);
        };
        Ok((tableau.class() == YoungTableauClass::General).then_some((head, tableau)))
    }

    fn tensor(
        &self,
        tensor: &SymbolicTensor<Aind>,
        form: YoungTensorForm,
        preserve_equivalence: bool,
    ) -> Result<Option<YoungCanonicalValue<Aind>>, CanonicalizationError> {
        let AtomView::Fun(function) = tensor.expression.as_view() else {
            return Ok(None);
        };
        let Some((head, tableau)) = self.general_young_tableau(tensor)? else {
            return Ok(None);
        };

        let preflight = YoungProjectorPreflight::new(&tableau).map_err(|error| {
            CanonicalizationError::YoungProjectorPlanning {
                head,
                reason: error.to_string(),
            }
        })?;
        if preflight.full_action_count() > MAX_YOUNG_PROJECTOR_ACTIONS {
            return Err(CanonicalizationError::YoungProjectorSizeLimit {
                head,
                shape: tableau.shape().to_vec(),
                requested_actions: preflight.full_action_count(),
                action_limit: MAX_YOUNG_PROJECTOR_ACTIONS,
            });
        }
        let plan = YoungProjectorPlan::from_preflight(&tableau, preflight).map_err(|error| {
            CanonicalizationError::YoungProjectorPlanning {
                head,
                reason: error.to_string(),
            }
        })?;
        debug_assert_eq!(plan.hook_product(), preflight.hook_product());
        let mut slots = Vec::<LibrarySlot<Aind>>::with_capacity(preflight.rank());
        let mut arguments = Vec::with_capacity(preflight.rank());
        for argument in function.iter() {
            slots.push(Slot::<LibraryRep, Aind>::try_from(argument).map_err(|_| {
                CanonicalizationError::StructureMismatch {
                    expression: tensor.expression.clone(),
                }
            })?);
            arguments.push(argument.to_owned());
        }
        let mut first_occurrences = BTreeMap::<SemanticAtomKey, usize>::new();
        let argument_classes: Vec<usize> = arguments
            .iter()
            .enumerate()
            .map(|(position, argument)| {
                *first_occurrences
                    .entry(SemanticAtomKey::new(argument.as_view()))
                    .or_insert(position)
            })
            .collect();
        let carrier_symmetry = YoungCarrierSymmetry::new(&tableau);
        let carrier_symbol =
            matches!(form, YoungTensorForm::Carrier).then(YoungColumnCarrier::symbol);
        let logical_terms = plan.terms().len();
        let mut combined = BTreeMap::<Vec<usize>, (Vec<usize>, Rational)>::new();
        for term in plan.terms() {
            let selected_classes = term
                .selector()
                .iter()
                .map(|&position| argument_classes[position])
                .collect::<Vec<_>>();
            let Some((source_for_target, column_sign)) =
                carrier_symmetry.canonicalize(&selected_classes)
            else {
                continue;
            };
            let selected_classes = source_for_target
                .iter()
                .map(|&source| selected_classes[source])
                .collect();
            let negative = column_sign < 0;
            let coefficient = if negative {
                -term.weight().clone()
            } else {
                term.weight().clone()
            };
            combined
                .entry(selected_classes)
                .and_modify(|(_, existing)| *existing += &coefficient)
                .or_insert_with(|| {
                    let selected_positions = source_for_target
                        .iter()
                        .map(|&source| term.selector()[source])
                        .collect();
                    (selected_positions, coefficient)
                });
        }
        let terms = combined
            .into_values()
            .filter(|(_, coefficient)| coefficient != &0)
            .map(|(selected_positions, coefficient)| {
                let selected_slots = selected_positions
                    .iter()
                    .map(|&position| slots[position])
                    .collect();
                let expression = selected_positions
                    .into_iter()
                    .map(|position| arguments[position].clone())
                    .fold(FunctionBuilder::new(head), FunctionBuilder::add_arg)
                    .finish();
                let tensor = SymbolicTensor {
                    structure: OrderedStructure::new(selected_slots).structure,
                    is_metric: tensor.is_metric,
                    is_composite: tensor.is_composite,
                    expression,
                };
                let equivalence = (preserve_equivalence && carrier_symbol.is_some())
                    .then(|| Network::from_tensor(tensor.clone()));
                let tensor = match carrier_symbol {
                    Some(carrier_symbol) => YoungColumnCarrier::tensor_with_symbol(
                        tensor,
                        &tableau,
                        &carrier_symmetry.columns,
                        carrier_symbol,
                    )?,
                    None => tensor,
                };
                Ok(WeightedNetwork {
                    coefficient,
                    network: Network::from_tensor(tensor),
                    equivalence,
                })
            })
            .collect::<Result<Vec<_>, CanonicalizationError>>()?;

        Ok(Some(YoungCanonicalValue {
            terms,
            logical_terms,
            distributable: true,
            contains_young: true,
        }))
    }

    fn product<F>(
        &self,
        values: Vec<YoungCanonicalValue<Aind>>,
        request: &mut CanonicalRequest<'_, Aind, F>,
    ) -> Result<YoungCanonicalValue<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        let mut values = values.into_iter().peekable();
        let Some(mut accumulator) = values.next() else {
            return Ok(YoungCanonicalValue::ordinary(Network::one()));
        };
        while let Some(value) = values.next() {
            accumulator =
                self.product_pair(accumulator, value, values.peek().is_some(), request)?;
        }
        Ok(accumulator)
    }

    fn product_pair<F>(
        &self,
        accumulator: YoungCanonicalValue<Aind>,
        value: YoungCanonicalValue<Aind>,
        has_more_factors: bool,
        request: &mut CanonicalRequest<'_, Aind, F>,
    ) -> Result<YoungCanonicalValue<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        let contains_young = accumulator.contains_young || value.contains_young;
        if accumulator.distributable || value.distributable {
            let logical_terms = (if accumulator.distributable {
                accumulator.logical_terms
            } else {
                1
            })
            .saturating_mul(if value.distributable {
                value.logical_terms
            } else {
                1
            });
            Self::check_live_terms(logical_terms)?;
            let left_terms = accumulator.into_factor_terms();
            let right_terms = value.into_factor_terms();
            let requested_terms = left_terms.len().saturating_mul(right_terms.len());
            let mut candidates = Vec::new();
            candidates.try_reserve_exact(requested_terms).map_err(|_| {
                CanonicalizationError::YoungExpansionSizeLimit {
                    requested_terms,
                    term_limit: MAX_YOUNG_EXPANSION_TERMS,
                }
            })?;
            for left in left_terms {
                for right in &right_terms {
                    candidates.push(WeightedNetwork {
                        coefficient: &left.coefficient * &right.coefficient,
                        network: left.network.clone() * right.network.clone(),
                        equivalence: None,
                    });
                }
            }
            let terms = if has_more_factors {
                self.canonicalize_terms(candidates, request, CanonicalStage::Intermediate)?
            } else {
                candidates
            };
            Ok(YoungCanonicalValue {
                terms,
                logical_terms,
                distributable: true,
                contains_young,
            })
        } else {
            Ok(YoungCanonicalValue::boundary(
                accumulator.into_network() * value.into_network(),
                contains_young,
            ))
        }
    }

    fn canonicalize_terms<F>(
        &self,
        candidates: Vec<WeightedNetwork<Aind>>,
        request: &mut CanonicalRequest<'_, Aind, F>,
        stage: CanonicalStage,
    ) -> Result<Vec<WeightedNetwork<Aind>>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        let temporary_namespace = match stage {
            CanonicalStage::Intermediate => {
                Some(request.reserve_temporary_namespace(self.source_indices().len())?)
            }
            CanonicalStage::Complete => None,
        };
        let mut combined = BTreeMap::<SemanticAtomKey, (Atom, Rational)>::new();
        for candidate in candidates {
            if candidate.coefficient == 0 {
                continue;
            }
            let atom = execute_atom(candidate.network)?;
            let policy = CanonicalPolicyNet::parse(atom)?;
            let canonical = match stage {
                CanonicalStage::Intermediate => request.canonize_temporary_graph(
                    policy,
                    self.source_indices(),
                    temporary_namespace.expect("intermediate stage reserved a namespace"),
                    ExternalLineMode::Preserve,
                )?,
                CanonicalStage::Complete => request.canonize_graph(policy)?,
            }
            .into_atom();
            if canonical.as_view().is_zero() {
                continue;
            }
            let positive_key = SemanticAtomKey::new(canonical.as_view());
            let negative = -canonical.clone();
            let negative_key = SemanticAtomKey::new(negative.as_view());
            let (key, atom, coefficient) = if negative_key < positive_key {
                (negative_key, negative, -candidate.coefficient)
            } else {
                (positive_key, canonical, candidate.coefficient)
            };
            combined
                .entry(key)
                .and_modify(|(_, existing)| *existing += &coefficient)
                .or_insert((atom, coefficient));
        }

        combined
            .into_values()
            .filter(|(_, coefficient)| coefficient != &0)
            .map(|(atom, coefficient)| {
                let policy = CanonicalPolicyNet::parse(atom)?;
                Ok(WeightedNetwork {
                    coefficient,
                    network: policy.network().clone(),
                    equivalence: None,
                })
            })
            .collect()
    }

    fn finish<F>(
        &self,
        mut value: YoungCanonicalValue<Aind>,
        request: &mut CanonicalRequest<'_, Aind, F>,
    ) -> Result<CanonicalPolicyNet<Aind>, CanonicalizationError>
    where
        F: FnMut(usize) -> Aind,
    {
        let aggregate = execute_atom(value.clone().into_network())?;
        let policy = CanonicalPolicyNet::parse(aggregate)?;
        match request.canonize_graph(policy) {
            Ok(canonical) => return Ok(canonical),
            Err(
                CanonicalizationError::GraphSizeLimit { .. }
                | CanonicalizationError::WholeGraphSizeLimit { .. },
            ) => {}
            Err(error) => return Err(error),
        }
        value.terms = self.canonicalize_terms(value.terms, request, CanonicalStage::Complete)?;
        value.distributable = false;
        CanonicalPolicyNet::parse(execute_atom(value.into_network())?)
    }
}

// Keep the parser-proven policy inline: boxing it would allocate on every
// ordinary/non-Young request merely to shrink this short-lived result enum.
#[allow(clippy::large_enum_variant)]
pub(super) enum YoungStraightened<Aind: AbsInd> {
    Policy {
        policy: CanonicalPolicyNet<Aind>,
        graph_canonical: bool,
    },
    Terminal(Atom),
}

pub(super) fn straighten<Aind, F>(
    policy: CanonicalPolicyNet<Aind>,
    request: &mut CanonicalRequest<'_, Aind, F>,
) -> Result<YoungStraightened<Aind>, CanonicalizationError>
where
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
    F: FnMut(usize) -> Aind,
    SymbolicTensor<Aind>: Contract<SymbolicTensor<Aind>, (), LCM = SymbolicTensor<Aind>>,
{
    let fold = YoungNetworkFold::new(policy.network());
    let scope = fold.scope()?;
    if !scope.contains_young {
        return Ok(YoungStraightened::Policy {
            policy,
            graph_canonical: false,
        });
    }
    if !scope.contains_young_power
        && !scope.repeated_young_head
        && !scope.contains_custom_normalizer
    {
        // A lone declared tensor with distinct external lines needs no private
        // carrier to align head ordering or dummy frames. Canonicalize its exact
        // projected sum directly, then select the shared byte-stable numeric
        // section before returning the terminal Atom.
        if scope.distinct_young_tensor_root {
            let declared = fold.factored_policy(YoungTensorForm::Declared, false)?;
            match request.canonize_graph(declared) {
                Ok(canonical) => {
                    let normalized =
                        YoungColumnCarrier::normalize_numeric_content(canonical.into_atom());
                    return Ok(YoungStraightened::Terminal(normalized));
                }
                Err(
                    CanonicalizationError::ConvergenceCycle { .. }
                    | CanonicalizationError::GraphSizeLimit { .. }
                    | CanonicalizationError::WholeGraphSizeLimit { .. },
                ) => {}
                Err(error) => return Err(error),
            }
        }
        let carrier = fold.factored_policy(YoungTensorForm::Carrier, scope.contains_sum)?;
        #[cfg(test)]
        fold.record_factored_graph_outcome(true);
        match request.canonize_graph(carrier) {
            Ok(canonical) => {
                let decoded = YoungColumnCarrier::decode_policy(canonical)?;
                #[cfg(test)]
                fold.record_young_route(YoungRoute::Carrier);
                return Ok(YoungStraightened::Policy {
                    policy: decoded,
                    graph_canonical: true,
                });
            }
            Err(CanonicalizationError::ConvergenceCycle { .. }) => {
                let factored =
                    fold.factored_policy(YoungTensorForm::Declared, scope.contains_sum)?;
                match request.canonize_graph_modulo(factored, |graph| {
                    let fold = YoungNetworkFold::new(graph.network());
                    let scope = fold.scope()?;
                    fold.factored_policy(YoungTensorForm::Declared, scope.contains_sum)
                }) {
                    Ok(canonical) => {
                        #[cfg(test)]
                        fold.record_young_route(YoungRoute::Composite);
                        return Ok(YoungStraightened::Policy {
                            policy: canonical,
                            graph_canonical: true,
                        });
                    }
                    Err(
                        CanonicalizationError::GraphSizeLimit { .. }
                        | CanonicalizationError::WholeGraphSizeLimit { .. },
                    ) => {}
                    Err(error) => return Err(error),
                }
            }
            Err(
                CanonicalizationError::GraphSizeLimit { .. }
                | CanonicalizationError::WholeGraphSizeLimit { .. }
                | CanonicalizationError::SignedDecorationOrbitLimit { .. },
            ) => {}
            Err(error) => return Err(error),
        }
        #[cfg(test)]
        {
            fold.record_factored_graph_outcome(false);
            fold.record_young_route(YoungRoute::Staged);
        }
    }

    let value = fold.run(request)?;
    if !value.contains_young {
        return Ok(YoungStraightened::Policy {
            policy,
            graph_canonical: false,
        });
    }
    // Internal stages use source-disjoint temporary dummies. Only complete
    // root terms consume the request's caller-provided dummy allocator.
    Ok(YoungStraightened::Policy {
        policy: fold.finish(value, request)?,
        graph_canonical: true,
    })
}

#[cfg(test)]
mod tests {
    use spenso::{
        bracket, slot,
        structure::{abstract_index::AbstractIndex, slot::IsAbstractSlot},
        tensor_symbol,
    };
    use symbolica::{atom::AtomCore, function, symbol};

    use crate::{
        IndexTooling, reference_cases::young::YoungProjector, test_support::test_initialize,
    };

    use super::*;

    fn reduced_terms(plan: &YoungProjectorPlan) -> Vec<(Vec<usize>, Rational)> {
        plan.terms()
            .iter()
            .map(|term| (term.selector().to_vec(), term.weight().clone()))
            .collect()
    }

    fn brute_force_terms(tableau: &YoungTableau) -> Vec<(Vec<usize>, Rational)> {
        let preflight = YoungProjectorPreflight::new(tableau).unwrap();
        let rows = tableau.rows().map(<[usize]>::to_vec).collect::<Vec<_>>();
        let columns = tableau.columns().collect::<Vec<_>>();
        let row_group = YoungProjectorPlan::subgroup(
            tableau.rank(),
            &rows,
            preflight.row_action_count(),
            YoungProjectorQuantity::RowActionCount,
        )
        .unwrap();
        let column_group = YoungProjectorPlan::subgroup(
            tableau.rank(),
            &columns,
            preflight.column_action_count(),
            YoungProjectorQuantity::ColumnActionCount,
        )
        .unwrap();
        let mut coefficients = BTreeMap::<Vec<usize>, Integer>::new();

        for row in &row_group {
            for column in &column_group {
                let action = column.compose(row);
                let (representative, structural_sign) = column_group
                    .iter()
                    .map(|structural| {
                        let candidate = action.compose(structural);
                        (candidate.map().to_vec(), structural.sign())
                    })
                    .min_by(|(left, _), (right, _)| left.cmp(right))
                    .unwrap();
                let contribution = Integer::from((column.sign() * structural_sign) as isize);
                *coefficients
                    .entry(representative)
                    .or_insert_with(|| 0.into()) += contribution;
            }
        }

        let denominator = Integer::from(preflight.hook_product());
        coefficients
            .into_iter()
            .filter(|(_, coefficient)| coefficient != &0)
            .map(|(selector, coefficient)| {
                (selector, Rational::from((coefficient, denominator.clone())))
            })
            .collect()
    }

    #[test]
    fn sorted_coset_reduction_matches_full_right_action() {
        for tableau in [
            YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap(),
            YoungTableau::new(vec![3, 1], vec![0, 2, 3, 1]).unwrap(),
            YoungTableau::new(vec![2, 1, 1], vec![0, 3, 1, 2]).unwrap(),
            YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap(),
        ] {
            let plan = YoungProjectorPlan::new(&tableau).unwrap();
            assert_eq!(reduced_terms(&plan), brute_force_terms(&tableau));
        }
    }

    #[test]
    fn hook_projector_has_three_signed_cosets() {
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let plan = YoungProjectorPlan::new(&tableau).unwrap();

        assert_eq!(plan.hook_product(), 3);
        assert_eq!(
            reduced_terms(&plan),
            [
                (vec![0, 1, 2], Rational::from((2, 3))),
                (vec![0, 2, 1], Rational::from((1, 3))),
                (vec![1, 2, 0], Rational::from((-1, 3))),
            ]
        );
    }

    #[test]
    fn riemann_projector_has_six_signed_cosets() {
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let plan = YoungProjectorPlan::new(&tableau).unwrap();

        assert_eq!(plan.hook_product(), 12);
        assert_eq!(
            reduced_terms(&plan),
            [
                (vec![0, 1, 2, 3], Rational::from((4, 12))),
                (vec![0, 2, 1, 3], Rational::from((2, 12))),
                (vec![0, 3, 1, 2], Rational::from((-2, 12))),
                (vec![1, 2, 0, 3], Rational::from((-2, 12))),
                (vec![1, 3, 0, 2], Rational::from((2, 12))),
                (vec![2, 3, 0, 1], Rational::from((4, 12))),
            ]
        );
    }

    #[test]
    fn riemann_carrier_combines_the_six_reference_cosets_into_three() {
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let plan = YoungProjectorPlan::new(&tableau).unwrap();
        let symmetry = YoungCarrierSymmetry::new(&tableau);
        let mut combined = BTreeMap::<Vec<usize>, Rational>::new();

        for term in plan.terms() {
            let (source_for_target, sign) = symmetry.canonicalize(term.selector()).unwrap();
            let representative = source_for_target
                .into_iter()
                .map(|source| term.selector()[source])
                .collect::<Vec<_>>();
            let coefficient = if sign < 0 {
                -term.weight().clone()
            } else {
                term.weight().clone()
            };
            combined
                .entry(representative)
                .and_modify(|existing| *existing += &coefficient)
                .or_insert(coefficient);
        }

        assert_eq!(
            combined.into_iter().collect::<Vec<_>>(),
            [
                (vec![0, 1, 2, 3], Rational::from((2, 3))),
                (vec![0, 2, 1, 3], Rational::from((1, 3))),
                (vec![0, 3, 1, 2], Rational::from((-1, 3))),
            ]
        );
    }

    #[test]
    fn carrier_normal_form_exchanges_only_whole_equal_height_columns_unsigned() {
        let tableau = YoungTableau::new(vec![2, 2], vec![2, 0, 3, 1]).unwrap();
        let symmetry = YoungCarrierSymmetry::new(&tableau);
        let direct = [2, 3, 0, 1];
        let exchanged = [0, 1, 2, 3];

        let normalize = |keys: &[usize]| {
            symmetry.canonicalize(keys).map(|(sources, sign)| {
                (
                    sources
                        .into_iter()
                        .map(|source| keys[source])
                        .collect::<Vec<_>>(),
                    sign,
                )
            })
        };
        assert_eq!(normalize(&direct), Some((direct.to_vec(), 1)));
        assert_eq!(normalize(&exchanged), Some((direct.to_vec(), 1)));
        assert_eq!(normalize(&[2, 3, 1, 0]), Some((direct.to_vec(), -1)));
        assert_eq!(normalize(&[2, 3, 0, 0]), None);
    }

    #[test]
    fn preflight_reports_action_overflow_before_allocation() {
        let tableau = YoungTableau::canonical(vec![usize::BITS as usize]).unwrap();
        assert!(matches!(
            YoungProjectorPreflight::new(&tableau),
            Err(YoungProjectorPlanError::ArithmeticOverflow {
                quantity: YoungProjectorQuantity::RowActionCount,
                ..
            })
        ));
    }

    fn contains_young_column_carrier(value: AtomView<'_>) -> bool {
        matches!(value, AtomView::Fun(function)
            if YoungColumnCarrier::declaration(function).unwrap().is_some())
            || value.children().any(contains_young_column_carrier)
    }

    fn direct_carrier_atom(value: Atom, tableau: &YoungTableau) -> Atom {
        let (input_negative, value) = match value.as_view() {
            AtomView::Fun(_) => (false, value),
            AtomView::Mul(product) => {
                let mut factors = product.iter();
                let (Some(first), Some(second), None) =
                    (factors.next(), factors.next(), factors.next())
                else {
                    panic!("a signed carrier test tensor must have two factors");
                };
                match (first, second) {
                    (coefficient @ AtomView::Num(_), function @ AtomView::Fun(_))
                    | (function @ AtomView::Fun(_), coefficient @ AtomView::Num(_))
                        if Rational::try_from(coefficient).ok() == Some(Rational::from(-1)) =>
                    {
                        (true, function.to_owned())
                    }
                    _ => panic!("a carrier test tensor may only have a minus sign"),
                }
            }
            _ => panic!("a carrier test tensor must be a function"),
        };
        let AtomView::Fun(function) = value.as_view() else {
            unreachable!("the carrier test sign was stripped")
        };
        let slots = function
            .iter()
            .map(Slot::<LibraryRep, AbstractIndex>::try_from)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let projected = SymbolicTensor {
            structure: OrderedStructure::new(slots).structure,
            is_metric: false,
            is_composite: false,
            expression: value,
        };
        let carrier = YoungColumnCarrier::tensor(projected, tableau).unwrap();
        let carrier_structure = carrier.structure.clone();
        let carrier = if input_negative {
            -carrier.expression
        } else {
            carrier.expression
        };
        let reparsed = CanonicalPolicyNet::<AbstractIndex>::parse(carrier.clone()).unwrap();
        assert_eq!(
            reparsed.network().store.tensors[0].structure,
            carrier_structure
        );
        carrier
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap()
    }

    #[test]
    fn carrier_layout_promotes_ordered_bundles_without_duplicate_argument_sites() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let tensor = tensor_symbol!(young_carrier_layout, young_tableau = tableau.clone());
        let [a, b, c, d] = [
            slot!(rep, young_carrier_layout_a),
            slot!(rep, young_carrier_layout_b),
            slot!(rep, young_carrier_layout_c),
            slot!(rep, young_carrier_layout_d),
        ]
        .map(|slot| slot.to_atom());
        let carrier = direct_carrier_atom(function!(tensor, a, b, c, d), &tableau);
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(carrier).unwrap();
        let layout = super::super::TensorLayout::scan(&policy.network().store.tensors[0]).unwrap();

        assert_eq!(layout.young_columns, [vec![0, 1], vec![2, 3]]);
        assert_eq!(
            layout.exchangeable_young_column_blocks(),
            BTreeMap::from([(2, vec![0, 1])])
        );
        assert_eq!(
            layout.antisymmetric_groups().collect::<Vec<_>>(),
            [
                (super::super::GroupKey::YoungColumn(0), true),
                (super::super::GroupKey::YoungColumn(1), true),
            ]
        );
    }

    #[test]
    fn fixed_carrier_preserves_declared_head_color_order() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let left = tensor_symbol!(young_carrier_color_left, young_tableau = tableau.clone());
        let right = tensor_symbol!(young_carrier_color_right, young_tableau = tableau.clone());
        let arguments = [
            slot!(rep, young_carrier_color_a),
            slot!(rep, young_carrier_color_b),
            slot!(rep, young_carrier_color_c),
            slot!(rep, young_carrier_color_d),
        ]
        .map(|slot| slot.to_atom());
        let component = |head| {
            let slots = arguments
                .iter()
                .map(|argument| {
                    Slot::<LibraryRep, AbstractIndex>::try_from(argument.as_view()).unwrap()
                })
                .collect::<Vec<_>>();
            SymbolicTensor {
                structure: OrderedStructure::new(slots).structure,
                is_metric: false,
                is_composite: false,
                expression: arguments
                    .iter()
                    .fold(FunctionBuilder::new(head), |builder, value| {
                        builder.add_arg(value.clone())
                    })
                    .finish(),
            }
        };
        let left = component(left);
        let right = component(right);
        // Equal slot layouts make this exactly the original SemanticSymbolKey order.
        let declared_order = super::super::TensorLayout::scan(&left)
            .unwrap()
            .color
            .cmp(&super::super::TensorLayout::scan(&right).unwrap().color);
        let left =
            super::super::TensorLayout::scan(&YoungColumnCarrier::tensor(left, &tableau).unwrap())
                .unwrap();
        let right =
            super::super::TensorLayout::scan(&YoungColumnCarrier::tensor(right, &tableau).unwrap())
                .unwrap();

        assert_eq!(left.head, right.head);
        assert_ne!(left.color, right.color);
        assert_eq!(left.color.cmp(&right.color), declared_order);
    }

    #[test]
    fn carrier_three_factor_network_is_dummy_invariant() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let tensor = tensor_symbol!(young_carrier_three_factor, young_tableau = tableau.clone());
        let first = [
            slot!(rep, young_carrier_three_factor_a),
            slot!(rep, young_carrier_three_factor_b),
            slot!(rep, young_carrier_three_factor_c),
            slot!(rep, young_carrier_three_factor_d),
            slot!(rep, young_carrier_three_factor_e),
            slot!(rep, young_carrier_three_factor_f),
        ]
        .map(|slot| slot.to_atom());
        let second = [
            slot!(rep, young_carrier_three_factor_i),
            slot!(rep, young_carrier_three_factor_j),
            slot!(rep, young_carrier_three_factor_k),
            slot!(rep, young_carrier_three_factor_l),
            slot!(rep, young_carrier_three_factor_m),
            slot!(rep, young_carrier_three_factor_n),
        ]
        .map(|slot| slot.to_atom());
        let build = |slots: &[Atom; 6]| {
            [[0, 1, 2, 3], [0, 1, 4, 5], [2, 3, 4, 5]]
                .into_iter()
                .map(|positions| {
                    direct_carrier_atom(
                        positions
                            .into_iter()
                            .fold(FunctionBuilder::new(tensor), |builder, position| {
                                builder.add_arg(slots[position].clone())
                            })
                            .finish(),
                        &tableau,
                    )
                })
                .fold(Atom::num(1), |product, factor| product * factor)
        };

        let first = build(&first)
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        let second = build(&second)
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn malformed_reserved_carrier_is_not_an_ordinary_tensor() {
        let tableau = YoungTableau::canonical(vec![2, 1]).unwrap();
        let tensor = tensor_symbol!(young_carrier_malformed, young_tableau = tableau);
        let malformed = FunctionBuilder::new(YoungColumnCarrier::symbol())
            .add_arg(Atom::var(tensor))
            .finish();
        let AtomView::Fun(function) = malformed.as_view() else {
            panic!("a carrier with its opaque head is a function")
        };

        assert!(matches!(
            YoungColumnCarrier::declaration(function),
            Err(CanonicalizationError::Projection(reason))
                if reason.contains("column arguments")
        ));
    }

    #[test]
    fn structural_carrier_round_trip_preserves_column_sign_and_zero() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![2, 0, 1]).unwrap();
        let tensor = tensor_symbol!(young_carrier_round_trip, young_tableau = tableau.clone());
        let [a, b, c] = [
            slot!(rep, young_carrier_round_trip_a),
            slot!(rep, young_carrier_round_trip_b),
            slot!(rep, young_carrier_round_trip_c),
        ]
        .map(|slot| slot.to_atom());
        let component = |first: &Atom, second: &Atom, third: &Atom| {
            function!(tensor, first.clone(), second.clone(), third.clone())
        };

        let encoded = direct_carrier_atom(component(&a, &b, &c), &tableau);
        let swapped = direct_carrier_atom(component(&a, &c, &b), &tableau);
        assert_eq!(swapped, -encoded.clone());
        assert!(
            component(&a, &b, &b)
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
                .is_zero()
        );

        let decoded = YoungColumnCarrier::decode::<AbstractIndex>(encoded.clone()).unwrap();
        assert!(!contains_young_column_carrier(decoded.as_view()));
        assert_eq!(direct_carrier_atom(decoded, &tableau), encoded);
    }

    #[test]
    fn distinct_declared_heads_use_the_carrier_fast_path_without_leaking_it() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let left = tensor_symbol!(young_carrier_distinct_left, young_tableau = tableau);
        let right = tensor_symbol!(young_carrier_distinct_right, young_tableau = tableau);
        let [a, b, c, d] = [
            slot!(rep, young_carrier_distinct_a),
            slot!(rep, young_carrier_distinct_b),
            slot!(rep, young_carrier_distinct_c),
            slot!(rep, young_carrier_distinct_d),
        ]
        .map(|slot| slot.to_atom());
        let component = |head, arguments: [&Atom; 4]| {
            arguments
                .into_iter()
                .fold(FunctionBuilder::new(head), |builder, argument| {
                    builder.add_arg(argument.clone())
                })
                .finish()
        };
        let input = component(left, [&a, &b, &c, &d]) * component(right, [&a, &c, &b, &d]);

        reset_factored_graph_outcomes(left);
        reset_factored_graph_outcomes(right);
        reset_young_route_outcomes(left);
        reset_young_route_outcomes(right);
        let canonical = input
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();

        assert!(!canonical.is_zero());
        assert!(!contains_young_column_carrier(canonical.as_view()));
        assert_eq!(factored_graph_outcomes(left), (1, 0));
        assert_eq!(factored_graph_outcomes(right), (1, 0));
        assert_eq!(
            young_route_outcomes(left),
            YoungRouteOutcomes {
                carrier: 1,
                composite: 0,
                staged: 0,
            }
        );
        assert_eq!(young_route_outcomes(right), young_route_outcomes(left));
        assert_eq!(
            canonical
                .clone()
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            canonical
        );
    }

    #[test]
    fn carrier_decoration_limit_routes_to_the_exact_staged_projector() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let tensor = tensor_symbol!(
            young_carrier_decoration_limit,
            young_tableau = tableau.clone()
        );
        let arguments = [
            slot!(rep, young_carrier_decoration_limit_a),
            slot!(rep, young_carrier_decoration_limit_b),
        ]
        .map(|slot| slot.to_atom());
        let component_arguments = [
            arguments[0].clone(),
            arguments[1].clone(),
            arguments[0].clone(),
            arguments[1].clone(),
        ];
        let input = component_arguments
            .iter()
            .cloned()
            .fold(FunctionBuilder::new(tensor), FunctionBuilder::add_arg)
            .finish();
        let expected = YoungProjector::new(tableau)
            .project(tensor, &component_arguments)
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();

        reset_young_route_outcomes(tensor);
        let canonical = super::super::projection::with_decoration_orbit_limit(1, || {
            input.try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        })
        .unwrap();

        assert_eq!(
            young_route_outcomes(tensor),
            YoungRouteOutcomes {
                carrier: 0,
                composite: 0,
                staged: 1,
            }
        );
        assert_eq!(canonical, expected);
        assert!(!contains_young_column_carrier(canonical.as_view()));
        assert_eq!(
            canonical
                .clone()
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            canonical
        );
    }

    #[test]
    fn carrier_encoding_leaves_valid_and_malformed_opaque_young_siblings_unchanged() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let exposed = tensor_symbol!(
            young_carrier_opaque_exposed,
            young_tableau = tableau.clone()
        );
        let hidden = tensor_symbol!(young_carrier_opaque_hidden, young_tableau = tableau.clone());
        let wrapper = symbol!("young_carrier_opaque_wrapper");
        let [a, b, c] = [
            slot!(rep, young_carrier_opaque_a),
            slot!(rep, young_carrier_opaque_b),
            slot!(rep, young_carrier_opaque_c),
        ]
        .map(|slot| slot.to_atom());
        let exposed = function!(exposed, a.clone(), b.clone(), c.clone());
        let exposed_canonical = exposed
            .clone()
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        let valid_hidden = function!(hidden, c.clone(), b.clone(), a.clone());
        let encoded_hidden = direct_carrier_atom(valid_hidden.clone(), &tableau);
        let hidden_calls = [valid_hidden, function!(hidden, c, b), encoded_hidden];

        for hidden_call in hidden_calls {
            let opaque = function!(wrapper, bracket!(hidden_call));
            let canonical = (exposed.clone() * opaque.clone())
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap();
            assert_eq!(canonical, exposed_canonical.clone() * opaque);
        }
    }

    #[test]
    fn declared_hook_obeys_bianchi_and_is_idempotent() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let tensor = tensor_symbol!(young_straightener_hook, young_tableau = tableau);
        let [a, b, c] = [
            slot!(rep, young_straightener_hook_a),
            slot!(rep, young_straightener_hook_b),
            slot!(rep, young_straightener_hook_c),
        ]
        .map(|slot| slot.to_atom());
        let component = |left: &Atom, middle: &Atom, right: &Atom| {
            function!(tensor, left.clone(), middle.clone(), right.clone())
        };
        let bianchi = component(&a, &b, &c) + component(&b, &c, &a) + component(&c, &a, &b);
        assert!(
            bianchi
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
                .is_zero()
        );

        let once = component(&a, &b, &c)
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(
            once.try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            once
        );
    }

    #[test]
    fn custom_normalized_young_heads_use_the_staged_path() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let tensor = tensor_symbol!(
            young_straightener_custom_normalized,
            young_tableau = tableau,
            norm = |_input, _out| {}
        );
        assert!(tensor.get_normalization_function().is_some());
        let [a, b, c] = [
            slot!(rep, young_straightener_custom_a),
            slot!(rep, young_straightener_custom_b),
            slot!(rep, young_straightener_custom_c),
        ]
        .map(|slot| slot.to_atom());
        let value = function!(tensor, a, b, c);
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(value.clone()).unwrap();
        assert!(
            YoungNetworkFold::new(policy.network())
                .scope()
                .unwrap()
                .contains_custom_normalizer
        );

        reset_factored_graph_outcomes(tensor);
        let once = value
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(factored_graph_outcomes(tensor), (0, 0));
        assert_eq!(
            once.try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            once
        );
    }

    #[test]
    fn custom_normalized_wrappers_around_young_use_the_staged_path() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let tensor = tensor_symbol!(
            young_straightener_custom_wrapper_tensor,
            young_tableau = tableau
        );
        let wrapper = symbol!(
            "young_straightener_custom_wrapper",
            norm = |_input, _out| {}
        );
        assert!(wrapper.get_normalization_function().is_some());
        let [a, b, c] = [
            slot!(rep, young_straightener_custom_wrapper_a),
            slot!(rep, young_straightener_custom_wrapper_b),
            slot!(rep, young_straightener_custom_wrapper_c),
        ]
        .map(|slot| slot.to_atom());
        let value = function!(wrapper, function!(tensor, a, b, c));
        let policy = CanonicalPolicyNet::<AbstractIndex>::parse(value.clone()).unwrap();
        assert!(
            YoungNetworkFold::new(policy.network())
                .scope()
                .unwrap()
                .contains_custom_normalizer
        );

        reset_factored_graph_outcomes(tensor);
        let once = value
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert_eq!(factored_graph_outcomes(tensor), (0, 0));
        assert_eq!(
            once.try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            once
        );
    }

    #[test]
    fn custom_normalized_siblings_of_young_use_the_staged_path() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let young = tensor_symbol!(
            young_straightener_custom_sibling_young,
            young_tableau = tableau
        );
        let normalized_tensor = tensor_symbol!(
            young_straightener_custom_sibling_tensor,
            norm = |_input, _out| {}
        );
        let ordinary_tensor = tensor_symbol!(young_straightener_custom_sibling_ordinary);
        let wrapper = symbol!(
            "young_straightener_custom_sibling_wrapper",
            norm = |_input, _out| {}
        );
        let [a, b, c, x] = [
            slot!(rep, young_straightener_custom_sibling_a),
            slot!(rep, young_straightener_custom_sibling_b),
            slot!(rep, young_straightener_custom_sibling_c),
            slot!(rep, young_straightener_custom_sibling_x),
        ]
        .map(|slot| slot.to_atom());
        let component = function!(young, a, b, c);
        let siblings = [
            function!(normalized_tensor, x.clone()),
            function!(wrapper, function!(ordinary_tensor, x)),
        ];

        for sibling in siblings {
            let value = component.clone() * sibling;
            let policy = CanonicalPolicyNet::<AbstractIndex>::parse(value.clone()).unwrap();
            assert!(
                YoungNetworkFold::new(policy.network())
                    .scope()
                    .unwrap()
                    .contains_custom_normalizer
            );

            reset_factored_graph_outcomes(young);
            let once = value
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap();
            assert_eq!(factored_graph_outcomes(young), (0, 0));
            assert_eq!(
                once.try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                    .unwrap(),
                once
            );
        }
    }

    #[test]
    fn declared_riemann_obeys_pair_and_bianchi_identities() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let tensor = tensor_symbol!(young_straightener_riemann, young_tableau = tableau);
        let [a, b, c, d] = [
            slot!(rep, young_straightener_riemann_a),
            slot!(rep, young_straightener_riemann_b),
            slot!(rep, young_straightener_riemann_c),
            slot!(rep, young_straightener_riemann_d),
        ]
        .map(|slot| slot.to_atom());
        let component = |first: &Atom, second: &Atom, third: &Atom, fourth: &Atom| {
            function!(
                tensor,
                first.clone(),
                second.clone(),
                third.clone(),
                fourth.clone()
            )
        };
        let canonize = |expression: Atom| {
            expression
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
        };
        let canonical = component(&a, &b, &c, &d);

        assert!(canonize(canonical.clone() + component(&b, &a, &c, &d)).is_zero());
        assert!(canonize(canonical.clone() - component(&c, &d, &a, &b)).is_zero());
        assert!(
            canonize(canonical + component(&a, &c, &d, &b) + component(&a, &d, &b, &c)).is_zero()
        );
        assert!(canonize(component(&a, &a, &c, &d)).is_zero());
    }

    #[test]
    fn declared_riemann_straightens_quadratic_identity() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let tensor = tensor_symbol!(
            young_straightener_riemann_quadratic,
            young_tableau = tableau
        );
        let [a, b, c, d] = [
            slot!(rep, young_straightener_riemann_quadratic_a),
            slot!(rep, young_straightener_riemann_quadratic_b),
            slot!(rep, young_straightener_riemann_quadratic_c),
            slot!(rep, young_straightener_riemann_quadratic_d),
        ]
        .map(|slot| slot.to_atom());
        let component = |first: &Atom, second: &Atom, third: &Atom, fourth: &Atom| {
            function!(
                tensor,
                first.clone(),
                second.clone(),
                third.clone(),
                fourth.clone()
            )
        };
        let direct = component(&a, &b, &c, &d);
        let crossed = component(&a, &c, &b, &d);
        let identity = Atom::num(2) * direct.clone() * crossed - direct.clone() * direct;

        reset_factored_graph_outcomes(tensor);
        assert!(
            identity
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
                .is_zero()
        );
        // Symbolica rewrites the repeated `direct` product to a Power, so this
        // expression takes the staged Power path without a carrier attempt.
        assert_eq!(factored_graph_outcomes(tensor), (0, 0));
    }

    #[test]
    fn negative_power_of_projected_scalar_contraction_remains_supported() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let tensor = tensor_symbol!(young_straightener_inverse, young_tableau = tableau);
        let [a, b, c] = [
            slot!(rep, young_straightener_inverse_a),
            slot!(rep, young_straightener_inverse_b),
            slot!(rep, young_straightener_inverse_c),
        ]
        .map(|slot| slot.to_atom());
        let component = function!(tensor, a, b, c);
        let inverse = (component.clone() * component).pow(Atom::num(-1));
        let once = inverse
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();

        assert!(!once.is_zero());
        assert_eq!(
            once.try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            once
        );
    }

    #[test]
    fn internally_contracted_young_powers_have_independent_copy_scopes() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let tensor = tensor_symbol!(young_straightener_internal_power, young_tableau = tableau);
        let [a, b, c, d] = [
            slot!(rep, young_straightener_internal_power_a),
            slot!(rep, young_straightener_internal_power_b),
            slot!(rep, young_straightener_internal_power_c),
            slot!(rep, young_straightener_internal_power_d),
        ]
        .map(|slot| slot.to_atom());
        let component = |first: &Atom, second: &Atom, third: &Atom, fourth: &Atom| {
            function!(
                tensor,
                first.clone(),
                second.clone(),
                third.clone(),
                fourth.clone()
            )
        };
        let canonize = |expression: Atom| {
            expression
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
        };

        for exponent in [2, -2] {
            let canonical = canonize(component(&a, &b, &a, &b).pow(Atom::num(exponent)));
            let renamed = canonize(component(&c, &d, &c, &d).pow(Atom::num(exponent)));
            assert!(!canonical.is_zero());
            assert_eq!(canonical, renamed);
            assert_eq!(canonize(canonical.clone()), canonical);
        }

        let repeated_column = component(&a, &a, &c, &d).pow(Atom::num(-2));
        assert_eq!(canonize(repeated_column), Atom::Zero.pow(Atom::num(-2)));
    }

    #[test]
    fn young_projection_does_not_distribute_an_unrelated_sum() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let tensor = tensor_symbol!(young_straightener_factored_sum, young_tableau = tableau);
        let [a, b, c] = [
            slot!(rep, young_straightener_factored_sum_a),
            slot!(rep, young_straightener_factored_sum_b),
            slot!(rep, young_straightener_factored_sum_c),
        ]
        .map(|slot| slot.to_atom());
        let scalar_sum = Atom::var(symbol!("young_straightener_scalar_left"))
            + Atom::var(symbol!("young_straightener_scalar_right"));
        let input = function!(tensor, a, b, c) * scalar_sum;
        let canonical = input
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();

        assert!(canonical.nterms() < canonical.expand().nterms());
    }

    #[test]
    fn power_of_an_ordinary_sum_of_young_tensors_stays_factored_and_stable() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let left = tensor_symbol!(young_straightener_power_sum_left, young_tableau = tableau);
        let right = tensor_symbol!(young_straightener_power_sum_right, young_tableau = tableau);
        let intrinsic = tensor_symbol!(young_straightener_power_sum_intrinsic; Antisymmetric);
        let symmetric = tensor_symbol!(young_straightener_power_sum_symmetric; Symmetric);
        let [a, b, c, d, e, f] = [
            slot!(rep, young_straightener_power_sum_a),
            slot!(rep, young_straightener_power_sum_b),
            slot!(rep, young_straightener_power_sum_c),
            slot!(rep, young_straightener_power_sum_d),
            slot!(rep, young_straightener_power_sum_e),
            slot!(rep, young_straightener_power_sum_f),
        ]
        .map(|slot| slot.to_atom());
        let component = |head, first: &Atom, second: &Atom, third: &Atom| {
            function!(head, first.clone(), second.clone(), third.clone())
        };
        for companion in [right, intrinsic, symmetric] {
            for exponent in [2, -2] {
                let input = (component(left, &a, &b, &c) + component(companion, &a, &b, &c))
                    .pow(Atom::num(exponent));
                let renamed = (component(left, &d, &e, &f) + component(companion, &d, &e, &f))
                    .pow(Atom::num(exponent))
                    .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                    .unwrap();
                let swapped = (component(left, &b, &a, &c) + component(companion, &b, &a, &c))
                    .pow(Atom::num(exponent))
                    .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                    .unwrap();
                let canonical = input
                    .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                    .unwrap();

                if exponent > 0 {
                    assert!(canonical.nterms() < canonical.expand().nterms());
                }
                assert!(!canonical.is_zero());
                assert_eq!(canonical, renamed);
                assert_eq!(canonical, swapped);
                assert_eq!(
                    canonical
                        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                        .unwrap(),
                    canonical
                );
            }
        }
    }

    #[test]
    fn nonlinear_function_retains_a_zero_young_argument() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 1], vec![0, 2, 1]).unwrap();
        let tensor = tensor_symbol!(young_straightener_nonlinear_hook, young_tableau = tableau);
        let wrapper = symbol!("young_straightener_nonlinear_wrapper");
        let [a, b, c] = [
            slot!(rep, young_straightener_nonlinear_a),
            slot!(rep, young_straightener_nonlinear_b),
            slot!(rep, young_straightener_nonlinear_c),
        ]
        .map(|slot| slot.to_atom());
        let component = |left: &Atom, middle: &Atom, right: &Atom| {
            function!(tensor, left.clone(), middle.clone(), right.clone())
        };
        let bianchi = component(&a, &b, &c) + component(&b, &c, &a) + component(&c, &a, &b);
        let canonical = function!(wrapper, bianchi)
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();

        assert_eq!(canonical, function!(wrapper, Atom::Zero));
    }

    #[test]
    fn young_projector_and_live_term_limits_are_typed() {
        let rep = test_initialize().mink4;
        let oversized_tableau = YoungTableau::canonical(vec![7, 1]).unwrap();
        let oversized = tensor_symbol!(
            young_straightener_action_limit,
            young_tableau = oversized_tableau
        );
        let oversized_expression = (0..8)
            .map(|index| slot!(rep, AbstractIndex::Dummy(index)).to_atom())
            .fold(FunctionBuilder::new(oversized), FunctionBuilder::add_arg)
            .finish();
        assert_eq!(
            oversized_expression
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap_err(),
            CanonicalizationError::YoungProjectorSizeLimit {
                head: oversized,
                shape: vec![7, 1],
                requested_actions: 10_080,
                action_limit: MAX_YOUNG_PROJECTOR_ACTIONS,
            }
        );

        let product_tableau = YoungTableau::canonical(vec![4, 2]).unwrap();
        let left = tensor_symbol!(
            young_straightener_product_limit_left,
            young_tableau = product_tableau.clone()
        );
        let right = tensor_symbol!(
            young_straightener_product_limit_right,
            young_tableau = product_tableau
        );
        let left = (0..6)
            .map(|index| slot!(rep, AbstractIndex::Dummy(index)).to_atom())
            .fold(FunctionBuilder::new(left), FunctionBuilder::add_arg)
            .finish();
        let right = (6..12)
            .map(|index| slot!(rep, AbstractIndex::Dummy(index)).to_atom())
            .fold(FunctionBuilder::new(right), FunctionBuilder::add_arg)
            .finish();
        assert_eq!(
            (left * right)
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap_err(),
            CanonicalizationError::YoungExpansionSizeLimit {
                requested_terms: 19_600,
                term_limit: MAX_YOUNG_EXPANSION_TERMS,
            }
        );

        let sum_tableau = YoungTableau::canonical(vec![4, 2]).unwrap();
        let sum_slots = (0..6)
            .map(|position| {
                let index = if position == 4 { 0 } else { position };
                slot!(rep, AbstractIndex::Dummy(index)).to_atom()
            })
            .collect::<Vec<_>>();
        let oversized_sum = (0..30)
            .map(|term| {
                let head = tensor_symbol!(
                    format!("young_straightener_sum_limit_{term}"),
                    young_tableau = sum_tableau.clone()
                );
                sum_slots
                    .iter()
                    .cloned()
                    .fold(FunctionBuilder::new(head), FunctionBuilder::add_arg)
                    .finish()
            })
            .fold(Atom::Zero, |sum, term| sum + term);
        assert_eq!(
            oversized_sum
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap_err(),
            CanonicalizationError::YoungExpansionSizeLimit {
                requested_terms: 4_200,
                term_limit: MAX_YOUNG_EXPANSION_TERMS,
            }
        );
    }

    #[test]
    fn projected_three_factor_network_is_dummy_invariant_and_stable() {
        let rep = test_initialize().mink4;
        let tableau = YoungTableau::new(vec![2, 2], vec![0, 2, 1, 3]).unwrap();
        let tensor = tensor_symbol!(young_straightener_three_factor, young_tableau = tableau);
        let indices = [
            slot!(rep, young_straightener_three_factor_a),
            slot!(rep, young_straightener_three_factor_b),
            slot!(rep, young_straightener_three_factor_c),
            slot!(rep, young_straightener_three_factor_d),
            slot!(rep, young_straightener_three_factor_e),
            slot!(rep, young_straightener_three_factor_f),
        ]
        .map(|slot| slot.to_atom());
        let renamed = [
            slot!(rep, young_straightener_three_factor_i),
            slot!(rep, young_straightener_three_factor_j),
            slot!(rep, young_straightener_three_factor_k),
            slot!(rep, young_straightener_three_factor_l),
            slot!(rep, young_straightener_three_factor_m),
            slot!(rep, young_straightener_three_factor_n),
        ]
        .map(|slot| slot.to_atom());
        let build = |slots: &[Atom; 6]| {
            function!(
                tensor,
                slots[0].clone(),
                slots[1].clone(),
                slots[2].clone(),
                slots[3].clone()
            ) * function!(
                tensor,
                slots[0].clone(),
                slots[1].clone(),
                slots[4].clone(),
                slots[5].clone()
            ) * function!(
                tensor,
                slots[2].clone(),
                slots[3].clone(),
                slots[4].clone(),
                slots[5].clone()
            )
        };
        let canonical = build(&indices)
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        let renamed = build(&renamed)
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert!(!canonical.is_zero());
        assert_eq!(canonical, renamed);
        assert_eq!(
            canonical
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            canonical
        );
    }
}
