use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use spenso::algebra::complex::Complex;
use symbolica::prelude::{
    AliasedAtom, Atom, AtomCore, AtomView, Rational, Replacement, SingleFloat, Symbol, function,
    symbol,
};

use crate::{
    GammaLoopContext,
    integrands::evaluation::EvaluationMetaData,
    utils::{F, FloatLike, hyperdual_utils::DualOrNot},
    uv::ParametricIntegrandTerm,
};

use super::{EvaluatorStack, GenericEvaluator, GenericEvaluatorFloat};

/// Exact evaluator-local inputs. Only the single-parametric program receives
/// these appended ports; all other programs specialize them after preprocessing.
#[derive(Clone, Debug, Default, PartialEq, Eq, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct ParametricResidueRows {
    pub(crate) parameters: Vec<Atom>,
    pub(crate) rows: Vec<Vec<Rational>>,
    pub(crate) input_columns: Vec<usize>,
    pub(crate) auxiliary_parameters: Vec<Atom>,
    pub(crate) auxiliary_values: Vec<Rational>,
    pub(crate) auxiliary_shared_outputs: usize,
    pub(crate) auxiliary_gathers: Vec<Vec<Option<usize>>>,
}

#[derive(Clone)]
enum ScalarAuxiliaryKind {
    Leaf,
    Scalar {
        definitions: std::collections::BTreeSet<usize>,
        columns: std::collections::BTreeSet<usize>,
    },
    Selector,
    Numerator(usize, Vec<Option<usize>>, std::collections::BTreeSet<usize>),
}

/// Scope and activity are established on the original caller DAG. Numerator
/// bodies are inspected for captures once and remain opaque during extraction.
struct ScalarAuxiliaryAnalysis<'a> {
    rows: &'a ParametricResidueRows,
    definitions: &'a [super::FnMapEntry],
    generated_start: usize,
    physical: std::collections::HashSet<Atom>,
    lookup_keys: std::collections::HashSet<Atom>,
    columns: std::collections::HashMap<Atom, usize>,
    functions: std::collections::HashMap<Atom, usize>,
    families: std::collections::HashMap<(Symbol, Vec<Atom>), usize>,
    tag_counts: std::collections::HashMap<Symbol, usize>,
    boolean: std::collections::HashSet<Atom>,
    closed: std::collections::HashMap<usize, Option<std::collections::BTreeSet<usize>>>,
    kinds: std::collections::HashMap<(usize, Atom), Option<ScalarAuxiliaryKind>>,
}

impl ScalarAuxiliaryAnalysis<'_> {
    fn function(&self, view: AtomView<'_>) -> Option<(usize, Vec<Atom>)> {
        let AtomView::Fun(call) = view else {
            return None;
        };
        let count = *self.tag_counts.get(&call.get_symbol())?;
        let tags = call.iter().take(count).map(|arg| arg.to_owned()).collect();
        let index = *self.families.get(&(call.get_symbol(), tags))?;
        (call.get_nargs() == count + self.definitions[index].args.len()).then(|| {
            (
                index,
                call.iter().skip(count).map(|arg| arg.to_owned()).collect(),
            )
        })
    }

    // Ordinary function calls reset lexical scope. Caller aliases never certify
    // a free symbol captured by a numerator or another ordinary function.
    fn closed_definition(&mut self, index: usize) -> Option<std::collections::BTreeSet<usize>> {
        if let Some(result) = self.closed.get(&index) {
            return result.clone();
        }
        self.closed.insert(index, None); // cycles are not movable
        let formals = self.definitions[index]
            .args
            .iter()
            .cloned()
            .map(Atom::from)
            .collect();
        let body = self.definitions[index].rhs.clone();
        let mut used = std::collections::BTreeSet::from([index]);
        let result = self
            .closed_body(body.as_view(), &formals, &mut used)
            .then_some(used);
        self.closed.insert(index, result.clone());
        result
    }

    fn closed_body(
        &mut self,
        view: AtomView<'_>,
        formals: &std::collections::HashSet<Atom>,
        used: &mut std::collections::BTreeSet<usize>,
    ) -> bool {
        if formals.contains(view.get_data()) || self.physical.contains(view.get_data()) {
            return true;
        }
        if self.columns.contains_key(view.get_data()) {
            return false;
        }
        if let Some((index, args)) = self.function(view) {
            if self.definitions[index].is_alias {
                return false;
            }
            if !args
                .iter()
                .all(|arg| self.closed_body(arg.as_view(), formals, used))
            {
                return false;
            }
            let Some(dependencies) = self.closed_definition(index) else {
                return false;
            };
            used.extend(dependencies);
            return true;
        }
        match view {
            AtomView::Num(_) => true,
            AtomView::Var(var) => [Symbol::E, Symbol::PI].contains(&var.get_symbol()),
            AtomView::Add(sum) => sum.iter().all(|part| self.closed_body(part, formals, used)),
            AtomView::Mul(product) => product
                .iter()
                .all(|part| self.closed_body(part, formals, used)),
            AtomView::Pow(power) => {
                self.closed_body(power.get_base(), formals, used)
                    && self.closed_body(power.get_exp(), formals, used)
            }
            AtomView::Fun(call) if call.get_symbol().is_builtin() => call
                .iter()
                .all(|part| self.closed_body(part, formals, used)),
            _ => false,
        }
    }

    fn scalar_dependencies(
        &mut self,
        view: AtomView<'_>,
        owner: Option<&AliasedAtom>,
        active_aliases: &mut std::collections::HashSet<Atom>,
        active_functions: &mut std::collections::HashSet<usize>,
        used: &mut std::collections::BTreeSet<usize>,
        columns: &mut std::collections::BTreeSet<usize>,
    ) -> Option<(bool, bool)> {
        if self.physical.contains(view.get_data()) {
            return Some((true, false));
        }
        if let Some(body) = owner.and_then(|owner| owner.get_aliases().get(view.get_data())) {
            let key = view.to_owned();
            if !active_aliases.insert(key.clone()) {
                return None;
            }
            let result = self.scalar_dependencies(
                body.as_view(),
                owner,
                active_aliases,
                active_functions,
                used,
                columns,
            );
            active_aliases.remove(&key);
            return result;
        }
        // Row values belong to this caller. Ordinary function bodies below
        // have their own lexical scope and cannot capture these bindings.
        if owner.is_some()
            && let Some(column) = self.columns.get(view.get_data())
        {
            columns.insert(*column);
            return Some((false, false));
        }
        if let Some(&index) = self.functions.get(view.get_data()) {
            if !active_functions.insert(index) {
                return None;
            }
            let body = self.definitions[index].rhs.clone();
            let result = self.scalar_dependencies(
                body.as_view(),
                None,
                active_aliases,
                active_functions,
                used,
                columns,
            );
            active_functions.remove(&index);
            used.insert(index);
            return result;
        }
        if let Some((index, args)) = self.function(view)
            && index < self.generated_start
            && !args.is_empty()
            && !self.definitions[index].is_alias
        {
            // Certify the ordinary body in its own formal/global scope, then
            // inspect actual arguments in this caller. Keep the call intact;
            // generated numerator definitions remain a separate boundary.
            let mut tagged_binding = false;
            for tag in &self.definitions[index].tags {
                tag.visitor(&mut |part| {
                    tagged_binding |= self.columns.contains_key(part.get_data())
                        || owner
                            .is_some_and(|owner| owner.get_aliases().contains_key(part.get_data()));
                    !tagged_binding
                });
            }
            // Role rebinding visits the scalar AST. Metadata must not contain
            // a caller binding that would change this exact function identity.
            if tagged_binding {
                return None;
            }
            let dependencies = self.closed_definition(index)?;
            if dependencies
                .iter()
                .any(|index| *index >= self.generated_start)
            {
                return None;
            }
            let mut physical = false;
            for arg in &args {
                physical |= self
                    .scalar_dependencies(
                        arg.as_view(),
                        owner,
                        active_aliases,
                        active_functions,
                        used,
                        columns,
                    )?
                    .0;
            }
            used.extend(dependencies);
            return Some((physical, true));
        }
        // Keep nested conditionals at the existing selector/gather boundary;
        // the enclosing arithmetic must retain distinct opaque selector ports.
        let mut recurse = |part| {
            self.scalar_dependencies(part, owner, active_aliases, active_functions, used, columns)
        };
        match view {
            AtomView::Num(_) => Some((false, false)),
            AtomView::Var(var) if [Symbol::E, Symbol::PI].contains(&var.get_symbol()) => {
                Some((false, false))
            }
            AtomView::Add(sum) => sum.iter().try_fold((false, true), |(depends, _), part| {
                Some((recurse(part)?.0 || depends, true))
            }),
            AtomView::Mul(product) => product
                .iter()
                .try_fold((false, true), |(depends, _), part| {
                    Some((recurse(part)?.0 || depends, true))
                }),
            AtomView::Pow(power) => Some((
                recurse(power.get_base())?.0 | recurse(power.get_exp())?.0,
                true,
            )),
            AtomView::Fun(call)
                if [
                    Symbol::EXP,
                    Symbol::LOG,
                    Symbol::SIN,
                    Symbol::COS,
                    Symbol::SQRT,
                    Symbol::CONJ,
                    Symbol::ABS,
                ]
                .contains(&call.get_symbol()) =>
            {
                call.iter().try_fold((false, true), |(depends, _), part| {
                    Some((recurse(part)?.0 || depends, true))
                })
            }
            _ => None,
        }
    }

    fn guard(&self, view: AtomView<'_>, owner: &AliasedAtom, row: Option<usize>) -> Option<bool> {
        if self.physical.contains(view.get_data())
            || owner.get_aliases().contains_key(view.get_data())
            || self.lookup_keys.contains(view.get_data())
        {
            return None;
        }
        if self.boolean.contains(view.get_data()) {
            return Some(
                row.is_none_or(|row| !self.rows.rows[row][self.columns[view.get_data()]].is_zero()),
            );
        }
        // Exact native row inputs shadow their apparent syntax, including IF
        // expressions whose stored values are not Boolean.
        if self.columns.contains_key(view.get_data()) {
            return None;
        }
        if view.is_zero() {
            return Some(false);
        }
        if view.to_owned() == Atom::one() {
            return Some(true);
        }
        let AtomView::Fun(call) = view else {
            return None;
        };
        if call.get_symbol() != Symbol::IF || call.get_nargs() != 3 {
            return None;
        }
        let mut args = call.iter();
        let condition = self.guard(args.next().unwrap(), owner, row)?;
        let yes = args.next().unwrap();
        let no = args.next().unwrap();
        if row.is_none() {
            self.guard(yes, owner, None)?;
            self.guard(no, owner, None)?;
            Some(true)
        } else {
            self.guard(if condition { yes } else { no }, owner, row)
        }
    }

    fn kind(
        &mut self,
        view: AtomView<'_>,
        owner_index: usize,
        owner: &AliasedAtom,
    ) -> Option<ScalarAuxiliaryKind> {
        if self.physical.contains(view.get_data())
            || matches!(view, AtomView::Num(_))
            || matches!(view,AtomView::Var(var) if [Symbol::E,Symbol::PI].contains(&var.get_symbol()))
        {
            return Some(ScalarAuxiliaryKind::Leaf);
        }
        let key = (owner_index, view.to_owned());
        if let Some(kind) = self.kinds.get(&key) {
            return kind.clone();
        }
        self.kinds.insert(key.clone(), None);
        // Preserve direct selector gathering before considering compound scalar
        // leaves. Aliases inherit only scalar/selector kinds in their own scope;
        // a numerator call behind an alias still keeps its original call port.
        let result = if let Some(body) = owner.get_aliases().get(view.get_data()) {
            match self.kind(body.as_view(), owner_index, owner) {
                Some(
                    kind @ (ScalarAuxiliaryKind::Leaf
                    | ScalarAuxiliaryKind::Scalar { .. }
                    | ScalarAuxiliaryKind::Selector),
                ) => Some(kind),
                _ => None,
            }
        } else if let AtomView::Fun(call) = view
            && call.get_symbol() == Symbol::IF
            && !self.columns.contains_key(view.get_data())
        {
            if call.get_nargs() == 3 {
                let mut args = call.iter();
                let condition = args.next().unwrap();
                let yes = args.next().unwrap();
                let no = args.next().unwrap();
                if self.guard(view, owner, None).is_some()
                    || (self.guard(condition, owner, None).is_some()
                        && [yes, no].into_iter().all(|part| {
                            matches!(
                                self.kind(part, owner_index, owner),
                                Some(
                                    ScalarAuxiliaryKind::Leaf
                                        | ScalarAuxiliaryKind::Scalar { .. }
                                        | ScalarAuxiliaryKind::Selector
                                )
                            )
                        }))
                {
                    Some(ScalarAuxiliaryKind::Selector)
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            let mut used = std::collections::BTreeSet::new();
            let mut columns = std::collections::BTreeSet::new();
            if let Some((_, arithmetic)) = self.scalar_dependencies(
                view,
                Some(owner),
                &mut Default::default(),
                &mut Default::default(),
                &mut used,
                &mut columns,
            ) {
                Some(if arithmetic || !columns.is_empty() {
                    ScalarAuxiliaryKind::Scalar {
                        definitions: used,
                        columns,
                    }
                } else {
                    ScalarAuxiliaryKind::Leaf
                })
            } else if let Some((index, args)) = self.function(view) {
                if index < self.generated_start
                    || self.definitions[index].is_alias
                    || self.definitions[index].inlining
                        == symbolica::evaluate::InliningPolicy::Never
                {
                    None
                } else {
                    let columns = args
                        .iter()
                        .map(|arg| {
                            (!self.physical.contains(arg.as_view().get_data()))
                                .then(|| self.columns.get(arg.as_view().get_data()).copied())
                                .flatten()
                        })
                        .collect::<Vec<_>>();
                    if !columns.iter().any(Option::is_some)
                        || !args.iter().zip(&columns).all(|(arg, column)| {
                            column.is_some() || matches!(arg.as_view(), AtomView::Num(_))
                        })
                    {
                        None
                    } else {
                        self.closed_definition(index)
                            .map(|used| ScalarAuxiliaryKind::Numerator(index, columns, used))
                    }
                }
            } else {
                None
            }
        };
        self.kinds.insert(key, result.clone());
        result
    }

    fn collect(
        &mut self,
        view: AtomView<'_>,
        owner_index: usize,
        owner: &AliasedAtom,
        active: &[usize],
        active_aliases: &mut std::collections::HashSet<Atom>,
        selected: &mut std::collections::BTreeMap<(usize, Atom), std::collections::BTreeSet<usize>>,
    ) {
        if active.is_empty() {
            return;
        }
        match self.kind(view, owner_index, owner) {
            Some(ScalarAuxiliaryKind::Leaf) => return,
            Some(_) => {
                selected
                    .entry((owner_index, view.to_owned()))
                    .or_default()
                    .extend(active.iter().copied());
                return;
            }
            None => {}
        }
        if let Some(body) = owner.get_aliases().get(view.get_data()) {
            let key = view.to_owned();
            if active_aliases.insert(key.clone()) {
                self.collect(
                    body.as_view(),
                    owner_index,
                    owner,
                    active,
                    active_aliases,
                    selected,
                );
                active_aliases.remove(&key);
            }
            return;
        }
        match view {
            AtomView::Add(sum) => {
                for part in sum.iter() {
                    self.collect(part, owner_index, owner, active, active_aliases, selected);
                }
            }
            AtomView::Mul(product) => {
                for part in product.iter() {
                    self.collect(part, owner_index, owner, active, active_aliases, selected);
                }
            }
            AtomView::Pow(power) => {
                self.collect(
                    power.get_base(),
                    owner_index,
                    owner,
                    active,
                    active_aliases,
                    selected,
                );
                self.collect(
                    power.get_exp(),
                    owner_index,
                    owner,
                    active,
                    active_aliases,
                    selected,
                );
            }
            AtomView::Fun(call) if call.get_symbol() == Symbol::IF && call.get_nargs() == 3 => {
                let mut args = call.iter();
                let condition = args.next().unwrap();
                let yes = args.next().unwrap();
                let no = args.next().unwrap();
                if self.guard(condition, owner, None).is_some() {
                    // A nested Boolean selector used as an outer condition can
                    // itself be gathered. The outer lazy branch remains intact.
                    self.collect(
                        condition,
                        owner_index,
                        owner,
                        active,
                        active_aliases,
                        selected,
                    );
                    let (positive, negative): (Vec<_>, Vec<_>) =
                        active.iter().copied().partition(|row| {
                            self.guard(condition, owner, Some(*row))
                                .expect("certified Boolean selector")
                        });
                    self.collect(yes, owner_index, owner, &positive, active_aliases, selected);
                    self.collect(no, owner_index, owner, &negative, active_aliases, selected);
                }
            }
            _ => {} // Unknown predicates, ordinary N calls and actuals stay opaque.
        }
    }

    fn leaf<'a>(
        &self,
        mut view: AtomView<'a>,
        owner_index: usize,
        owner: &'a AliasedAtom,
        row: usize,
    ) -> Atom {
        loop {
            // Exact native inputs stay opaque. Only aliases already certified
            // as pure selectors are followed; scalar aliases keep their DAG.
            if self.physical.contains(view.get_data()) {
                break;
            }
            if let Some(body) = owner.get_aliases().get(view.get_data()) {
                if matches!(
                    self.kinds.get(&(owner_index, view.to_owned())),
                    Some(Some(ScalarAuxiliaryKind::Selector))
                ) {
                    view = body.as_view();
                    continue;
                }
                break;
            }
            if self.columns.contains_key(view.get_data())
                || self.lookup_keys.contains(view.get_data())
            {
                break;
            }
            let AtomView::Fun(call) = view else {
                break;
            };
            if call.get_symbol() != Symbol::IF {
                break;
            }
            let mut args = call.iter();
            let condition = args.next().unwrap();
            let yes = args.next().unwrap();
            let no = args.next().unwrap();
            view = if self
                .guard(condition, owner, Some(row))
                .expect("certified Boolean selector")
            {
                yes
            } else {
                no
            };
        }
        // Boolean selector leaves may themselves be original row ports.
        // Resolving their certified 0/1 values does not rebuild arithmetic or
        // replace a physical lookup/alias, which guard() treats as a boundary.
        self.guard(view, owner, Some(row))
            .map_or_else(|| view.to_owned(), |value| Atom::num(i64::from(value)))
    }
}

impl ParametricResidueRows {
    fn select(values: &[(usize, Atom)], bits: &[Atom]) -> Atom {
        if values.is_empty() {
            return Atom::Zero;
        }
        let Some((bit, remaining)) = bits.split_last() else {
            debug_assert_eq!(values.len(), 1);
            return values[0].1.clone();
        };
        let mask = 1usize << remaining.len();
        let split = values.partition_point(|(index, _)| index & mask == 0);
        let positive = Self::select(&values[split..], remaining);
        let negative = Self::select(&values[..split], remaining);
        if positive == negative {
            return positive;
        }
        // Keep complete carriers at selector leaves. Factoring through a child IF
        // would make its parent arm compound arithmetic, blocking direct gather
        // selection and retaining native branches in the row evaluator.
        Symbol::IF.call_args([bit.clone(), positive, negative])
    }

    pub(super) fn lower(
        integrand: &Atom,
        terms: &[ParametricIntegrandTerm],
    ) -> Result<(Vec<Atom>, Vec<Atom>, Self)> {
        // Physical maps, cuts and threshold localization are complete here.
        // Only now may equal full coefficient rows share their scalar sum;
        // earlier UV owners retain each carrier's native map provenance.
        let term_rows = terms
            .iter()
            .map(|term| {
                let mut grouped = std::collections::BTreeMap::<Vec<Rational>, Vec<&Atom>>::new();
                for row in &term.rows {
                    if row.coefficients.len() != term.numerator.args.len() {
                        return Err(eyre!(
                            "parametric numerator row has {} coefficients for {} formals",
                            row.coefficients.len(),
                            term.numerator.args.len()
                        ));
                    }
                    grouped
                        .entry(row.coefficients.clone())
                        .or_default()
                        .push(&row.carrier);
                }
                Ok(grouped
                    .into_iter()
                    .filter_map(|(coefficients, carriers)| {
                        let carrier = Atom::add_many(carriers);
                        (!carrier.is_zero()).then_some((coefficients, carrier))
                    })
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>>>()?;
        let has_scalar = integrand != &Atom::Zero;
        let row_count = term_rows
            .iter()
            .map(Vec::len)
            .max()
            .unwrap_or(0)
            .max(usize::from(has_scalar));
        // A cut whose terms all vanish still has one executable zero, never an
        // empty reduction or an extra contribution from another source.
        if row_count == 0 {
            return Ok((
                vec![Atom::Zero],
                vec![Atom::one()],
                Self {
                    parameters: Vec::new(),
                    rows: vec![Vec::new()],
                    auxiliary_gathers: vec![Vec::new()],
                    ..Default::default()
                },
            ));
        }
        // Zip independent additive source sums by their local row ordinal.
        // Rows within a family already describe its complete Cartesian product
        // of UV/cograph sources; their coefficients and carriers stay together.
        // Inactive tails repeat that family's final row. Exact column equality
        // permits shared input ports without identifying numerator definitions;
        // constant columns need no port at all.
        let term_columns = terms
            .iter()
            .zip(&term_rows)
            .map(|(term, rows)| {
                if rows.is_empty() {
                    return Vec::new();
                }
                (0..term.numerator.args.len())
                    .map(|column| {
                        (0..row_count)
                            .map(|row| rows[row.min(rows.len() - 1)].0[column].clone())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut columns = term_columns
            .iter()
            .flatten()
            .filter(|column| column.iter().any(|value| value != &column[0]))
            .map(|column| (column.clone(), 0))
            .collect::<std::collections::BTreeMap<_, _>>();
        for (index, port) in columns.values_mut().enumerate() {
            *port = index;
        }
        let coefficient_count = columns.len();
        let row_bits = usize::BITS as usize - row_count.saturating_sub(1).leading_zeros() as usize;
        let coefficient = symbol!("gammalooprs::parametric_residue::coefficient"; Scalar);
        let selector = symbol!("gammalooprs::parametric_residue::selector"; Scalar);
        let parameters = (0..coefficient_count)
            .map(|i| function!(coefficient, i))
            .chain((0..row_bits).map(|i| function!(selector, i)))
            .collect::<Vec<_>>();
        let mut atoms = Vec::with_capacity(terms.len() + usize::from(has_scalar));
        let mut guards = Vec::with_capacity(atoms.capacity());
        for ((term, source_rows), term_columns) in terms.iter().zip(&term_rows).zip(&term_columns) {
            if source_rows.is_empty() {
                continue;
            }
            let replacements = term
                .numerator
                .args
                .iter()
                .zip(term_columns)
                .map(|(formal, column)| {
                    let value = columns.get(column).map_or_else(
                        || Atom::num(column[0].clone()),
                        |port| parameters[*port].clone(),
                    );
                    Replacement::new(Atom::from(formal.clone()).to_pattern(), value)
                })
                .collect::<Vec<_>>();
            let carriers = source_rows
                .iter()
                .enumerate()
                .map(|(row, (_, carrier))| (row, carrier.replace_multiple(&replacements)))
                .collect::<Vec<_>>();
            // Keep the numerator product exposed during tensor preprocessing.
            // Its whole-product guard is attached to the resulting scalar.
            // Inactive rows are excluded there, so a single carrier needs no
            // inner conditional hiding its open tensor interface.
            let carrier = match carriers.as_slice() {
                [(_, carrier)] => carrier.clone(),
                _ => Self::select(&carriers, &parameters[coefficient_count..]),
            };
            atoms.push(term.numerator.lhs.replace_multiple(&replacements) * carrier);
            guards.push(if source_rows.len() == row_count {
                // Unused selector bit patterns are never supplied by the table.
                Atom::one()
            } else {
                let active = (0..source_rows.len())
                    .map(|row| (row, Atom::one()))
                    .collect::<Vec<_>>();
                Self::select(&active, &parameters[coefficient_count..])
            });
        }
        if has_scalar {
            atoms.push(integrand.clone());
            guards.push(Self::select(
                &[(0, Atom::one())],
                &parameters[coefficient_count..],
            ));
        }
        let rows: Vec<Vec<Rational>> = (0..row_count)
            .map(|row| {
                columns
                    .keys()
                    .map(|column| column[row].clone())
                    .chain((0..row_bits).map(|bit| Rational::from(((row >> bit) & 1) as i64)))
                    .collect()
            })
            .collect();
        Ok((
            atoms,
            guards,
            Self {
                input_columns: (0..parameters.len()).collect(),
                auxiliary_gathers: vec![Vec::new(); rows.len()],
                parameters,
                rows,
                ..Default::default()
            },
        ))
    }

    pub(super) fn join_preprocessed(
        &self,
        atoms: Vec<AliasedAtom>,
        guards: &[Atom],
    ) -> AliasedAtom {
        assert_eq!(atoms.len(), guards.len());
        // Another family's inputs can make an inactive numerator singular.
        // Guard the complete scalar product only after its tensors contracted.
        // Sources with the same activity guard execute together. Keep their
        // scalar sum in one lazy branch so evaluator CSE can share its factors.
        let mut groups = std::collections::BTreeMap::<&Atom, Vec<&Atom>>::new();
        for (atom, guard) in atoms.iter().zip(guards) {
            groups.entry(guard).or_default().push(atom.get_root());
        }
        let mut joined =
            AliasedAtom::from(Atom::add_many(groups.into_iter().map(|(guard, roots)| {
                Symbol::IF.call_args([guard.clone(), Atom::add_many(roots), Atom::Zero])
            })));
        for atom in atoms {
            for (alias, body) in atom.into_inner_with_aliases().1 {
                joined.register_alias(alias, body);
            }
        }
        joined.prune();
        joined
    }

    pub(super) fn specialize(
        &self,
        atom: &AliasedAtom,
        row: &[Rational],
        row_id: usize,
    ) -> AliasedAtom {
        let replacements = self
            .parameters
            .iter()
            .zip(row)
            .map(|(parameter, value)| {
                Replacement::new(parameter.to_pattern(), Atom::num(value.clone()))
            })
            .collect::<Vec<_>>();
        EvaluatorStack::specialize_aliases(
            atom.get_root(),
            atom,
            |expression| expression.replace_multiple(&replacements),
            Some(row_id),
        )
    }

    /// Evaluate shared physical scalars and active numerator argument rows once.
    /// Each caller keeps its own port; only auxiliary outputs and dynamic input
    /// values share identities. No numerator body is substituted or expanded.
    pub(super) fn extract_physical_scalars(
        &mut self,
        atoms: &mut [AliasedAtom],
        physical_parameters: &[Atom],
        definitions: &[super::FnMapEntry],
        generated_definition_start: usize,
        port_symbol: Symbol,
    ) -> Option<(Vec<Atom>, Vec<AliasedAtom>, Vec<super::FnMapEntry>)> {
        use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
        self.input_columns = (0..self.parameters.len()).collect();
        self.auxiliary_parameters.clear();
        self.auxiliary_values.clear();
        self.auxiliary_shared_outputs = 0;
        self.auxiliary_gathers = vec![Vec::new(); self.rows.len()];
        let physical = physical_parameters.iter().cloned().collect::<HashSet<_>>();
        let eligible_owners = atoms
            .iter()
            .map(|atom| {
                !physical_parameters
                    .iter()
                    .chain(&self.parameters)
                    .chain(definitions.iter().map(|entry| &entry.lhs))
                    .any(|key| {
                        let mut collision = false;
                        key.visitor(&mut |view| {
                            collision |= atom.get_aliases().contains_key(view.get_data());
                            !collision
                        });
                        collision
                    })
            })
            .collect::<Vec<_>>();
        let mut tag_counts = HashMap::new();
        let mut ambiguous_heads = HashSet::new();
        for definition in definitions {
            if let Some(call) = definition.lhs.as_fun_view()
                && let Some(previous) = tag_counts.insert(call.get_symbol(), definition.tags.len())
                && previous != definition.tags.len()
            {
                ambiguous_heads.insert(call.get_symbol());
            }
        }
        tag_counts.retain(|head, _| !ambiguous_heads.contains(head));
        let mut analysis = ScalarAuxiliaryAnalysis {
            rows: self,
            definitions,
            generated_start: generated_definition_start,
            physical: physical.clone(),
            lookup_keys: physical_parameters
                .iter()
                .chain(definitions.iter().map(|entry| &entry.lhs))
                .cloned()
                .collect(),
            columns: self
                .parameters
                .iter()
                .enumerate()
                .map(|(i, p)| (p.clone(), i))
                .collect(),
            functions: definitions
                .iter()
                .enumerate()
                .take(generated_definition_start)
                .filter(|(_, entry)| entry.args.is_empty() && !entry.is_alias)
                .map(|(i, entry)| (entry.lhs.clone(), i))
                .collect(),
            families: definitions
                .iter()
                .enumerate()
                .filter_map(|(i, entry)| {
                    entry
                        .lhs
                        .as_fun_view()
                        .map(|call| ((call.get_symbol(), entry.tags.clone()), i))
                })
                .collect(),
            tag_counts,
            boolean: self
                .parameters
                .iter()
                .enumerate()
                .filter(|(i, _)| {
                    self.rows
                        .iter()
                        .all(|row| row[*i].is_zero() || row[*i] == Rational::one())
                })
                .map(|(_, p)| p.clone())
                .collect(),
            closed: HashMap::new(),
            kinds: HashMap::new(),
        };
        let mut candidates = BTreeMap::new();
        let all_rows = (0..self.rows.len()).collect::<Vec<_>>();
        for (owner, atom) in atoms.iter().enumerate() {
            if eligible_owners[owner] {
                analysis.collect(
                    atom.get_root().as_view(),
                    owner,
                    atom,
                    &all_rows,
                    &mut HashSet::new(),
                    &mut candidates,
                );
            }
        }
        let mut selected = candidates
            .into_iter()
            .filter(|(_, active)| active.len() > 1)
            .map(|((owner, call), active)| {
                let kind = analysis
                    .kind(call.as_view(), owner, &atoms[owner])
                    .expect("selected scalar boundary");
                (owner, call, active, kind)
            })
            .collect::<Vec<_>>();
        if selected.is_empty() {
            return None;
        }
        // Fixed outputs form a prefix copied once, followed by original distinct
        // row-dependent caller ports. Stable source ordering fixes both layouts.
        selected.sort_by_key(|(_, _, _, kind)| {
            !matches!(kind, ScalarAuxiliaryKind::Scalar { columns, .. } if columns.is_empty())
        });
        let ports = (0..selected.len())
            .map(|i| function!(port_symbol, i))
            .collect::<Vec<_>>();
        let mut originals = Vec::with_capacity(atoms.len());
        for (owner, atom) in atoms.iter_mut().enumerate() {
            let replacements = selected
                .iter()
                .zip(&ports)
                .filter(|((index, _, _, _), _)| *index == owner)
                .map(|((_, call, _, _), port)| (call.clone(), port.clone()))
                .collect::<HashMap<_, _>>();
            let original = std::mem::take(atom);
            *atom = EvaluatorStack::specialize_aliases(
                original.get_root(),
                &original,
                |body| {
                    body.replace_map(|view, _, out| {
                        if physical.contains(view.get_data()) {
                            out.set_from_view(&view);
                        } else if let Some(port) = replacements.get(view.get_data()) {
                            out.set_from_view(&port.as_view());
                        } else if let AtomView::Fun(call) = view {
                            // Collection and rewriting must have the same lazy
                            // boundary, even when a selected call also occurs
                            // inside an unexamined physical branch.
                            if call.get_symbol() != Symbol::IF
                                || call.get_nargs() != 3
                                || analysis
                                    .guard(call.iter().next().unwrap(), &original, None)
                                    .is_none()
                            {
                                out.set_from_view(&view);
                            }
                        }
                    })
                },
                None,
            );
            originals.push(original);
        }
        // Selected parents retain their original local DAG. As in physical-only
        // extraction, remove dependency-only ports after tracing the final owner.
        let port_indices = ports
            .iter()
            .enumerate()
            .map(|(i, p)| (p.clone(), i))
            .collect::<HashMap<_, _>>();
        let mut remaining = BTreeSet::new();
        for atom in atoms.iter() {
            for body in std::iter::once(atom.get_root()).chain(atom.get_aliases().values()) {
                body.visitor(&mut |view| {
                    if let Some(i) = port_indices.get(view.get_data()) {
                        remaining.insert(*i);
                        false
                    } else {
                        true
                    }
                });
            }
        }
        let mut kept_ports = Vec::new();
        let mut roots = Vec::new();
        let mut used = BTreeSet::new();
        let mut auxiliary_parameters = Vec::new();
        let mut auxiliary_values = Vec::new();
        let mut value_ports = BTreeMap::<(usize, Rational), Atom>::new();
        let mut output_indices = HashMap::<(Option<usize>, Atom, Vec<Rational>), usize>::new();
        let mut gathers = vec![Vec::new(); self.rows.len()];
        let mut shared = 0;
        for (index, (owner, call, active, kind)) in selected.into_iter().enumerate() {
            if !remaining.contains(&index) {
                continue;
            }
            kept_ports.push(ports[index].clone());
            if let ScalarAuxiliaryKind::Scalar {
                definitions,
                columns,
            } = &kind
                && columns.is_empty()
            {
                used.extend(definitions.iter().copied());
                roots.push(EvaluatorStack::specialize_aliases(
                    &call,
                    &originals[owner],
                    Clone::clone,
                    None,
                ));
                shared += 1;
                continue;
            }
            if let ScalarAuxiliaryKind::Numerator(_, _, dependencies)
            | ScalarAuxiliaryKind::Scalar {
                definitions: dependencies,
                ..
            } = &kind
            {
                used.extend(dependencies.iter().copied());
            }
            for (row, gather) in gathers.iter_mut().enumerate() {
                if !active.contains(&row) {
                    gather.push(None);
                    continue;
                }
                let (scope, root, values, replacements) = match &kind {
                    ScalarAuxiliaryKind::Selector | ScalarAuxiliaryKind::Scalar { .. } => {
                        let leaf = if matches!(kind, ScalarAuxiliaryKind::Selector) {
                            analysis.leaf(call.as_view(), owner, &originals[owner], row)
                        } else {
                            call.clone()
                        };
                        if leaf.is_zero() {
                            gather.push(None);
                            continue;
                        }
                        let mut selected_columns = BTreeSet::new();
                        // A scalar retains the original caller DAG and its
                        // certified inputs. Only selectors choose a new leaf.
                        let columns = if let ScalarAuxiliaryKind::Scalar { columns, .. } = &kind {
                            columns
                        } else if let Some(ScalarAuxiliaryKind::Scalar {
                            definitions,
                            columns,
                        }) = analysis
                            .kinds
                            .get(&(owner, leaf.clone()))
                            .cloned()
                            .flatten()
                        {
                            // Reuse only an existing certificate; leave the
                            // selector-alias traversal cache unchanged.
                            used.extend(definitions);
                            selected_columns = columns;
                            &selected_columns
                        } else {
                            // A Leaf can still require a zero-argument function
                            // definition; retain the complete analysis fallback.
                            let mut dependencies = BTreeSet::new();
                            analysis
                                .scalar_dependencies(
                                    leaf.as_view(),
                                    Some(&originals[owner]),
                                    &mut HashSet::new(),
                                    &mut HashSet::new(),
                                    &mut dependencies,
                                    &mut selected_columns,
                                )
                                .expect("resolved Boolean selector leaf");
                            used.extend(dependencies);
                            &selected_columns
                        };
                        let mut values = Vec::with_capacity(columns.len());
                        let mut replacements = HashMap::new();
                        // Each distinct original column keeps its own role,
                        // including equal-valued columns in divisions. The
                        // values remain native inputs, never literal atoms.
                        for (role, &column) in columns.iter().enumerate() {
                            let value = self.rows[row][column].clone();
                            values.push(value.clone());
                            let port =
                                value_ports.entry((role, value.clone())).or_insert_with(|| {
                                    let port =
                                        function!(port_symbol, role, auxiliary_parameters.len());
                                    auxiliary_parameters.push(port.clone());
                                    auxiliary_values.push(value);
                                    port
                                });
                            replacements.insert(self.parameters[column].clone(), port.clone());
                        }
                        (Some(owner), leaf, values, replacements)
                    }
                    ScalarAuxiliaryKind::Numerator(definition, columns, _) => {
                        let actuals = call
                            .as_fun_view()
                            .expect("numerator call")
                            .iter()
                            .skip(definitions[*definition].tags.len())
                            .map(|arg| arg.to_owned())
                            .collect::<Vec<_>>();
                        let mut first_positions = HashMap::new();
                        let mut replacements = HashMap::new();
                        for (position, (argument, column)) in
                            actuals.iter().zip(columns).enumerate()
                        {
                            if let Some(column) = column {
                                let role = *first_positions.entry(*column).or_insert(position);
                                let value = self.rows[row][*column].clone();
                                let port =
                                    value_ports.entry((role, value.clone())).or_insert_with(|| {
                                        let port = function!(
                                            port_symbol,
                                            role,
                                            auxiliary_parameters.len()
                                        );
                                        auxiliary_parameters.push(port.clone());
                                        auxiliary_values.push(value);
                                        port
                                    });
                                replacements.insert(argument.clone(), port.clone());
                            }
                        }
                        // Actual arguments alone are renamed. The role retains
                        // this call's equality partition even for equal values
                        // in distinct original coefficient columns.
                        let function = call.as_fun_view().expect("numerator call");
                        let tag_count = definitions[*definition].tags.len();
                        let arguments = function
                            .iter()
                            .enumerate()
                            .map(|(position, argument)| {
                                if position < tag_count {
                                    argument.to_owned()
                                } else {
                                    replacements
                                        .get(argument.get_data())
                                        .cloned()
                                        .unwrap_or_else(|| argument.to_owned())
                                }
                            })
                            .collect::<Vec<_>>();
                        (
                            None,
                            function.get_symbol().call_args(arguments),
                            Vec::new(),
                            HashMap::new(),
                        )
                    }
                    _ => unreachable!("dynamic auxiliary kind"),
                };
                let next = roots.len();
                let output = *output_indices
                    // The alias handle alone does not identify a scalar
                    // output: its owned RHS also depends on this exact tuple.
                    .entry((scope, root.clone(), values))
                    .or_insert_with(|| {
                        roots.push(match scope {
                            Some(owner) => EvaluatorStack::specialize_aliases(
                                &root,
                                &originals[owner],
                                |body| {
                                    body.replace_map(|view, _, out| {
                                        if physical.contains(view.get_data())
                                            || originals[owner]
                                                .get_aliases()
                                                .contains_key(view.get_data())
                                        {
                                            out.set_from_view(&view);
                                        } else if let Some(port) = replacements.get(view.get_data())
                                        {
                                            out.set_from_view(&port.as_view());
                                        } else if analysis.functions.contains_key(view.get_data()) {
                                            out.set_from_view(&view);
                                        }
                                    })
                                },
                                None,
                            ),
                            None => AliasedAtom::from(root),
                        });
                        next
                    });
                gather.push(Some(output));
            }
        }
        let row_lookup = self
            .parameters
            .iter()
            .enumerate()
            .map(|(i, p)| (p.clone(), i))
            .collect::<HashMap<_, _>>();
        let nonrow = physical_parameters
            .iter()
            .chain(&kept_ports)
            .cloned()
            .collect::<HashSet<_>>();
        let mut opaque_function = false;
        let mut columns = BTreeSet::new();
        for atom in atoms.iter() {
            for body in std::iter::once(atom.get_root()).chain(atom.get_aliases().values()) {
                body.visitor(&mut |view| {
                    if nonrow.contains(view.get_data()) {
                        false
                    } else if let Some(column) = row_lookup.get(view.get_data()) {
                        columns.insert(*column);
                        false
                    } else if atom.get_aliases().contains_key(view.get_data()) {
                        false
                    } else {
                        if let AtomView::Fun(call) = view {
                            let known = analysis.function(view);
                            // Tags select the definition; they are not native
                            // arguments or captured globals. Its RHS is scanned
                            // separately below, including free row parameters.
                            if known.as_ref().is_some_and(|(_, args)| args.is_empty()) {
                                return false;
                            }
                            opaque_function |= !call.get_symbol().is_builtin() && known.is_none();
                        }
                        true
                    }
                });
            }
        }
        // Unsplit ordinary functions can capture global row ports that are
        // absent from their caller arguments. Keep those inputs independently
        // of the roots' reachable-column scan; formal bindings shadow globals.
        for definition in definitions {
            let formals = definition
                .args
                .iter()
                .cloned()
                .map(Atom::from)
                .collect::<HashSet<_>>();
            definition.rhs.visitor(&mut |view| {
                if formals.contains(view.get_data()) || nonrow.contains(view.get_data()) {
                    false
                } else if let Some(column) = row_lookup.get(view.get_data()) {
                    columns.insert(*column);
                    false
                } else {
                    if let AtomView::Fun(call) = view {
                        let known = analysis.function(view);
                        if known.as_ref().is_some_and(|(_, args)| args.is_empty()) {
                            return false;
                        }
                        opaque_function |= !call.get_symbol().is_builtin() && known.is_none();
                    }
                    true
                }
            });
        }
        // A directly registered overload can share a head with a recorded
        // definition while capturing different globals. Only an exact known
        // head/tags/arity signature permits pruning the original input list.
        drop(analysis);
        self.auxiliary_parameters = auxiliary_parameters;
        self.auxiliary_values = auxiliary_values;
        self.auxiliary_shared_outputs = shared;
        self.auxiliary_gathers = gathers;
        self.input_columns = if opaque_function {
            (0..self.parameters.len()).collect()
        } else {
            columns.into_iter().collect()
        };
        (!kept_ports.is_empty()).then(|| {
            (
                kept_ports,
                roots,
                used.into_iter().map(|i| definitions[i].clone()).collect(),
            )
        })
    }

    pub(super) fn evaluate<T: FloatLike>(
        &self,
        evaluator: &mut GenericEvaluator,
        physical_input: &[Complex<F<T>>],
        auxiliary_output: &[Complex<F<T>>],
        metadata: &mut EvaluationMetaData,
    ) -> Vec<DualOrNot<Complex<F<T>>>> {
        let shape_len = evaluator.dual_shape.as_ref().map_or(1, Vec::len);
        let prototype = physical_input
            .first()
            .map_or_else(|| F(T::new_zero()), |value| value.re.clone());
        let zero = Complex::new_re(prototype.zero());
        let mut input = physical_input.to_vec();
        let dynamic_ports = self.auxiliary_gathers.first().map_or(0, Vec::len);
        let dynamic_start = physical_input.len() - dynamic_ports * shape_len;
        // Inputs are read-only to the native evaluator. Keep the previous row's
        // gathered jets and copy only changed output identities; an inactive
        // source resets its full jet to zero before the next row executes.
        input[dynamic_start..].fill(zero.clone());
        let mut previous = vec![None; dynamic_ports];
        input.resize(
            input.len() + self.input_columns.len() * shape_len,
            zero.clone(),
        );
        let mut result: Option<Vec<DualOrNot<Complex<F<T>>>>> = None;
        let start = std::time::Instant::now();
        let mut evaluate = <T as GenericEvaluatorFloat>::get_evaluator(evaluator);
        metadata.evaluator_evaluation_time = metadata
            .evaluator_evaluation_time
            .saturating_add(start.elapsed());
        for (row_index, row) in self.rows.iter().enumerate() {
            if dynamic_ports != 0 {
                for (port, &output) in self.auxiliary_gathers[row_index].iter().enumerate() {
                    if output == previous[port] {
                        continue;
                    }
                    let target = &mut input
                        [dynamic_start + port * shape_len..dynamic_start + (port + 1) * shape_len];
                    if let Some(output) = output {
                        target.clone_from_slice(
                            &auxiliary_output[output * shape_len..(output + 1) * shape_len],
                        );
                    } else {
                        target.fill(zero.clone());
                    }
                    previous[port] = output;
                }
            }
            for (slot, &column) in self.input_columns.iter().enumerate() {
                input[physical_input.len() + slot * shape_len] =
                    Complex::new_re(prototype.from_rational(&row[column]));
            }
            let start = std::time::Instant::now();
            let output = evaluate(&input);
            metadata.evaluator_evaluation_time = metadata
                .evaluator_evaluation_time
                .saturating_add(start.elapsed());
            if let Some(result) = &mut result {
                for (sum, value) in result.iter_mut().zip(output) {
                    *sum += value;
                }
            } else {
                result = Some(output);
            }
        }
        result.expect("parametric residue table is nonempty")
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use symbolica::{evaluate::InliningPolicy, prelude::Indeterminate};

    use crate::{
        integrands::process::param_builder::FnMapEntry, utils::GS, uv::ParametricIntegrandRow,
    };

    use super::*;

    #[test]
    fn physical_function_extraction_respects_aliases_rows_and_lexical_boundaries() {
        crate::initialisation::test_initialise().unwrap();
        let x = function!(symbol!("parametric_row_test::physical_input"), 0);
        let q = Atom::var(symbol!("parametric_row_test::physical_row"));
        let head = symbol!("parametric_row_test::physical_function");
        let energy = function!(head, 0);
        let captured = function!(head, 1);
        let pole = function!(head, 2);
        let once = function!(head, 3);
        let cyclic = function!(head, 4);
        let unknown = function!(head, 5);
        let caller_alias = function!(head, 6);
        let formal_wrapper = function!(head, 7, q.clone());
        let via_formal = function!(head, 8);
        let opaque = function!(
            symbol!("parametric_row_test::opaque_numerator"),
            q.clone(),
            pole.clone()
        );
        let mut definitions = [
            (x.clone(), x.clone()), // Exact input lookup must precede self definitions.
            (energy.clone(), (x.pow(2) + 1).pow((1, 2))),
            (captured.clone(), &x + &q),
            (pole.clone(), x.pow(-1)),
            (once.clone(), (&x + 2).pow(-1)),
            (cyclic.clone(), cyclic.clone()),
            (
                unknown.clone(),
                function!(symbol!("parametric_row_test::unrecorded"), x.clone()),
            ),
            (caller_alias.clone(), x.clone()),
            (formal_wrapper.clone(), q.clone()),
            (via_formal.clone(), function!(head, 7, x.clone())),
        ]
        .into_iter()
        .map(|(lhs, rhs)| FnMapEntry {
            tags: lhs
                .as_fun_view()
                .unwrap()
                .iter()
                .map(|a| a.to_owned())
                .collect(),
            lhs,
            rhs,
            args: Vec::new(),
            inlining: InliningPolicy::Always,
            is_alias: false,
        })
        .collect::<Vec<_>>();
        definitions[7].is_alias = true;
        definitions[8].tags = vec![Atom::num(7)];
        definitions[8].args = vec![Indeterminate::try_from(q.clone()).unwrap()];
        let alias_head = symbol!("parametric_row_test::physical_alias");
        let outer = function!(alias_head, 0);
        let leaf = function!(alias_head, 1);
        let dead = function!(alias_head, 2);
        let mut source = AliasedAtom::from(
            Symbol::IF.call_args([q.clone(), outer.clone(), once.clone()])
                + Symbol::IF.call_args([x.clone(), pole.clone(), Atom::Zero])
                + Symbol::IF.call_args([&x + &q, pole.clone(), Atom::Zero])
                + &captured
                + &cyclic
                + &unknown
                + &caller_alias
                + &via_formal
                + &opaque,
        );
        source.register_alias(outer.clone(), &leaf + 1);
        source.register_alias(leaf.clone(), energy.clone());
        source.register_alias(dead, pole.clone());
        let mut rows = ParametricResidueRows {
            parameters: vec![q.clone()],
            input_columns: vec![0],
            rows: [0, 1, 1].map(|value| vec![Rational::from(value)]).to_vec(),
            ..Default::default()
        };
        let original = source.clone();
        let port_symbol = symbol!("parametric_row_test::physical_port");
        let (ports, roots, retained) = rows
            .extract_physical_scalars(
                std::slice::from_mut(&mut source),
                std::slice::from_ref(&x),
                &definitions,
                definitions.len(),
                port_symbol,
            )
            .unwrap();
        assert_eq!(roots.len(), 3);
        // The ordinary wrapper's formal q shadows the global row input. Its
        // unchanged call and both definitions now share one physical output.
        assert_eq!(roots[0].get_root(), &via_formal);
        assert!(roots[0].get_aliases().is_empty());
        assert_eq!(roots[1].get_root(), &once);
        assert!(roots[1].get_aliases().is_empty());
        assert_eq!(roots[2].get_root(), &outer);
        assert_eq!(roots[2].get_aliases()[&leaf], energy);
        assert_eq!(roots[2].get_aliases().len(), 2);
        assert_eq!(ports.len(), 2);
        assert_eq!(
            retained,
            [1, 4, 8, 9].map(|index| definitions[index].clone())
        );
        assert_eq!(rows.auxiliary_shared_outputs, 1);
        assert_eq!(
            rows.auxiliary_gathers,
            vec![vec![Some(1)], vec![Some(2)], vec![Some(2)]]
        );
        assert_eq!(
            rows.input_columns,
            vec![0],
            "the unsplit captured function still receives its global row input"
        );
        assert!(source.get_aliases().is_empty());
        assert!(source.get_root().contains(&ports[0]));
        assert!(
            source.get_root().contains(&opaque),
            "opaque actual arguments stay intact"
        );
        assert!(source.get_root().contains(&captured));
        assert!(
            source.get_root().contains(&pole),
            "physical branches retain their pole"
        );
        let mut repeated = original;
        let again = rows
            .extract_physical_scalars(
                std::slice::from_mut(&mut repeated),
                &[x],
                &definitions,
                definitions.len(),
                port_symbol,
            )
            .unwrap();
        assert_eq!(again.0, ports);
        assert_eq!(again.2, retained);
        assert_eq!(again.1.len(), roots.len());
        for (left, right) in again.1.iter().zip(&roots) {
            assert_eq!(left.get_root(), right.get_root());
            assert_eq!(left.get_aliases(), right.get_aliases());
        }
        assert_eq!(source.get_root(), repeated.get_root());
        assert_eq!(source.get_aliases(), repeated.get_aliases());
    }

    #[test]
    fn physical_alias_outputs_keep_only_reachable_owner_dags() {
        crate::initialisation::test_initialise().unwrap();
        let x = Atom::var(symbol!("parametric_row_test::large_dag_x"));
        let y = Atom::var(symbol!("parametric_row_test::large_dag_y"));
        let q = Atom::var(symbol!("parametric_row_test::large_dag_row"));
        let leaf = Atom::var(symbol!("parametric_row_test::large_dag_leaf"));
        let head = symbol!("parametric_row_test::large_dag_alias");
        let parent = function!(head, 0);
        let pole = function!(head, 1);
        let cycle_a = function!(head, 2);
        let cycle_b = function!(head, 3);
        let spokes = (0..128)
            .map(|index| function!(head, index + 4))
            .collect::<Vec<_>>();
        let mut sources = [&x, &y]
            .into_iter()
            .enumerate()
            .map(|(owner, input)| {
                let mut source = AliasedAtom::from(
                    (&q + 1) * &parent
                        + Symbol::IF.call_args([x.clone(), pole.clone(), Atom::Zero]),
                );
                // Identical handles in different owners deliberately have different
                // bodies. The non-function leaf also exercises unscoped selection.
                source.register_alias(leaf.clone(), input.pow(2) + Atom::num(owner + 1));
                source.register_alias(parent.clone(), Atom::add_many(&spokes));
                for (index, spoke) in spokes.iter().enumerate() {
                    source.register_alias(spoke.clone(), &leaf + Atom::num(index));
                }
                source.register_alias(pole.clone(), x.pow(-1));
                // Neither disconnected cycles nor unused singular definitions belong
                // to an auxiliary output's reachable DAG.
                source.register_alias(cycle_a.clone(), &cycle_b + 1);
                source.register_alias(cycle_b.clone(), &cycle_a + 1);
                for index in 0..512 {
                    source
                        .register_alias(function!(head, index + 132), x.pow(-1) + Atom::num(index));
                }
                source
            })
            .collect::<Vec<_>>();
        let originals = sources.clone();
        let mut rows = ParametricResidueRows {
            parameters: vec![q.clone()],
            input_columns: vec![0],
            rows: vec![vec![Rational::zero()], vec![Rational::one()]],
            ..Default::default()
        };
        let (ports, roots, definitions) = rows
            .extract_physical_scalars(
                &mut sources,
                &[x.clone(), y.clone()],
                &[],
                0,
                symbol!("parametric_row_test::large_dag_port"),
            )
            .unwrap();
        assert_eq!(ports.len(), 2, "dependency-only outputs are omitted");
        assert_eq!(roots.len(), 4);
        assert!(definitions.is_empty());
        assert_eq!(rows.auxiliary_shared_outputs, 0);
        assert!(rows.input_columns.is_empty());
        for (owner, (source, port)) in sources.iter().zip(&ports).enumerate() {
            assert_eq!(source.get_aliases().len(), 1);
            assert_eq!(source.get_aliases()[&pole], x.pow(-1));
            assert_eq!(
                source.get_root(),
                &(port + Symbol::IF.call_args([x.clone(), pole.clone(), Atom::Zero]))
            );
            for row in 0..2 {
                let root = &roots[rows.auxiliary_gathers[row][owner].unwrap()];
                let parameter = rows
                    .auxiliary_parameters
                    .iter()
                    .zip(&rows.auxiliary_values)
                    .find(|(_, value)| **value == row)
                    .unwrap()
                    .0;
                assert_eq!(root.get_root(), &((parameter + 1) * &parent));
                assert_eq!(root.get_aliases().len(), spokes.len() + 2);
                for (key, body) in root.get_aliases() {
                    assert_eq!(body, &originals[owner].get_aliases()[key]);
                }
                assert!(!root.get_aliases().contains_key(&pole));
                assert!(!root.get_aliases().contains_key(&cycle_a));
                assert!(!root.get_aliases().contains_key(&cycle_b));
                let input = if owner == 0 { &x } else { &y };
                assert_eq!(
                    root.clone().into_inner(),
                    (parameter + 1)
                        * (128 * (input.pow(2) + Atom::num(owner + 1)) + Atom::num(8128))
                );
            }
        }
        assert_ne!(
            roots[rows.auxiliary_gathers[0][0].unwrap()].get_aliases()[&leaf],
            roots[rows.auxiliary_gathers[0][1].unwrap()].get_aliases()[&leaf]
        );
        assert_ne!(ports[0], ports[1]);
    }

    #[test]
    fn selector_factors_preserve_sparse_rows_and_branch_local_poles() {
        crate::initialisation::test_initialise().unwrap();
        let x = Atom::var(symbol!("parametric_row_test::selector_factor_x"));
        let low = Atom::var(symbol!("parametric_row_test::selector_factor_low"));
        let high = Atom::var(symbol!("parametric_row_test::selector_factor_high"));
        let common = x.pow(-1) * (&x + 2).pow(-1);
        let selected = ParametricResidueRows::select(
            &[(0, &common * 2), (1, &common * 3)],
            &[low.clone(), high.clone()],
        );
        assert_eq!(
            selected,
            Symbol::IF.call_args([
                high.clone(),
                Atom::Zero,
                Symbol::IF.call_args([low.clone(), &common * 3, &common * 2]),
            ])
        );
        for (ordinal, expected) in [2, 3, 0, 0].into_iter().enumerate() {
            let concrete = selected.replace_multiple([
                Replacement::new(low.to_pattern(), Atom::num(ordinal & 1)),
                Replacement::new(high.to_pattern(), Atom::num(ordinal >> 1)),
            ]);
            assert_eq!(concrete, &common * expected);
            if ordinal >= 2 {
                assert_eq!(
                    concrete.replace(x.to_pattern()).with(Atom::Zero),
                    Atom::Zero
                );
            }
        }
        // A pole absent from another active branch cannot be factored through
        // that branch, even though rational collection could cancel it there.
        let local_pole = ParametricResidueRows::select(
            &[(0, x.pow(-1)), (1, Atom::one())],
            std::slice::from_ref(&low),
        );
        assert_eq!(
            local_pole,
            Symbol::IF.call_args([low.clone(), Atom::one(), x.pow(-1)])
        );

        // A deeper tree gathers complete original carriers, including an
        // ordinary denominator wrapper. Repeated leaves share outputs; absent
        // rows never evaluate their pole, even through three selector levels.
        let middle = Atom::var(symbol!("parametric_row_test::selector_factor_middle"));
        let bits = vec![low, middle, high];
        let argument = Atom::var(symbol!("parametric_row_test::selector_factor_argument"));
        let denominator = function!(GS.tree_denom_wrapper, common.clone());
        let carriers = [2, 3, 2, 4, 3]
            .into_iter()
            .enumerate()
            .map(|(row, weight)| (row, &denominator * weight))
            .collect::<Vec<_>>();
        let selected = ParametricResidueRows::select(&carriers, &bits);
        selected.visitor(&mut |view| {
            if let AtomView::Fun(call) = view
                && call.get_symbol() == Symbol::IF
            {
                for branch in call.iter().skip(1) {
                    assert!(
                        branch.is_zero()
                            || matches!(branch, AtomView::Fun(nested) if nested.get_symbol() == Symbol::IF)
                            || carriers.iter().any(|(_, carrier)| carrier.as_view() == branch),
                        "selector arms retain complete original carrier leaves"
                    );
                }
            }
            true
        });
        let definition = FnMapEntry {
            lhs: function!(GS.tree_denom_wrapper, argument.clone()),
            rhs: argument.clone(),
            args: vec![Indeterminate::try_from(argument).unwrap()],
            tags: Vec::new(),
            inlining: InliningPolicy::Always,
            is_alias: false,
        };
        let mut rows = ParametricResidueRows {
            parameters: bits,
            rows: (0..8)
                .map(|row| (0..3).map(|bit| Rational::from((row >> bit) & 1)).collect())
                .collect(),
            ..Default::default()
        };
        let mut sources = vec![AliasedAtom::from(selected)];
        let (ports, roots, definitions) = rows
            .extract_physical_scalars(
                &mut sources,
                std::slice::from_ref(&x),
                std::slice::from_ref(&definition),
                1,
                symbol!("parametric_row_test::selector_factor_output"),
            )
            .unwrap();
        assert_eq!(ports.len(), 1);
        assert_eq!(roots.len(), 3, "one output per unique original carrier");
        assert_eq!(definitions, vec![definition]);
        assert_eq!(sources[0].get_root(), &ports[0]);
        assert!(sources[0].get_aliases().is_empty());
        assert!(rows.input_columns.is_empty());
        assert_eq!(rows.auxiliary_shared_outputs, 0);
        assert!(rows.auxiliary_parameters.is_empty());
        for (row, gather) in rows.auxiliary_gathers.iter().enumerate() {
            assert_eq!(gather.len(), 1);
            if let Some((_, carrier)) = carriers.get(row) {
                let root = &roots[gather[0].unwrap()];
                assert_eq!(root.get_root(), carrier);
                assert!(root.get_aliases().is_empty());
                assert!(!root.get_root().contains_symbol(Symbol::IF));
            } else {
                assert_eq!(gather[0], None, "inactive pole has no output binding");
            }
        }
        assert_eq!(rows.auxiliary_gathers[0], rows.auxiliary_gathers[2]);
        assert_eq!(rows.auxiliary_gathers[1], rows.auxiliary_gathers[4]);
    }

    #[test]
    fn coefficient_rows_keep_exact_inputs_and_drop_only_cancelled_carriers() {
        crate::initialisation::test_initialise().unwrap();
        let q = Atom::var(symbol!("parametric_row_test::exact_coefficient"));
        let family = symbol!("parametric_row_test::exact_family");
        let precise = Rational::new(9_007_199_254_740_993i64, 9_007_199_254_740_992i64);
        let mut terms = [
            (q.clone(), vec![(Rational::one(), 1), (precise.clone(), 2)]),
            (q.pow(2), vec![(Rational::one(), 4)]),
            (
                q.pow(-1),
                vec![(Rational::zero(), 7), (Rational::zero(), -7)],
            ),
        ]
        .into_iter()
        .enumerate()
        .map(|(tag, (rhs, rows))| ParametricIntegrandTerm {
            numerator: Arc::new(FnMapEntry {
                lhs: function!(family, tag, q.clone()),
                rhs,
                tags: vec![Atom::num(tag)],
                args: vec![Indeterminate::try_from(q.clone()).unwrap()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            }),
            rows: rows
                .into_iter()
                .map(|(coefficient, carrier)| ParametricIntegrandRow {
                    coefficients: vec![coefficient],
                    carrier: Atom::num(carrier),
                    source_keys: Vec::new(),
                })
                .collect(),
        })
        .collect::<Vec<_>>();
        let (atoms, guards, rows) = ParametricResidueRows::lower(&Atom::Zero, &terms).unwrap();
        assert_eq!(atoms.len(), 2);
        assert_eq!(rows.rows.len(), 2, "binary64-equal rationals stay distinct");
        assert_eq!(rows.rows[0][0], Rational::one());
        assert_eq!(rows.rows[1][0], precise);
        for atom in &atoms {
            let mut calls = 0;
            atom.visitor(&mut |view| {
                if let Some(call) = view.as_fun_view()
                    && call.get_symbol() == family
                {
                    calls += 1;
                }
                true
            });
            assert_eq!(calls, 1, "one generic numerator call per family");
        }
        let joined = rows.join_preprocessed(atoms.into_iter().map(Into::into).collect(), &guards);
        for (index, expected) in [Atom::num(5), Atom::num(precise) * 2]
            .into_iter()
            .enumerate()
        {
            assert_eq!(
                rows.specialize(&joined, &rows.rows[index], index)
                    .into_inner()
                    .replace_multiple(terms.iter().map(|term| term.numerator.replacement())),
                expected
            );
        }
        let (atoms, guards, empty) =
            ParametricResidueRows::lower(&Atom::Zero, &terms[2..]).unwrap();
        assert_eq!(atoms, vec![Atom::Zero]);
        assert_eq!(guards, vec![Atom::one()]);
        assert!(empty.parameters.is_empty());
        assert_eq!(empty.rows, vec![Vec::<Rational>::new()]);

        terms.reverse();
        let (_, _, reversed) = ParametricResidueRows::lower(&Atom::Zero, &terms).unwrap();
        assert_eq!(
            rows, reversed,
            "catalogue order depends only on exact inputs"
        );
        terms[0].rows[0].coefficients.clear();
        assert!(ParametricResidueRows::lower(&Atom::Zero, &terms).is_err());
    }

    #[test]
    fn ordinal_rows_share_exact_columns_without_changing_formal_bindings() {
        crate::initialisation::test_initialise().unwrap();
        let q = Atom::var(symbol!("parametric_row_test::ordinal_q"));
        let r = Atom::var(symbol!("parametric_row_test::ordinal_r"));
        let family = symbol!("parametric_row_test::ordinal_family");
        let mut terms = [
            (
                vec![q.clone(), r.clone()],
                &q + &r * 2,
                vec![vec![1, 1], vec![3, 3], vec![7, 7]],
                1,
            ),
            // The formal order differs, and the two columns remain distinct.
            (
                vec![r.clone(), q.clone()],
                &r + &q * 3,
                vec![vec![2, 5], vec![5, 2]],
                1,
            ),
            (vec![q.clone()], &q * 5, vec![vec![1], vec![3], vec![7]], 2),
            // A singleton source is constant, with no numeric coefficient port.
            (vec![q.clone()], q, vec![vec![11]], 3),
        ]
        .into_iter()
        .enumerate()
        .map(
            |(tag, (args, rhs, rows, carrier))| ParametricIntegrandTerm {
                numerator: Arc::new(FnMapEntry {
                    lhs: family
                        .call_args(std::iter::once(Atom::num(tag)).chain(args.iter().cloned())),
                    rhs,
                    tags: vec![Atom::num(tag)],
                    args: args
                        .into_iter()
                        .map(|arg| Indeterminate::try_from(arg).unwrap())
                        .collect(),
                    inlining: InliningPolicy::Always,
                    is_alias: false,
                }),
                rows: rows
                    .into_iter()
                    .map(|coefficients| ParametricIntegrandRow {
                        coefficients: coefficients.into_iter().map(Rational::from).collect(),
                        carrier: Atom::num(carrier),
                        source_keys: Vec::new(),
                    })
                    .collect(),
            },
        )
        .collect::<Vec<_>>();
        let (atoms, guards, rows) = ParametricResidueRows::lower(&Atom::num(13), &terms).unwrap();
        assert_eq!(rows.rows.len(), 3);
        assert_eq!(
            rows.parameters.len(),
            5,
            "three distinct varying columns and two selector bits"
        );
        assert_eq!(guards[0], Atom::one());
        assert_eq!(guards[2], Atom::one());
        // Retain the independently preprocessed scalar bodies as native aliases.
        // Joining equal guards must neither expand their numerator calls nor
        // lose the alias ownership used by later coefficient specialization.
        let alias = symbol!("parametric_row_test::ordinal_join_alias");
        let aliased = atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| {
                let handle = function!(alias, index);
                let mut retained = AliasedAtom::from(handle.clone());
                retained.register_alias(handle, atom.clone());
                retained
            })
            .collect();
        let joined = rows.join_preprocessed(aliased, &guards);
        assert_eq!(joined.get_aliases().len(), atoms.len());
        for (index, atom) in atoms.iter().enumerate() {
            assert_eq!(joined.get_aliases()[&function!(alias, index)], *atom);
        }
        let expected_root = function!(alias, 0)
            + function!(alias, 2)
            + Symbol::IF.call_args([guards[1].clone(), function!(alias, 1), Atom::Zero])
            + Symbol::IF.call_args([
                guards[3].clone(),
                function!(alias, 3) + function!(alias, 4),
                Atom::Zero,
            ]);
        assert_eq!(joined.get_root(), &expected_root);
        for (index, expected) in [76, 50, 91].into_iter().enumerate() {
            assert_eq!(
                rows.specialize(&joined, &rows.rows[index], index)
                    .into_inner()
                    .replace_multiple(terms.iter().map(|term| term.numerator.replacement())),
                Atom::num(expected)
            );
        }
        terms.reverse();
        let (_, _, reversed) = ParametricResidueRows::lower(&Atom::num(13), &terms).unwrap();
        assert_eq!(
            rows, reversed,
            "exact column numbering is independent of family order"
        );
    }

    #[test]
    fn row_specialization_visits_only_reachable_aliases() {
        crate::initialisation::test_initialise().unwrap();
        let coefficient = Atom::var(symbol!("parametric_row_test::coefficient"));
        let selector = Atom::var(symbol!("parametric_row_test::selector"));
        let x = Atom::var(symbol!("parametric_row_test::x"));
        let head = symbol!("parametric_row_test::alias");
        let outer = function!(head, 0);
        let leaf = function!(head, 1);
        let dead = function!(head, 2);
        let mut source = AliasedAtom::from(Symbol::IF.call_args([
            selector.clone(),
            dead.clone(),
            outer.clone(),
        ]));
        source.register_alias(
            outer,
            Symbol::IF.call_args([selector.clone(), dead.clone(), leaf.clone()]) + &leaf,
        );
        source.register_alias(leaf, &coefficient * &x);
        source.register_alias(dead, selector.pow(-1));
        let catalog = ParametricResidueRows {
            input_columns: vec![0, 1],
            parameters: vec![coefficient, selector],
            rows: vec![
                vec![Rational::from(2), Rational::zero()],
                vec![Rational::from(7), Rational::from(1)],
            ],
            ..Default::default()
        };
        let selected = catalog.specialize(&source, &catalog.rows[0], 0);
        assert_eq!(selected.get_aliases().len(), 2);
        assert!(!selected.get_aliases().contains_key(&function!(head, 2, 0)));
        assert_eq!(selected.into_inner(), 4 * x);
        let selected = catalog.specialize(&source, &catalog.rows[1], 1);
        assert_eq!(selected.get_aliases().len(), 1);
        assert_eq!(selected.into_inner(), Atom::num(1));
        assert_eq!(source.get_aliases().len(), 3);
    }
}
