use std::collections::BTreeMap;

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use symbolica::prelude::{
    AliasedAtom, Atom, AtomCore, Rational, Replacement, Symbol, function, symbol,
};

use crate::GammaLoopContext;
use crate::uv::ParametricIntegrandTerm;

use super::EvaluatorStack;
#[cfg(test)]
use symbolica::prelude::AtomView;

/// Exact coefficients and selectors assembled from independent physical source sums.
/// Symbolica owns their evaluator specialization and execution.
#[derive(Clone, Debug, Default, PartialEq, Eq, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct ParametricResidueRows {
    pub(crate) parameters: Vec<Atom>,
    pub(crate) rows: Vec<Vec<Rational>>,
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
                let mut grouped = BTreeMap::<Vec<Rational>, Vec<&Atom>>::new();
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
            .collect::<BTreeMap<_, _>>();
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
        Ok((atoms, guards, Self { parameters, rows }))
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
        let mut groups = BTreeMap::<&Atom, Vec<&Atom>>::new();
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
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::{iter, slice};

    use symbolica::evaluate::{FunctionMap, InliningPolicy};
    use symbolica::prelude::Indeterminate;

    use crate::initialisation;
    use crate::integrands::process::param_builder::FnMapEntry;
    use crate::utils::GS;
    use crate::uv::ParametricIntegrandRow;

    use super::*;

    #[test]
    fn native_row_plan_respects_aliases_rows_and_lexical_boundaries() {
        initialisation::test_initialise().unwrap();
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
                + &caller_alias
                + &via_formal,
        );
        source.register_alias(outer.clone(), &leaf + 1);
        source.register_alias(leaf.clone(), energy.clone());
        source.register_alias(dead, pole.clone());
        let rows = ParametricResidueRows {
            parameters: vec![q.clone()],
            rows: [0, 1, 1].map(|value| vec![Rational::from(value)]).to_vec(),
        };
        let mut functions = FunctionMap::new();
        for definition in &definitions {
            definition.register(&mut functions).unwrap();
        }
        let inputs = [x, q];
        let build = || {
            source
                .get_root()
                .evaluator(&inputs)
                .function_map(functions.clone())
                .add_aliases(
                    source
                        .get_aliases()
                        .iter()
                        .map(|(key, body)| (key.clone(), body.clone())),
                )
                .unwrap()
        };
        let mut staged = build()
            .build_parameter_rows(&rows.parameters, &rows.rows)
            .unwrap()
            .map_coeff(&|value| value.re.to_f64());
        let mut summed = build()
            .sum_over_parameters(&rows.parameters, &rows.rows)
            .unwrap()
            .build()
            .unwrap()
            .map_coeff(&|value| value.re.to_f64());
        for point in [2.0_f64, 3., -3.] {
            let expected = 1. / (point + 2.)
                + 2. * ((point * point + 1.).sqrt() + 1.)
                + 6. / point
                + 9. * point
                + 2.;
            let mut output = [0.];
            staged
                .evaluate(&[point], &mut output, Rational::to_f64)
                .unwrap();
            assert!((output[0] - expected).abs() < 1e-12);
            assert!((summed.evaluate_single(&[point]) - expected).abs() < 1e-12);
        }
        // These were opaque to the old extraction pass. The complete native
        // builder owns their validation now; it must reject each unsupported
        // root instead of classifying it as a reusable scalar.
        for invalid in [cyclic, unknown, opaque] {
            assert!(
                invalid
                    .evaluator(&inputs)
                    .function_map(functions.clone())
                    .build_parameter_rows(&rows.parameters, &rows.rows)
                    .is_err()
            );
        }
    }

    #[test]
    fn physical_alias_outputs_keep_only_reachable_owner_dags() {
        initialisation::test_initialise().unwrap();
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
        let sources = [&x, &y]
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
        let sources = sources
            .into_iter()
            .map(AliasedAtom::scope_aliases)
            .collect::<Vec<_>>();
        let rows = ParametricResidueRows {
            parameters: vec![q.clone()],
            rows: vec![vec![Rational::zero()], vec![Rational::one()]],
        };
        let inputs = [x, y, q];
        let roots = sources
            .iter()
            .map(AliasedAtom::get_root)
            .collect::<Vec<_>>();
        let mut functions = FunctionMap::new();
        functions
            .add_aliases(sources.iter().flat_map(|source| {
                source
                    .get_aliases()
                    .iter()
                    .map(|(key, body)| (key.clone(), body.clone()))
            }))
            .unwrap();
        let original = Atom::evaluator_multiple(&roots, &inputs)
            .function_map(functions.clone())
            .build()
            .unwrap();
        let reachable = functions.retain_reachable(&roots, &inputs);
        let pruned = Atom::evaluator_multiple(&roots, &inputs)
            .function_map(reachable.clone())
            .build()
            .unwrap();
        assert_eq!(
            original.count_operations(),
            pruned.count_operations(),
            "disconnected cycles and unused poles never enter the native program"
        );
        let mut staged = Atom::evaluator_multiple(&roots, &inputs)
            .function_map(reachable)
            .build_parameter_rows(&rows.parameters, &rows.rows)
            .unwrap()
            .map_coeff(&|value| value.re.to_f64());
        for point in [[2.0_f64, 3.], [0., 1.]] {
            let mut output = [0.; 2];
            staged
                .evaluate(&point, &mut output, Rational::to_f64)
                .unwrap();
            for (owner, input) in point.iter().enumerate() {
                let pole = if point[0] == 0. { 0. } else { 2. / point[0] };
                let expected = 3. * (128. * (input * input + owner as f64 + 1.) + 8128.) + pole;
                assert_eq!(output[owner], expected);
            }
            assert_ne!(
                output[0], output[1],
                "identically named private aliases keep different owners"
            );
        }
    }

    #[test]
    fn selector_factors_preserve_sparse_rows_and_branch_local_poles() {
        initialisation::test_initialise().unwrap();
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
            slice::from_ref(&low),
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
        let rows = ParametricResidueRows {
            parameters: bits,
            rows: (0..8)
                .map(|row| (0..3).map(|bit| Rational::from((row >> bit) & 1)).collect())
                .collect(),
        };
        let inputs = iter::once(x)
            .chain(rows.parameters.iter().cloned())
            .collect::<Vec<_>>();
        let mut functions = FunctionMap::new();
        definition.register(&mut functions).unwrap();
        let build = || selected.evaluator(&inputs).function_map(functions.clone());
        for (row, weight) in [2, 3, 2, 4, 3, 0, 0, 0].into_iter().enumerate() {
            let mut evaluator = build()
                .sum_over_parameters(&rows.parameters, &rows.rows[row..row + 1])
                .unwrap()
                .build()
                .unwrap()
                .map_coeff(&|value| value.re.to_f64());
            assert_eq!(evaluator.evaluate_single(&[2.]), weight as f64 / 8.);
            if weight == 0 {
                assert_eq!(
                    evaluator.evaluate_single(&[0.]),
                    0.,
                    "inactive pole is never evaluated"
                );
            }
        }
        let mut staged = build()
            .build_parameter_rows(&rows.parameters, &rows.rows)
            .unwrap()
            .map_coeff(&|value| value.re.to_f64());
        let mut output = [0.];
        staged
            .evaluate(&[2.], &mut output, Rational::to_f64)
            .unwrap();
        assert_eq!(
            output,
            [14. / 8.],
            "repeated carriers retain their multiplicities"
        );
        let mut inactive = build()
            .build_parameter_rows(&rows.parameters, &rows.rows[5..])
            .unwrap()
            .map_coeff(&|value| value.re.to_f64());
        inactive
            .evaluate(&[0.], &mut output, Rational::to_f64)
            .unwrap();
        assert_eq!(output, [0.], "absent rows do not prepare the pole at zero");
    }

    #[test]
    fn coefficient_rows_keep_exact_inputs_and_drop_only_cancelled_carriers() {
        initialisation::test_initialise().unwrap();
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
        initialisation::test_initialise().unwrap();
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
                    lhs: family.call_args(iter::once(Atom::num(tag)).chain(args.iter().cloned())),
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
        initialisation::test_initialise().unwrap();
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
            parameters: vec![coefficient, selector],
            rows: vec![
                vec![Rational::from(2), Rational::zero()],
                vec![Rational::from(7), Rational::from(1)],
            ],
        };
        let selected = catalog.specialize(&source, &catalog.rows[0], 0);
        assert_eq!(selected.get_aliases().len(), 2);
        assert!(
            selected
                .get_aliases()
                .keys()
                .all(|key| !source.get_aliases().contains_key(key)),
            "selected rows have independent library-owned alias scopes"
        );
        assert_eq!(selected.into_inner(), 4 * x);
        let selected = catalog.specialize(&source, &catalog.rows[1], 1);
        assert_eq!(selected.get_aliases().len(), 1);
        assert_eq!(selected.into_inner(), Atom::num(1));
        assert_eq!(source.get_aliases().len(), 3);
    }
}
