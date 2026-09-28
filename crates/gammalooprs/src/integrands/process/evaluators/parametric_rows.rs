use std::collections::HashMap;

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use spenso::algebra::complex::Complex;
use symbolica::prelude::{
    AliasedAtom, Atom, AtomCore, FunctionBuilder, Rational, Replacement, SingleFloat, Symbol,
    function, symbol,
};

use crate::{
    GammaLoopContext,
    integrands::evaluation::EvaluationMetaData,
    utils::{F, FloatLike, hyperdual_utils::DualOrNot},
    uv::ParametricIntegrandTerm,
};

use super::{GenericEvaluator, evaluate_evaluator};

/// Exact evaluator-local inputs. Only the single-parametric program receives
/// these appended ports; all other programs specialize them after preprocessing.
#[derive(Clone, Debug, PartialEq, Eq, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct ParametricResidueRows {
    pub(crate) parameters: Vec<Atom>,
    pub(crate) rows: Vec<Vec<Rational>>,
}

impl ParametricResidueRows {
    fn select(values: &[Atom], bits: &[Atom]) -> Atom {
        if values.len() <= 1 {
            return values.first().cloned().unwrap_or(Atom::Zero);
        }
        let (bit, remaining) = bits.split_last().expect("insufficient dispatch bits");
        let split = (1usize << remaining.len()).min(values.len());
        Symbol::IF.call_args([
            bit.clone(),
            Self::select(&values[split..], remaining),
            Self::select(&values[..split], remaining),
        ])
    }

    pub(super) fn lower(
        integrand: &Atom,
        terms: &[ParametricIntegrandTerm],
    ) -> Result<(Vec<Atom>, Self)> {
        // Physical maps, cuts and threshold localization are complete here.
        // Only now may equal full coefficient rows share their scalar sum;
        // earlier UV owners retain each carrier's native map provenance.
        let term_rows = terms
            .iter()
            .map(|term| {
                let mut grouped = std::collections::BTreeMap::<Vec<Rational>, Vec<&Atom>>::new();
                for row in &term.rows {
                    grouped
                        .entry(row.coefficients.clone())
                        .or_default()
                        .push(&row.carrier);
                }
                grouped
                    .into_iter()
                    .filter_map(|(coefficients, carriers)| {
                        let carrier = Atom::add_many(carriers);
                        (!carrier.is_zero()).then_some((coefficients, carrier))
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let coefficient_count = terms
            .iter()
            .map(|term| term.numerator.args.len())
            .max()
            .unwrap_or(0);
        let row_count = term_rows.iter().map(Vec::len).max().unwrap_or(1).max(1);
        let has_scalar = integrand != &Atom::Zero;
        let term_count = terms.len() + usize::from(has_scalar);
        let bit_count =
            |count: usize| usize::BITS as usize - count.saturating_sub(1).leading_zeros() as usize;
        let row_bits = bit_count(row_count);
        let term_bits = bit_count(term_count);
        let coefficient = symbol!("gammalooprs::parametric_residue::coefficient"; Scalar);
        let selector = symbol!("gammalooprs::parametric_residue::selector"; Scalar);
        // Different additive terms reuse the same port positions. Their lazy
        // dispatch gives each position its term-local meaning and bounds arity
        // by the widest product, rather than the sum of all forest signatures.
        let parameters = (0..coefficient_count)
            .map(|i| function!(coefficient, i))
            .chain((0..row_bits + term_bits).map(|i| function!(selector, i)))
            .collect::<Vec<_>>();
        let mut atoms = Vec::with_capacity(term_count);
        let mut rows = Vec::new();
        for (term_id, (term, source_rows)) in terms.iter().zip(&term_rows).enumerate() {
            let replacements = term
                .numerator
                .args
                .iter()
                .zip(&parameters)
                .map(|(formal, port)| {
                    Replacement::new(Atom::from(formal.clone()).to_pattern(), port.clone())
                })
                .collect::<Vec<_>>();
            let carriers = source_rows
                .iter()
                .map(|(_, carrier)| carrier.replace_multiple(&replacements))
                .collect::<Vec<_>>();
            atoms.push(
                term.numerator.lhs.replace_multiple(&replacements)
                    * Self::select(
                        &carriers,
                        &parameters[coefficient_count..coefficient_count + row_bits],
                    ),
            );
            for (row_id, (coefficients, _)) in source_rows.iter().enumerate() {
                if coefficients.len() != term.numerator.args.len() {
                    return Err(eyre!(
                        "parametric numerator row has {} coefficients for {} formals",
                        coefficients.len(),
                        term.numerator.args.len()
                    ));
                }
                let mut values = coefficients.clone();
                values.resize(coefficient_count, Rational::zero());
                values
                    .extend((0..row_bits).map(|bit| Rational::from(((row_id >> bit) & 1) as i64)));
                values.extend(
                    (0..term_bits).map(|bit| Rational::from(((term_id >> bit) & 1) as i64)),
                );
                rows.push(values);
            }
        }
        if has_scalar {
            atoms.push(integrand.clone());
            let mut values = vec![Rational::zero(); coefficient_count + row_bits];
            values.extend(
                (0..term_bits).map(|bit| Rational::from(((terms.len() >> bit) & 1) as i64)),
            );
            rows.push(values);
        }
        // A cut whose terms all vanish still has one executable zero, never an
        // empty reduction or an extra contribution from another source.
        if rows.is_empty() {
            atoms = vec![Atom::Zero];
            rows.push(vec![Rational::zero(); parameters.len()]);
        }
        Ok((atoms, Self { parameters, rows }))
    }

    pub(super) fn join_preprocessed(&self, atoms: Vec<AliasedAtom>) -> AliasedAtom {
        let bit_count =
            usize::BITS as usize - atoms.len().saturating_sub(1).leading_zeros() as usize;
        let roots = atoms
            .iter()
            .map(|atom| atom.get_root().clone())
            .collect::<Vec<_>>();
        let mut joined = AliasedAtom::from(Self::select(
            &roots,
            &self.parameters[self.parameters.len() - bit_count..],
        ));
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
        let mut replacements = self
            .parameters
            .iter()
            .zip(row)
            .map(|(parameter, value)| {
                Replacement::new(parameter.to_pattern(), Atom::num(value.clone()))
            })
            .collect::<Vec<_>>();
        let aliases = atom
            .get_aliases()
            .keys()
            .map(|alias| {
                (
                    alias.clone(),
                    FunctionBuilder::from_atom(alias).add_arg(row_id).finish(),
                )
            })
            .collect::<HashMap<_, _>>();
        replacements.extend(
            aliases
                .iter()
                .map(|(old, new)| Replacement::new(old.to_pattern(), new.clone())),
        );
        let mut result = AliasedAtom::from(atom.get_root().replace_multiple(&replacements));
        for (alias, body) in atom.get_aliases() {
            result.register_alias(aliases[alias].clone(), body.replace_multiple(&replacements));
        }
        result.prune();
        result
    }

    pub(super) fn evaluate<T: FloatLike>(
        &self,
        evaluator: &mut GenericEvaluator,
        physical_input: &[Complex<F<T>>],
        metadata: &mut EvaluationMetaData,
    ) -> Vec<DualOrNot<Complex<F<T>>>> {
        let shape_len = evaluator.dual_shape.as_ref().map_or(1, Vec::len);
        let prototype = physical_input
            .first()
            .map_or_else(|| F(T::new_zero()), |value| value.re.clone());
        let mut input = physical_input.to_vec();
        input.resize(
            input.len() + self.parameters.len() * shape_len,
            Complex::new_re(prototype.zero()),
        );
        let mut result: Option<Vec<DualOrNot<Complex<F<T>>>>> = None;
        for row in &self.rows {
            for (column, value) in row.iter().enumerate() {
                input[physical_input.len() + column * shape_len] =
                    Complex::new_re(prototype.from_rational(value));
            }
            let output = evaluate_evaluator(evaluator, &input, metadata);
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
