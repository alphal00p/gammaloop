//! Scalar-only oracles for retained coefficients across distinct Taylor operations.

use std::{collections::BTreeMap, sync::Arc};

use color_eyre::Result;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    evaluate::InliningPolicy,
    function, symbol,
};

use super::branches::{DirectResidueBranches, DirectResidueKey};
use crate::{
    cff::{CutCFFIndex, expression::OrientationID},
    initialisation::test_initialise,
    integrands::process::param_builder::FnMapEntry,
    utils::GS,
    uv::Integrands,
};

struct ScalarFixture {
    definition: Arc<FnMapEntry>,
}

impl ScalarFixture {
    fn new(body: impl FnOnce(&Atom, &Atom) -> Atom) -> Self {
        let x = symbol!("retained_scalar_test::x"; Scalar);
        let y = symbol!("retained_scalar_test::y"; Scalar);
        let tags = vec![DirectResidueBranches::numerator_scope().1, Atom::num(0)];
        Self {
            definition: Arc::new(FnMapEntry {
                lhs: Integrands::scalar_symbol()
                    .call_args(tags.iter().cloned().chain([Atom::var(x), Atom::var(y)])),
                rhs: body(&Atom::var(x), &Atom::var(y)),
                args: vec![x.into(), y.into()],
                tags,
                inlining: InliningPolicy::Always,
                is_alias: false,
            }),
        }
    }

    fn surface() -> Self {
        Self::new(|x, y| {
            Atom::add_many((1..=4).map(|power| x.pow(-power) * (x + y).pow(-power - 1)))
        })
    }

    fn call(&self, x: Atom, y: Atom) -> Atom {
        Integrands::scalar_symbol().call_args(self.definition.tags.iter().cloned().chain([x, y]))
    }

    fn resolve(&self, call: &Atom) -> Atom {
        call.replace_multiple(&[self.definition.replacement()])
    }

    fn integrands(
        &self,
        roots: impl IntoIterator<Item = (CutCFFIndex, Atom)>,
    ) -> Result<Integrands> {
        Integrands::from_iter(roots).with_scalar_definitions(vec![Arc::clone(&self.definition)])
    }
}

fn taylor(
    branches: &DirectResidueBranches,
    variable: Symbol,
    depth: i64,
) -> Result<DirectResidueBranches> {
    branches.series_preserving_numerators(
        variable,
        Atom::Zero.as_view(),
        depth,
        DirectResidueBranches::numerator_scope().1,
    )
}

fn scalar_equal(actual: &Atom, expected: &Atom) {
    // Every fixture in this module is scalar denominator algebra. There is no
    // graph numerator whose factorization this exact comparison could destroy.
    assert!(
        (actual - expected).together().cancel().is_zero(),
        "scalar Taylor mismatch: actual {actual}, expected {expected}"
    );
}

fn assert_scalar_boundary(branches: &DirectResidueBranches) {
    for (_, integrands) in branches.iter_keys() {
        assert!(
            integrands.numerators().is_empty(),
            "temporary scalar families escaped"
        );
        for (_, root) in integrands.iter() {
            assert!(!root.contains_symbol(Symbol::DERIVATIVE));
            assert!(!root.contains_symbol(symbol!("gammalooprs::uv::numerator_family")));
        }
        for entry in integrands.scalar_definitions() {
            assert!(
                !entry.rhs.contains_symbol(Integrands::scalar_symbol()),
                "scalar definitions must remain flat"
            );
            assert!(!entry.rhs.contains_symbol(Symbol::DERIVATIVE));
        }
    }
}

#[test]
fn retained_scalars_survive_successive_momentum_taylor_operations() -> Result<()> {
    test_initialise()?;
    let t = symbol!("retained_nested_test::t"; Scalar);
    let u = symbol!("retained_nested_test::u"; Scalar);
    let q = Atom::var(symbol!("retained_nested_test::q"; Scalar));
    let p = Atom::var(symbol!("retained_nested_test::p"; Scalar));
    let r = Atom::var(symbol!("retained_nested_test::r"; Scalar));
    let mass = Atom::var(symbol!("retained_nested_test::mass"; Scalar));
    let shift = Atom::var(symbol!("retained_nested_test::shift"; Scalar));
    let energy = function!(GS.on_shell_energy, 7, q.pow(2) + mass.pow(2));
    let other = function!(GS.on_shell_energy, 8, (&q + &p).pow(2) + mass.pow(2));
    let fixture = ScalarFixture::surface();
    let call = fixture.call(energy, GS.wrap_esurface(&(other + shift)));
    let source = fixture.resolve(&call);
    let cut = CutCFFIndex::new_all_none();
    let branches =
        DirectResidueBranches::production(OrientationID(0), fixture.integrands([(cut, call)])?)?;
    let deform_t = |atom: &Atom| {
        atom.replace(q.to_pattern())
            .with((&q + Atom::var(t) * &p).to_pattern())
    };
    let first = taylor(&branches.map(deform_t), t, 1)?;
    assert_scalar_boundary(&first);
    assert!(
        first
            .iter_keys()
            .all(|(_, roots)| !roots.scalar_definitions().is_empty())
    );
    let first = first.map(|atom| atom.replace(t).with(1));
    let deform_u = |atom: &Atom| {
        atom.replace(q.to_pattern())
            .with((&q + Atom::var(u) * &r).to_pattern())
    };
    // The enclosing operation receives calls and shared definitions, not a
    // restored copy of the first operation's scalar coefficients.
    let second = taylor(&first.map(deform_u), u, 1)?;
    assert_scalar_boundary(&second);
    assert!(
        second
            .iter_keys()
            .all(|(_, roots)| !roots.scalar_definitions().is_empty())
    );
    let expected_first = deform_t(&source)
        .series(t, Atom::Zero, 1)?
        .to_atom()
        .replace(t)
        .with(1);
    let expected = deform_u(&expected_first)
        .series(u, Atom::Zero, 1)?
        .to_atom();
    let resolved = second.iter_keys().next().unwrap().1.resolved_scalars()?;
    scalar_equal(resolved.iter().next().unwrap().1, &expected);
    Ok(())
}

#[test]
fn retained_scalars_specialize_valuation_for_each_argument_tuple() -> Result<()> {
    test_initialise()?;
    let t = symbol!("retained_valuation_test::t"; Scalar);
    let delta = Atom::var(t);
    let root = (Atom::one() + &delta).pow(Atom::num((1, 2)));
    let fixture = ScalarFixture::new(|x, y| (x - y).pow(-1));
    let arguments = [
        (root.clone(), Atom::one()),
        (root.clone(), Atom::num(2)),
        (&root - &delta / Atom::num(2), Atom::one()),
        (delta.pow(-1), Atom::one()),
    ];
    let mut expected = BTreeMap::new();
    let mut roots = Vec::new();
    for (index, (x, y)) in arguments.into_iter().enumerate() {
        let cut = CutCFFIndex {
            left_threshold_order: Some(index),
            ..CutCFFIndex::new_all_none()
        };
        let call = fixture.call(x, y);
        let expression = (&call + call.pow(2) * (Atom::one() + &delta)) / delta.pow(2);
        // These four specializations have different leading orders, including
        // cancellation through the linear term and a Laurent-valued argument.
        expected.insert(
            cut,
            fixture
                .resolve(&expression)
                .series(t, Atom::Zero, 0)?
                .to_atom(),
        );
        roots.push((cut, expression));
    }
    let branches = DirectResidueBranches::production(OrientationID(0), fixture.integrands(roots)?)?;
    let actual = taylor(&branches, t, 0)?;
    assert_scalar_boundary(&actual);
    let resolved = actual.iter_keys().next().unwrap().1.resolved_scalars()?;
    for (cut, atom) in resolved.iter() {
        scalar_equal(atom, &expected[cut]);
    }
    Ok(())
}

#[test]
fn retained_scalars_expose_energy_owners_to_enclosing_mass_rearrangement() -> Result<()> {
    test_initialise()?;
    let t = symbol!("retained_owner_test::t"; Scalar);
    let u = symbol!("retained_owner_test::u"; Scalar);
    let q = Atom::var(symbol!("retained_owner_test::q"; Scalar));
    let p = Atom::var(symbol!("retained_owner_test::p"; Scalar));
    let mass = Atom::var(symbol!("retained_owner_test::mass"; Scalar));
    let vacuum = Atom::var(symbol!("retained_owner_test::vacuum_mass"; Scalar));
    let original = function!(GS.on_shell_energy, 7, q.pow(2) + mass.pow(2));
    let rearranged = function!(GS.on_shell_energy, 7, q.pow(2) + vacuum.pow(2));
    let spectator = function!(GS.on_shell_energy, 8, q.pow(2) + mass.pow(2));
    let fixture = ScalarFixture::surface();
    let call = fixture.call(original.clone(), spectator.clone());
    let source = fixture.resolve(&call);
    let deform = |atom: &Atom| {
        atom.replace(q.to_pattern())
            .with((&q + Atom::var(t) * &p).to_pattern())
    };
    let branches = DirectResidueBranches::production(
        OrientationID(0),
        fixture.integrands([(CutCFFIndex::new_all_none(), call)])?,
    )?;
    let first = taylor(&branches.map(deform), t, 1)?.map(|atom| atom.replace(t).with(1));
    assert!(
        first
            .iter_keys()
            .all(|(_, roots)| !roots.scalar_definitions().is_empty())
    );
    let transformed = first.map(|atom| {
        atom.replace(original.to_pattern())
            .with(rearranged.to_pattern())
    });
    let root = transformed
        .iter_keys()
        .next()
        .unwrap()
        .1
        .iter()
        .next()
        .unwrap()
        .1;
    assert!(root.contains_symbol(Integrands::scalar_symbol()));
    assert!(root.contains(rearranged.as_view()));
    assert!(root.contains(spectator.as_view()));
    assert!(!root.contains(original.as_view()));
    let expected = deform(&source)
        .series(t, Atom::Zero, 1)?
        .to_atom()
        .replace(t)
        .with(1)
        .replace(original.to_pattern())
        .with(rearranged.to_pattern());
    let deform_outer = |atom: &Atom| {
        atom.replace(q.to_pattern())
            .with((&q + Atom::var(u) * &p).to_pattern())
    };
    let actual = taylor(&transformed.map(deform_outer), u, 1)?;
    assert_scalar_boundary(&actual);
    let resolved = actual.iter_keys().next().unwrap().1.resolved_scalars()?;
    scalar_equal(
        resolved.iter().next().unwrap().1,
        &deform_outer(&expected).series(u, Atom::Zero, 1)?.to_atom(),
    );
    Ok(())
}

#[test]
fn retained_scalar_ids_remain_distinct_across_residue_rows_and_cuts() -> Result<()> {
    test_initialise()?;
    let t = symbol!("retained_rows_test::t"; Scalar);
    let q = Atom::var(symbol!("retained_rows_test::q"; Scalar));
    let mass = Atom::var(symbol!("retained_rows_test::mass"; Scalar));
    let shift = Atom::var(symbol!("retained_rows_test::shift"; Scalar));
    let fixture = ScalarFixture::surface();
    let mut rows = Vec::new();
    let mut expected = BTreeMap::new();
    for row in 0..2 {
        let mut roots = Vec::new();
        for column in 0..2 {
            let cut = CutCFFIndex {
                left_threshold_order: Some(column),
                ..CutCFFIndex::new_all_none()
            };
            let momentum = &q + Atom::num(2 * row + column) + Atom::var(t);
            let energy = function!(GS.on_shell_energy, 7 + row, momentum.pow(2) + mass.pow(2));
            let call = fixture.call(energy.clone(), GS.wrap_esurface(&(energy + &shift)));
            let root = call / Atom::var(t);
            expected.insert(
                (row, cut),
                fixture.resolve(&root).series(t, Atom::Zero, 0)?.to_atom(),
            );
            roots.push((cut, root));
        }
        rows.push((
            DirectResidueKey::production(OrientationID(row)),
            fixture.integrands(roots)?,
        ));
    }
    let actual = taylor(&DirectResidueBranches::from_keyed(rows)?, t, 0)?;
    assert_scalar_boundary(&actual);
    let mut definitions = BTreeMap::new();
    for (key, integrands) in actual.iter_keys() {
        assert!(!integrands.scalar_definitions().is_empty());
        for entry in integrands.scalar_definitions() {
            if let Some(previous) = definitions.insert(entry.lhs.clone(), entry.as_ref()) {
                assert_eq!(
                    previous,
                    entry.as_ref(),
                    "distinct roots reused a scalar id for different bodies"
                );
            }
        }
        for (cut, atom) in integrands.resolved_scalars()?.iter() {
            scalar_equal(atom, &expected[&(key.selector_host.0, *cut)]);
        }
    }
    assert!(!definitions.is_empty());
    Ok(())
}

#[test]
fn retained_scalars_compose_with_factorized_parameterized_numerators() -> Result<()> {
    test_initialise()?;
    let t = symbol!("retained_mixed_test::t"; Scalar);
    let s = symbol!("retained_mixed_test::s"; Scalar);
    let q = Atom::var(symbol!("retained_mixed_test::q"; Scalar));
    let mass = Atom::var(symbol!("retained_mixed_test::mass"; Scalar));
    let shift = Atom::var(symbol!("retained_mixed_test::shift"; Scalar));
    let a = Atom::var(symbol!("retained_mixed_test::a"; Scalar));
    let b = Atom::var(symbol!("retained_mixed_test::b"; Scalar));
    let n1 = Atom::var(symbol!("retained_mixed_test::n1"; Scalar));
    let n2 = Atom::var(symbol!("retained_mixed_test::n2"; Scalar));
    let numerator = (&a + &b).pow(8);
    let jet = Atom::var(s) + Atom::var(t) * &n1 + Atom::var(t).pow(2) * &n2;
    let tag = DirectResidueBranches::numerator_scope().1;
    let family = symbol!("gammalooprs::uv::numerator_family");
    let definition = Arc::new(FnMapEntry {
        lhs: function!(family, &tag, s),
        rhs: &numerator * &jet,
        args: vec![s.into()],
        tags: vec![tag.clone()],
        inlining: InliningPolicy::Always,
        is_alias: false,
    });
    let energy = function!(
        GS.on_shell_energy,
        7,
        (&q + Atom::var(t)).pow(2) + mass.pow(2)
    );
    let fixture = ScalarFixture::surface();
    let scalar = fixture.call(energy, shift);
    let scalar_body = fixture.resolve(&scalar);
    let cut = CutCFFIndex::new_all_none();
    let branches = DirectResidueBranches::from_keyed(
        [0, 1]
            .into_iter()
            .map(|parameter| {
                let root = &scalar * function!(family, &tag, parameter) / Atom::var(t).pow(3);
                Ok((
                    DirectResidueKey::production(OrientationID(parameter)),
                    fixture
                        .integrands([(cut, root)])?
                        .with_numerators([Arc::clone(&definition)])?,
                ))
            })
            .collect::<Result<Vec<_>>>()?,
    )?;
    let expanded = taylor(&branches, t, 0)?;
    for (key, integrands) in expanded.iter_keys() {
        assert!(!integrands.numerators().is_empty());
        assert!(!integrands.scalar_definitions().is_empty());
        for entry in integrands.numerators() {
            assert!(
                !entry
                    .lhs
                    .contains_symbol(symbol!("gammalooprs::uv::scalar_series_source"))
            );
            assert!(entry.rhs.contains(numerator.as_view()));
            assert!(!entry.rhs.contains_symbol(Symbol::DERIVATIVE));
        }
        let resolved = integrands.resolved()?;
        let actual = resolved.iter().next().unwrap().1;
        // Remove this exact retained numerator factor before scalar algebra;
        // never distribute its eighth power for the comparison.
        assert!(actual.contains(numerator.as_view()));
        let actual = actual.replace(numerator.to_pattern()).with(1);
        let specialized = jet.replace(s).with(Atom::num(key.selector_host.0));
        let expected = (&scalar_body * specialized / Atom::var(t).pow(3))
            .series(t, Atom::Zero, 0)?
            .to_atom();
        scalar_equal(&actual, &expected);
    }
    Ok(())
}

#[test]
fn retained_scalars_reject_nonpolynomial_root_occurrences() -> Result<()> {
    test_initialise()?;
    let t = symbol!("retained_polynomial_test::t"; Scalar);
    let fixture = ScalarFixture::new(|x, y| x + y);
    let call = fixture.call(Atom::one() + Atom::var(t), Atom::one());
    let cut = CutCFFIndex::new_all_none();
    for root in [
        (&call - Atom::num(2)).pow(-1),
        call.pow(Atom::num((1, 2))),
        function!(symbol!("retained_polynomial_test::nonlinear"), &call),
    ] {
        let branches = DirectResidueBranches::production(
            OrientationID(0),
            fixture.integrands([(cut, root)])?,
        )?;
        assert!(
            taylor(&branches, t, 0)
                .unwrap_err()
                .to_string()
                .contains("polynomial use")
        );
    }
    Ok(())
}

#[test]
fn retained_scalars_normalize_large_and_fractional_positive_valuations() -> Result<()> {
    test_initialise()?;
    let t = symbol!("retained_high_valuation_test::t"; Scalar);
    let q = Atom::var(symbol!("retained_high_valuation_test::q"; Scalar));
    let mass = Atom::var(symbol!("retained_high_valuation_test::mass"; Scalar));
    let fixture = ScalarFixture::new(|x, y| x * (y + Atom::one()).pow(-3));
    let cut = CutCFFIndex::new_all_none();
    for (valuation, center) in [
        (Atom::num(14), Atom::Zero),
        (Atom::num((3, 2)), Atom::Zero),
        (Atom::num(14), Atom::num(2)),
    ] {
        let delta = Atom::var(t) - &center;
        let energy = function!(GS.on_shell_energy, 7, (&q + &delta).pow(2) + mass.pow(2));
        let call = fixture.call(delta.pow(&valuation), energy.clone());
        // Despite the body's large absolute valuation, this root needs only
        // three regular coefficients. The fractional case also checks that
        // shifting the endpoint preserves the native series lattice.
        let root = call / delta.pow(&valuation + Atom::num(2));
        let branches = DirectResidueBranches::production(
            OrientationID(0),
            fixture.integrands([(cut, root)])?,
        )?;
        let actual = branches.series_preserving_numerators(
            t,
            center.as_view(),
            0,
            DirectResidueBranches::numerator_scope().1,
        )?;
        assert_scalar_boundary(&actual);
        let resolved = actual.iter_keys().next().unwrap().1.resolved_scalars()?;
        // Cancel the known monomial analytically before building the oracle;
        // computing a depth-16 native series would repeat the performance bug.
        let expected = ((energy + Atom::one()).pow(-3) / delta.pow(2))
            .series(t, center.as_view(), 0)?
            .to_atom();
        scalar_equal(resolved.iter().next().unwrap().1, &expected);
    }
    Ok(())
}
