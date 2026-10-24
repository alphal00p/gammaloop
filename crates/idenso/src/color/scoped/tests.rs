use super::*;
use crate::{
    IndexTooling,
    representations::{ColorAdjoint, ColorFundamental},
};
use spenso::structure::{abstract_index::AbstractIndex, slot::IsAbstractSlot};
use symbolica::{function, symbol};

type Tensor = SymbolicTensor<PartialStructure>;

fn family(name: &str) -> ScopedColorFamily {
    crate::representations::initialize();
    ColorFamilyDeclaration::su_n(
        LibraryRep::new_dual(&format!("scoped_color::{name}::fund")).unwrap(),
        LibraryRep::new_self_dual(&format!("scoped_color::{name}::adj")).unwrap(),
        [
            format!("scoped_color::{name}::t"),
            format!("scoped_color::{name}::f"),
        ],
        symbol!(&format!("scoped_color::{name}::Nc"); Real),
    )
    .unwrap()
    .register()
    .unwrap()
}

fn native_generator(a: usize, i: usize, j: usize) -> Atom {
    let adj = LibraryRep::from(ColorAdjoint {}).new_rep(8);
    let fund = LibraryRep::from(ColorFundamental {}).new_rep(3);
    CS.explicit_t(
        adj.slot::<AbstractIndex, _>(AbstractIndex::Normal(a))
            .to_atom(),
        fund.slot::<AbstractIndex, _>(AbstractIndex::Normal(i))
            .to_atom(),
        fund.dual()
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(j))
            .to_atom(),
    )
}

fn native_trace() -> Tensor {
    Tensor::infer(native_generator(1, 2, 3) * native_generator(4, 3, 2)).unwrap()
}

fn native_reduce(value: &Tensor) -> Tensor {
    value
        .simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
}

#[test]
fn scoped_color_equal_dimension_families_reduce_independently() {
    let a = family("reduction_a");
    let b = family("reduction_b");
    let native = native_trace();
    let first = native.rename_symbols(&a.native_to_scoped()).unwrap();
    let second = native.rename_symbols(&b.native_to_scoped()).unwrap();
    let both = first.multiply(&second).unwrap();
    assert_eq!(both.rank(), 4);
    let result = b
        .reduce(
            &a.reduce(&both, ColorSimplifySettings::default()).unwrap(),
            ColorSimplifySettings::default(),
        )
        .unwrap();
    let expected = native_reduce(&native)
        .rename_symbols(&a.native_to_scoped())
        .unwrap()
        .multiply(
            &native_reduce(&native)
                .rename_symbols(&b.native_to_scoped())
                .unwrap(),
        )
        .unwrap();
    assert_eq!(result.rank(), 4);
    assert_eq!(
        result.canonize().unwrap().expression(),
        expected.canonize().unwrap().expression()
    );
}

#[test]
fn scoped_color_conjugates_mixed_families_and_scalar_once() {
    let a = family("adjoint_a");
    let b = family("adjoint_b");
    let native = Tensor::infer(native_generator(1, 2, 3)).unwrap();
    let first = native.rename_symbols(&a.native_to_scoped()).unwrap();
    let second = native.rename_symbols(&b.native_to_scoped()).unwrap();
    let source = Tensor::infer(Atom::i() * first.expression() * second.expression()).unwrap();
    let expected_first = native
        .dirac_adjoint(true)
        .unwrap()
        .rename_symbols(&a.native_to_scoped())
        .unwrap();
    let expected_second = native
        .dirac_adjoint(true)
        .unwrap()
        .rename_symbols(&b.native_to_scoped())
        .unwrap();
    let expected =
        Tensor::infer(-Atom::i() * expected_first.expression() * expected_second.expression())
            .unwrap();
    let adjoint = source.dirac_adjoint(true).unwrap();
    assert_eq!(adjoint.expression(), expected.expression());
    assert_eq!(
        adjoint.dirac_adjoint(true).unwrap().expression(),
        source.expression()
    );
}

#[test]
fn scoped_color_compact_traces_use_native_conjugation() {
    let a = family("trace_adjoint_a");
    let b = family("trace_adjoint_b");
    let native = native_trace().contract(Default::default()).unwrap();
    let first = native.rename_symbols(&a.native_to_scoped()).unwrap();
    let second = native.rename_symbols(&b.native_to_scoped()).unwrap();
    let source = first.multiply(&second).unwrap();
    let native_adjoint = native.dirac_adjoint(true).unwrap();
    let expected = native_adjoint
        .rename_symbols(&a.native_to_scoped())
        .unwrap()
        .multiply(
            &native_adjoint
                .rename_symbols(&b.native_to_scoped())
                .unwrap(),
        )
        .unwrap();
    assert_eq!(
        source.dirac_adjoint(true).unwrap().expression(),
        expected.expression()
    );
}

#[test]
fn scoped_color_dimension_and_casimir_payloads_keep_group_identity() {
    let a = family("symbolic_a");
    let b = family("symbolic_b");
    let fund = LibraryRep::from(ColorFundamental {}).new_rep(CS.nc);
    let adj_dimension = symbol!("scoped_color::native_adj_dimension");
    let adj = LibraryRep::from(ColorAdjoint {}).new_rep(adj_dimension);
    let generator = |a, i, j| {
        CS.explicit_t(
            adj.slot::<AbstractIndex, _>(AbstractIndex::Normal(a))
                .to_atom(),
            fund.slot::<AbstractIndex, _>(AbstractIndex::Normal(i))
                .to_atom(),
            fund.dual()
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(j))
                .to_atom(),
        )
    };
    let native = Tensor::infer(generator(1, 2, 3) * generator(4, 3, 2)).unwrap();
    let mut first_map = a.native_to_scoped();
    first_map.insert(adj_dimension, symbol!("scoped_color::symbolic_a::dim_adj"));
    let mut second_map = b.native_to_scoped();
    second_map.insert(adj_dimension, symbol!("scoped_color::symbolic_b::dim_adj"));
    let both = native
        .rename_symbols(&first_map)
        .unwrap()
        .multiply(&native.rename_symbols(&second_map).unwrap())
        .unwrap();
    let result = b
        .reduce(
            &a.reduce(&both, ColorSimplifySettings::default()).unwrap(),
            ColorSimplifySettings::default(),
        )
        .unwrap();
    let expected = native_reduce(&native)
        .rename_symbols(&first_map)
        .unwrap()
        .multiply(&native_reduce(&native).rename_symbols(&second_map).unwrap())
        .unwrap();
    assert_eq!(
        result.canonize().unwrap().expression(),
        expected.canonize().unwrap().expression()
    );
    assert!(
        result
            .expression()
            .contains(Atom::var(a.declaration.nc).as_view())
    );
    assert!(
        result
            .expression()
            .contains(Atom::var(b.declaration.nc).as_view())
    );
    let cas = CS.cas(2, fund.to_symbolic([]));
    assert_ne!(
        cas.rename_symbols(&first_map),
        cas.rename_symbols(&second_map)
    );
}

#[test]
fn scoped_color_remaps_declaration_and_expression_together() {
    let a = family("remap_source");
    let b = family("remap_target");
    let symbols = a.symbol_map_to(&b).unwrap();
    let declaration = a
        .declaration
        .remapped(&symbols, b.declaration.heads.clone())
        .unwrap();
    assert_eq!(declaration, b.declaration);
    let restored = ColorFamilyDeclaration::from_user_data(&declaration.to_user_data())
        .unwrap()
        .register()
        .unwrap();
    assert_eq!(restored, b);
    let value = native_trace()
        .rename_symbols(&a.native_to_scoped())
        .unwrap();
    let mapped = value.rename_symbols(&symbols).unwrap();
    assert_eq!(
        mapped.expression(),
        native_trace()
            .rename_symbols(&b.native_to_scoped())
            .unwrap()
            .expression()
    );
    assert_eq!(ScopedColorFamily::from_symbol(b.heads[0]).unwrap(), Some(b));
}

#[test]
fn scoped_color_native_collision_and_unsupported_action_diagnostics() {
    let a = family("diagnostics");
    let native = native_trace();
    let scoped = native.rename_symbols(&a.native_to_scoped()).unwrap();
    let collision = native.multiply(&scoped).unwrap();
    assert!(
        a.reduce(&collision, ColorSimplifySettings::default())
            .unwrap_err()
            .contains("existing unmapped symbol")
    );
    let dimension_collision = Tensor::infer(scoped.expression() * Atom::var(CS.nc)).unwrap();
    assert!(
        a.reduce(&dimension_collision, ColorSimplifySettings::default())
            .unwrap_err()
            .contains("existing unmapped symbol")
    );
    let foreign = LibraryRep::new_dual("scoped_color::unassigned").unwrap();
    assert!(
        a.validate_representation(foreign)
            .unwrap_err()
            .contains("explicit native generator")
    );
    let bad = function!(a.heads[0], 1, 2, 3);
    assert!(
        Tensor::infer(bad)
            .unwrap_err()
            .to_string()
            .contains("native slots")
    );
    let mut metadata = a.heads[0].get_data().clone();
    let UserData::List(fields) = &mut metadata else {
        unreachable!()
    };
    fields[0] = UserData::Integer(5);
    let wrong = SymbolBuilder::new(NamespacedSymbol::parse("scoped_color::invalid_role"))
        .with_attributes(a.heads[0].get_attributes().to_vec())
        .with_tags(a.heads[0].get_tags())
        .with_user_data(metadata)
        .build()
        .unwrap();
    assert!(
        ScopedColorFamily::from_symbol(wrong)
            .unwrap_err()
            .contains("role")
    );
}

#[test]
fn scoped_color_native_signatures_validate_dimensions_and_compact_channels() {
    let a = family("signature_a");
    let b = family("signature_b");
    let adjoint = a.declaration.adjoint.new_rep(8);
    let slot = |rep: Representation<LibraryRep>, index| {
        rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    };
    assert!(
        a.structure_constant(
            slot(adjoint, 1),
            slot(a.declaration.adjoint.new_rep(7), 2),
            slot(adjoint, 3),
        )
        .is_err()
    );
    assert!(
        Tensor::infer(function!(
            a.structure_symbol(),
            slot(adjoint, 1),
            Atom::var(SPENSO_TAG.chain_in),
            slot(adjoint, 2)
        ))
        .is_err()
    );
    let compact = |family: &ScopedColorFamily, reverse: bool| {
        function!(
            family.generator_symbol(),
            slot(family.declaration.adjoint.new_rep(8), 1),
            Atom::var(if reverse {
                SPENSO_TAG.chain_out
            } else {
                SPENSO_TAG.chain_in
            }),
            Atom::var(if reverse {
                SPENSO_TAG.chain_in
            } else {
                SPENSO_TAG.chain_out
            })
        )
    };
    let channel = a.declaration.fundamental.new_rep(3).to_symbolic([]);
    assert!(Tensor::infer(spenso::shadowing::trace(&channel, [compact(&a, false)])).is_ok());
    assert!(Tensor::infer(spenso::shadowing::trace(&channel, [compact(&b, false)])).is_err());
    assert!(Tensor::infer(spenso::shadowing::trace(channel, [compact(&a, true)])).is_err());
}

#[test]
fn scoped_color_explicit_su_n_normalization_keeps_fixed_ranks_exact() {
    let family = family("normalization");
    for n in [2, 3] {
        let adj = family.declaration.adjoint.new_rep(n * n - 1);
        let fund = family.declaration.fundamental.new_rep(n);
        let slot = |rep: Representation<LibraryRep>, index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        };
        let a = slot(adj, 1);
        let b = slot(adj, 4);
        let first = family
            .fundamental_generator(a.clone(), slot(fund, 2), slot(fund.dual(), 3))
            .unwrap();
        let second = family
            .fundamental_generator(b.clone(), slot(fund, 3), slot(fund.dual(), 2))
            .unwrap();
        let trace = first.multiply(&second).unwrap();
        let reduced = family
            .reduce(
                &trace,
                ColorSimplifySettings::default().with_cof_dimension_invariants(),
            )
            .unwrap();
        let expected = Tensor::infer(
            spenso::network::library::symbolic::ETS.metric_literal(a, b) / Atom::num(2),
        )
        .unwrap();
        assert_eq!(
            reduced.canonize().unwrap().expression(),
            expected.canonize().unwrap().expression()
        );
        assert!(!reduced.expression().contains_symbol(CS.nc));
        assert!(!reduced.expression().contains_symbol(family.declaration.nc));
    }
}

#[test]
fn scoped_color_dual_and_adjoint_generator_signs_follow_native_conventions() {
    let family = family("generator_actions");
    let fund = family.declaration.fundamental;
    let adj = family.declaration.adjoint;
    let slot = |rep: LibraryRep, dimension, index| {
        rep.new_rep(dimension)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    };
    let a = slot(adj, 3, 1);
    let i = slot(fund, 2, 2);
    let j = slot(fund.dual(), 2, 3);
    let direct = family
        .generator(fund, a.clone(), i.clone(), j.clone())
        .unwrap();
    let dual = family.generator(fund.dual(), a.clone(), j, i).unwrap();
    assert_eq!(dual.expression(), &(-direct.expression()));
    let b = slot(adj, 3, 2);
    let c = slot(adj, 3, 3);
    let action = family
        .generator(adj, a.clone(), b.clone(), c.clone())
        .unwrap();
    let transpose = family
        .generator(adj, a.clone(), c.clone(), b.clone())
        .unwrap();
    let structure = family.structure_constant(a, b, c).unwrap();
    assert_eq!(action.expression(), &(-Atom::i() * structure.expression()));
    assert_eq!(action.expression().spenso_conj(), *transpose.expression());
}

#[test]
fn scoped_color_bifundamental_actions_keep_both_spaces_and_conjugate_once() {
    let a = family("bifundamental_a");
    let b = family("bifundamental_b");
    let native = Tensor::infer(native_generator(1, 2, 3)).unwrap();
    let first = native.rename_symbols(&a.native_to_scoped()).unwrap();
    let second = native.rename_symbols(&b.native_to_scoped()).unwrap();
    let slot = |rep: LibraryRep| {
        rep.new_rep(3)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(3))
            .to_atom()
    };
    let phi = Tensor::infer(function!(
        spenso::tensor_symbol!("scoped_color::bifundamental_phi"),
        slot(a.declaration.fundamental),
        slot(b.declaration.fundamental)
    ))
    .unwrap();
    let source = first.multiply(&second).unwrap().multiply(&phi).unwrap();
    assert_eq!(source.rank(), 4);
    let expected = first
        .dirac_adjoint(true)
        .unwrap()
        .multiply(&second.dirac_adjoint(true).unwrap())
        .unwrap()
        .multiply(&phi.dirac_adjoint(true).unwrap())
        .unwrap();
    assert_eq!(
        source
            .dirac_adjoint(true)
            .unwrap()
            .canonize()
            .unwrap()
            .expression(),
        expected.canonize().unwrap().expression()
    );
    assert_eq!(
        source
            .dirac_adjoint(true)
            .unwrap()
            .dirac_adjoint(true)
            .unwrap()
            .expression(),
        source.expression()
    );
}

#[test]
fn scoped_color_three_generator_traces_keep_scalar_spectators_under_adjoint() {
    let a = family("three_trace_a");
    let b = family("three_trace_b");
    let native = Tensor::infer(
        native_generator(1, 4, 5) * native_generator(2, 5, 6) * native_generator(3, 6, 4),
    )
    .unwrap()
    .contract(Default::default())
    .unwrap();
    let first = native.rename_symbols(&a.native_to_scoped()).unwrap();
    let second = native.rename_symbols(&b.native_to_scoped()).unwrap();
    let coefficient = Atom::i()
        * Atom::var(symbol!("scoped_color::complex_coefficient"))
        * CS.cas(2, b.declaration.fundamental.new_rep(3).to_symbolic([]));
    let source = Tensor::infer(&coefficient * first.expression() * second.expression()).unwrap();
    let expected = Tensor::infer(
        coefficient.spenso_conj()
            * first.dirac_adjoint(true).unwrap().expression()
            * second.dirac_adjoint(true).unwrap().expression(),
    )
    .unwrap();
    assert_eq!(
        source
            .dirac_adjoint(true)
            .unwrap()
            .canonize()
            .unwrap()
            .expression(),
        expected.canonize().unwrap().expression()
    );
    let reduced = b
        .reduce(
            &a.reduce(&source, ColorSimplifySettings::default()).unwrap(),
            ColorSimplifySettings::default(),
        )
        .unwrap();
    assert_eq!(
        reduced
            .dirac_adjoint(true)
            .unwrap()
            .dirac_adjoint(true)
            .unwrap()
            .canonize()
            .unwrap()
            .expression(),
        reduced.canonize().unwrap().expression()
    );
}

#[test]
fn scoped_color_malformed_head_is_diagnosed_without_disabling_valid_raw_families() {
    let a = family("malformed_spectator_a");
    let source = Tensor::infer(native_generator(1, 2, 3))
        .unwrap()
        .rename_symbols(&a.native_to_scoped())
        .unwrap();
    let malformed =
        SymbolBuilder::new(NamespacedSymbol::parse("scoped_color::malformed_spectator"))
            .with_tags(vec![
                TAG.into(),
                spenso::network::tags::SPENSO_TAG.tensor.clone(),
            ])
            .with_user_data(UserData::None)
            .build()
            .unwrap();
    let malformed = function!(malformed, 1);
    let combined = source.expression() * &malformed;
    assert!(Tensor::infer(combined.clone()).is_err());
    assert_eq!(
        combined.spenso_conj(),
        source.expression().spenso_conj() * malformed.spenso_conj()
    );
}

#[test]
fn scoped_color_invariant_only_archive_persists_its_explicit_declaration() {
    use base64::Engine;
    let a = family("invariant_archive");
    let reduced = a
        .reduce(
            &native_trace()
                .rename_symbols(&a.native_to_scoped())
                .unwrap(),
            ColorSimplifySettings::default(),
        )
        .unwrap();
    assert!(!reduced.expression().contains_symbol(a.generator_symbol()));
    // Shared invariant heads retain their scoped representation arguments.
    // A Theory archive stores the declaration explicitly, without modifying
    // Spenso-owned immutable representation metadata or relying on a registry.
    let carrier = SymbolBuilder::new(NamespacedSymbol::parse(
        "scoped_color::archived_declaration",
    ))
    .with_user_data(a.declaration.to_user_data())
    .build()
    .unwrap();
    let payload = function!(
        symbol!("scoped_color::archive_payload"),
        Atom::var(carrier),
        reduced.expression()
    );
    let mut bytes = Vec::new();
    payload.export(&mut bytes).unwrap();
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "color::scoped::tests::scoped_color_archive_child",
            "--ignored",
            "--nocapture",
        ])
        .env(
            "IDENSO_SCOPED_COLOR_ARCHIVE",
            base64::engine::general_purpose::STANDARD.encode(bytes),
        )
        .env("IDENSO_SCOPED_COLOR_INVARIANT_ARCHIVE", "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
#[ignore = "fresh-process archive helper"]
fn scoped_color_archive_child() {
    use base64::Engine;
    let payload = base64::engine::general_purpose::STANDARD
        .decode(std::env::var("IDENSO_SCOPED_COLOR_ARCHIVE").unwrap())
        .unwrap();
    let value = Atom::import(&mut payload.as_slice(), None).unwrap();
    if std::env::var_os("IDENSO_SCOPED_COLOR_INVARIANT_ARCHIVE").is_some() {
        let AtomView::Fun(payload) = value.as_view() else {
            panic!("archive payload")
        };
        let args = payload.iter().collect::<Vec<_>>();
        let AtomView::Var(carrier) = args[0] else {
            panic!("declaration carrier")
        };
        let family = ColorFamilyDeclaration::from_user_data(carrier.get_symbol().get_data())
            .unwrap()
            .register()
            .unwrap();
        let source = Tensor::infer(args[1].to_owned()).unwrap();
        assert_eq!(
            family
                .reduce(&source, ColorSimplifySettings::default())
                .unwrap()
                .expression(),
            source.expression()
        );
        assert_eq!(
            source
                .dirac_adjoint(true)
                .unwrap()
                .dirac_adjoint(true)
                .unwrap()
                .expression(),
            source.expression()
        );
        return;
    }
    let source = Tensor::infer(value).unwrap();
    let families = ScopedColorFamily::in_expression(source.expression().as_view()).unwrap();
    assert_eq!(families.len(), 1);
    let reduced = families[0]
        .reduce(&source, ColorSimplifySettings::default())
        .unwrap();
    assert_ne!(source.expression(), reduced.expression());
    assert_eq!(
        source
            .dirac_adjoint(true)
            .unwrap()
            .dirac_adjoint(true)
            .unwrap()
            .expression(),
        source.expression()
    );
}

#[test]
fn scoped_color_fresh_process_archive_restores_native_action() {
    use base64::Engine;
    let a = family("archive");
    let source = native_trace()
        .rename_symbols(&a.native_to_scoped())
        .unwrap();
    let mut payload = Vec::new();
    source.expression().export(&mut payload).unwrap();
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "color::scoped::tests::scoped_color_archive_child",
            "--ignored",
            "--nocapture",
        ])
        .env(
            "IDENSO_SCOPED_COLOR_ARCHIVE",
            base64::engine::general_purpose::STANDARD.encode(payload),
        )
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn scoped_color_symmetric_invariant_representation_is_metadata() {
    let family = family("symmetric_metadata");
    let representation = family.declaration.fundamental.new_rep(3);
    let adjoint = family.declaration.adjoint.new_rep(8);
    let slots = (1..=3)
        .map(|index| {
            adjoint
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        })
        .collect();
    let source = Tensor::infer(CS.symmetric_d(representation.to_symbolic([]), slots)).unwrap();
    assert_eq!(source.rank(), 3, "the leading representation is metadata");
    assert_eq!(source.canonize().unwrap().rank(), 3);
    assert_eq!(source.contract(Default::default()).unwrap().rank(), 3);
    assert_eq!(source.dirac_adjoint(true).unwrap().rank(), 3);
    let slot = |index| {
        adjoint
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    };
    let metric =
        Tensor::infer(spenso::network::library::symbolic::ETS.metric_literal(slot(3), slot(4)))
            .unwrap();
    let contracted = source
        .multiply(&metric)
        .unwrap()
        .contract(Default::default())
        .unwrap();
    let expected = Tensor::infer(CS.symmetric_d(
        representation.to_symbolic([]),
        vec![slot(1), slot(2), slot(4)],
    ))
    .unwrap();
    assert_eq!(
        contracted.canonize().unwrap().expression(),
        expected.canonize().unwrap().expression()
    );
    let formal = CS.symmetric_d(
        Atom::var(symbol!("scoped_color::formal_symmetric_label")),
        vec![slot(1), slot(2), slot(3)],
    );
    let AtomView::Fun(function) = formal.as_view() else {
        unreachable!()
    };
    assert!(CS.invariant_representations(function).is_err());
    assert_eq!(Tensor::infer(formal).unwrap().rank(), 3);
}

#[test]
fn scoped_color_native_invariant_metadata_validates_arity_and_compact_payloads() {
    let owner = family("invariant_metadata");
    let fundamental = owner.declaration.fundamental.new_rep(3);
    let other = family("invariant_metadata_other");
    let other_fundamental = other.declaration.fundamental.new_rep(3);
    let metadata = |expression: &Atom| {
        let AtomView::Fun(function) = expression.as_view() else {
            unreachable!()
        };
        CS.invariant_representations(function)
    };
    let degree = Atom::var(symbol!("scoped_color::invariant_degree"));
    for invariant in [
        CS.cas(&degree, fundamental.to_symbolic([])),
        CS.idx(2, fundamental.to_symbolic([])),
    ] {
        assert_eq!(metadata(&invariant).unwrap(), Some(vec![fundamental]));
        assert_eq!(Tensor::infer(invariant).unwrap().rank(), 0);
    }
    let mixed = CS.gram(
        4,
        fundamental.to_symbolic([]),
        other_fundamental.to_symbolic([]),
    );
    assert_eq!(
        metadata(&mixed).unwrap(),
        Some(vec![fundamental, other_fundamental])
    );
    let wrong_dimension = owner.declaration.fundamental.new_rep(999);
    let invariant = CS.cas(2, wrong_dimension.to_symbolic([]));
    assert_eq!(metadata(&invariant).unwrap(), Some(vec![wrong_dimension]));
    let formal = Atom::var(symbol!("scoped_color::formal_invariant_label"));
    for invariant in [
        CS.cas(2, &formal),
        CS.gram(4, &formal, other_fundamental.to_symbolic([])),
    ] {
        assert!(metadata(&invariant).is_err());
        assert_eq!(Tensor::infer(invariant).unwrap().rank(), 0);
    }
    for invariant in [
        function!(CS.cas, fundamental.to_symbolic([])),
        function!(CS.cas, 2, fundamental.to_symbolic([]), 99),
        function!(CS.idx, 2),
        function!(CS.gram, 4, fundamental.to_symbolic([])),
        CS.cas(
            2,
            fundamental
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(1))
                .to_atom(),
        ),
        CS.cas(
            fundamental
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(1))
                .to_atom(),
            fundamental.to_symbolic([]),
        ),
        CS.cas(
            function!(
                spenso::tensor_symbol!("scoped_color::invariant_tensor_degree"),
                fundamental
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(1))
                    .to_atom()
            ),
            fundamental.to_symbolic([]),
        ),
    ] {
        assert!(metadata(&invariant).is_err(), "{invariant}");
        assert!(Tensor::infer(invariant.clone()).is_err(), "{invariant}");
    }
}
