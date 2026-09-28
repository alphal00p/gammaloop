use crate::representations::ColorSextet;
use insta::assert_snapshot;
use spenso::network::parsing::{ParseSettings, StructureFromAtom};
use spenso::network::tags::SPENSO_TAG;
use spenso::shadowing::{Collectable, TensorCollectExt, TensorCollectFilter};
use spenso::structure::{Canonicalized, IndexlessNamedStructure, TensorStructure};
use symbolica_utils::AtomPrintExt;

static _CF: LazyLock<Canonicalized<IndexlessNamedStructure<Symbol, ()>>> = LazyLock::new(|| {
    IndexlessNamedStructure::from_iter(
        [
            ColorAdjoint {}.new_rep(8),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(4),
            ColorAdjoint {}.new_rep(2),
        ],
        CS.f,
        None,
    )
});

use spenso::{antisym, chain, network::parsing::ShadowedStructure, s, slot, sym, trace};
use symbolica::{function, parse, parse_lit};

use crate::shorthands::schoonschip::Schoonschip;
use crate::tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor};
use crate::{Cookable, IndexTooling};
use crate::{color_cas, color_f, color_gram, color_idx, color_str_t, color_t, f};

use super::*;

use crate::test_support::{TestReps, test_initialize};

fn assert_color_zero(expr: Atom) {
    assert!(
        crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .is_zero()
    );
}

#[test]
fn test_color_structures() {
    test_initialize();
    let f = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            ColorAdjoint {}.new_rep(8),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(4),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(7),
        ],
        symbol!("test"),
        None,
    );
    let logical_indices: [AbstractIndex; 6] =
        [5.into(), 4.into(), 2.into(), 3.into(), 1.into(), 0.into()];
    let storage_indices = f.layout().logical_to_canonical(&logical_indices);
    let input_layout = f.layout().clone();
    let order = f.canonical().order();
    let f = f
        .into_canonical()
        .reindex_storage(&storage_indices)
        .unwrap()
        .map_target(|a| SymbolicTensor::from_named(&a).unwrap());

    let f_p = f.apply();

    // println!("{}", f_p);
    let simplified = f_p.expression.schoonschip();
    // println!("{}", simplified);
    let f_parsed = ShadowedStructure::<AbstractIndex>::parse(simplified.as_view()).unwrap();

    assert_eq!(order, f_parsed.canonical().order());
    let logical_positions = (0..order).collect::<Vec<_>>();
    let canonical_positions = input_layout.logical_to_canonical(&logical_positions);
    assert_eq!(
        input_layout.canonical_to_logical(&canonical_positions),
        logical_positions
    );
}

#[test]
fn test_color_simplification() {
    test_initialize();

    let atom = parse_lit!(
        f(
            coad(Nc ^ 2 - 1, 1),
            coad(Nc ^ 2 - 1, 2),
            coad(Nc ^ 2 - 1, 3)
        ) ^ 2,
        default_namespace = "spenso"
    );
    let dimension_value = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let dimension = Atom::var(symbol!("test_color_simplification_dimension"));
    let admitted = atom
        .replace(dimension_value.to_pattern())
        .with(dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(dimension.to_pattern())
            .with(dimension_value.to_pattern()),
        atom
    );
    let simplified = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .replace(dimension.to_pattern())
        .with(dimension_value.to_pattern());

    assert_snapshot!(simplified.to_bare_ordered_string(), @"(-1+Nc^2)*cas(2,coad(-1+Nc^2))");
}

#[test]
fn concrete_qcd_representations_can_be_made_parametric() {
    test_initialize();
    let nc = Atom::var(CS.nc);
    let adjoint_index = Atom::var(s!(a));
    let fundamental_index = Atom::var(s!(i));
    let concrete = ColorAdjoint {}.to_symbolic([Atom::num(8), adjoint_index.clone()])
        * ColorFundamental {}.to_symbolic([Atom::num(3), fundamental_index.clone()]);
    let expected = ColorAdjoint {}
        .to_symbolic([nc.clone().pow(Atom::num(2)) - Atom::one(), adjoint_index])
        * ColorFundamental {}.to_symbolic([nc, fundamental_index]);

    assert_eq!(concrete.to_parametric_color(), expected);
}

#[test]
fn untyped_structure_constants_do_not_assume_an_adjoint_dimension() {
    test_initialize();

    let atom = f!(3, 1, 5) * f!(3, 5, 1);
    assert!(SymbolicTensor::infer(atom.clone()).is_err());
    assert_eq!(
        super::simplify::ColorAlgebraSimplifier {
            settings: ColorSimplifySettings::default(),
        }
        .step(atom.as_view(), true),
        atom
    );
}

#[test]
fn partially_typed_or_inconsistent_structure_constants_are_not_rewritten() {
    test_initialize();
    let partially_typed = parse_lit!(f(coad(dA, a), b, c) ^ 2, default_namespace = "spenso");
    assert!(SymbolicTensor::infer(partially_typed.clone()).is_err());
    assert_eq!(
        super::simplify::ColorAlgebraSimplifier {
            settings: ColorSimplifySettings::default(),
        }
        .step(partially_typed.as_view(), true),
        partially_typed
    );

    let inconsistent = parse_lit!(
        f(coad(dA, a), coad(dA, b), coad(other_dA, c)) ^ 2,
        default_namespace = "spenso"
    );
    assert!(SymbolicTensor::infer(inconsistent.clone()).is_err());
    assert_eq!(
        super::simplify::ColorAlgebraSimplifier {
            settings: ColorSimplifySettings::default(),
        }
        .step(inconsistent.as_view(), true),
        inconsistent
    );
}

#[test]
fn color_generator_macro_accepts_gamma_style_index_shorthands() {
    test_initialize();
    let r = TestReps::new();
    let cof_nc = ColorFundamental {}.new_rep(CS.nc);

    let a = slot!(r.coad_da, a);
    assert_eq!(color_t!(a), color_t!(slot!(r.coad_da, a)));
    assert_snapshot!(color_t!(1).to_bare_ordered_string(), @"t(1,in,out)");

    let a = slot!(r.coad_da, a);
    let i = slot!(cof_nc, i);
    let j = slot!(cof_nc.dual(), j);
    assert_eq!(
        color_t!(a, i, j),
        color_t!(
            slot!(r.coad_da, a),
            slot!(cof_nc, i),
            slot!(cof_nc.dual(), j)
        ),
    );
    assert_eq!(
        color_t!(RS.a__, RS.i__, RS.j__),
        CS.explicit_t(
            ColorAdjoint {}.to_symbolic([RS.a__]),
            ColorFundamental {}.to_symbolic([RS.i__]),
            ColorAntiFundamental {}.to_symbolic([RS.j__]),
        ),
    );
    assert_eq!(
        color_t!([CS.adj_, RS.a_], [CS.nc_, RS.i_], [CS.nc_, RS.j_]),
        CS.t_pattern(CS.nc_, CS.adj_, RS.a_, RS.i_, RS.j_),
    );
}

#[test]
fn color_structure_macro_accepts_gamma_style_index_shorthands() {
    test_initialize();
    let r = TestReps::new();

    let a = slot!(r.coad_da, a);
    let b = slot!(r.coad_da, b);
    let c = slot!(r.coad_da, c);
    assert_eq!(
        color_f!(a, b, c),
        color_f!(
            slot!(r.coad_da, a),
            slot!(r.coad_da, b),
            slot!(r.coad_da, c)
        ),
    );
    assert_eq!(
        color_f!(RS.a__, RS.b__, RS.c__),
        CS.structure_f(
            ColorAdjoint {}.to_symbolic([RS.a__]),
            ColorAdjoint {}.to_symbolic([RS.b__]),
            ColorAdjoint {}.to_symbolic([RS.c__]),
        ),
    );
    assert_eq!(
        color_f!([CS.adj_, RS.a_], [CS.adj_, RS.b_], [CS.adj_, RS.c_]),
        CS.f_pattern(CS.adj_, RS.a_, RS.b_, RS.c_),
    );
}

#[test]
fn color_invariant_macros_build_scalar_heads() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let default_adjoint_dimension = n.clone().pow(Atom::num(2)) - Atom::num(1);
    let cof_n = ColorFundamental {}.to_symbolic([n.clone()]);
    let a = ColorAdjoint {}.to_symbolic([default_adjoint_dimension.clone(), Atom::var(s!(a))]);
    let b = ColorAdjoint {}.to_symbolic([default_adjoint_dimension, Atom::var(s!(b))]);

    assert_eq!(
        color_cas!(2, cof_n.clone()),
        CS.cas(Atom::num(2), cof_n.clone())
    );
    assert_eq!(
        color_idx!(2, cof_n.clone()),
        CS.idx(Atom::num(2), cof_n.clone())
    );
    assert_eq!(
        color_gram!(3, cof_n.clone(), cof_n.clone()),
        CS.gram(Atom::num(3), cof_n.clone(), cof_n.clone())
    );
    assert_eq!(
        color_str_t!(cof_n.clone(), a.clone(), b.clone()),
        trace!(
            cof_n.clone(),
            sym!(color_t!(a.clone()), color_t!(b.clone()))
        )
    );

    let expr = color_cas!(2, cof_n.clone())
        * color_idx!(2, cof_n.clone())
        * color_gram!(3, cof_n.clone(), cof_n);
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert!(net.graph.dangling_indices().is_empty());
    assert_eq!(net.simple_execute::<()>().unwrap(), expr);
}

#[test]
fn color_invariant_print_special_cases_are_compact() {
    test_initialize();
    let cof_n = ColorFundamental {}.to_symbolic([Atom::var(CS.nc)]);

    let mut compact = SpensoPrintSettings::compact().nice_symbolica();
    compact.color_builtin_symbols = false;
    assert_eq!(
        color_cas!(2, cof_n.clone()).printer(compact).to_string(),
        "CF"
    );

    let typst = SpensoPrintSettings::typst_options();
    assert_eq!(color_cas!(2, cof_n).printer(typst).to_string(), "C_F");
}

#[test]
fn cof_dimension_invariant_rules_substitute_supported_fundamental_cases() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let default_adjoint_dimension = n.clone().pow(Atom::num(2)) - Atom::num(1);
    let cof_n = ColorFundamental {}.to_symbolic([n.clone()]);
    let coad_n = ColorAdjoint {}.to_symbolic([default_adjoint_dimension.clone()]);

    assert_eq!(
        color_idx!(2, cof_n.clone()).to_cof_dimension_invariants(),
        Atom::num(1) / Atom::num(2)
    );
    assert_eq!(
        color_cas!(2, cof_n.clone())
            .to_cof_dimension_invariants()
            .expand(),
        ((n.clone().pow(Atom::num(2)) - Atom::num(1)) / (Atom::num(2) * n.clone())).expand()
    );
    assert_eq!(
        color_cas!(2, coad_n).to_cof_dimension_invariants(),
        n.clone()
    );
    assert_eq!(
        color_cas!(2, ColorAdjoint {}.to_symbolic([Atom::num(8)])).to_cof_dimension_invariants(),
        Atom::num(3)
    );
    assert_eq!(
        (Atom::var(CS.ca) * Atom::var(CS.tr) * default_adjoint_dimension)
            .to_cof_dimension_invariants()
            .expand(),
        (n.clone() * (Atom::num(1) / Atom::num(2)) * (n.clone().pow(Atom::num(2)) - Atom::num(1)))
            .expand()
    );

    let gram_three = color_gram!(3, cof_n.clone(), cof_n.clone())
        .to_cof_dimension_invariants()
        .expand();
    let expected_three = ((n.clone().pow(Atom::num(2)) - Atom::num(1))
        * (n.clone().pow(Atom::num(2)) - Atom::num(4))
        / (Atom::num(16) * n.clone()))
    .expand();
    assert_eq!(gram_three, expected_three);

    let gram_four = color_gram!(4, cof_n.clone(), cof_n)
        .to_cof_dimension_invariants()
        .expand();
    let n_squared = n.clone().pow(Atom::num(2));
    let expected_four = ((n_squared.clone() - Atom::num(1))
        * (n_squared.clone().pow(Atom::num(2)) - Atom::num(6) * n_squared + Atom::num(18))
        / (Atom::num(96) * n.pow(Atom::num(2))))
    .expand();
    assert_eq!(gram_four, expected_four);
}

#[test]
fn cof_dimension_simplification_preserves_factorized_spectators() {
    test_initialize();
    let invariants = parse_lit!(
        cas(2, coad(8)) * idx(2, cof(3)),
        default_namespace = "spenso"
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((invariants).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        invariants
    );

    let spectator = parse_lit!((opaque(x) + opaque(y)) ^ 3 * (opaque(z) + opaque(w)));
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    let simplified = crate::tensor::SymbolicTensor::infer(
        (invariants * spectator.clone()).as_atom_view().to_owned(),
    )
    .unwrap()
    .simplify_color(settings)
    .unwrap()
    .resolved()
    .unwrap()
    .into_expression();
    // Exact Atom equality checks that the spectator's sums and power stay intact.
    assert_eq!(simplified, (Atom::num(3) / Atom::num(2)) * spectator);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((simplified).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        simplified
    );
}

#[test]
fn color_collection_preserves_powered_noncolor_trace_sums() {
    test_initialize();
    // Opaque word factors still declare their structural channels.
    let _ = spenso::tensor_symbol!("spenso::opaque_trace_a");
    let _ = spenso::tensor_symbol!("spenso::opaque_trace_b");
    let traces = parse!(
        "(trace(bis(4), opaque_trace_a(in, out)) + trace(bis(4), opaque_trace_b(in, out))) ^ 5",
        default_namespace = "spenso"
    );
    let color = parse_lit!(
        f(coad(8, a), coad(8, b), coad(8, c)) ^ 2,
        default_namespace = "spenso"
    );
    for settings in [
        ColorSimplifySettings::default().with_cof_dimension_invariants(),
        ColorSimplifySettings::default()
            .with_cof_dimension_invariants()
            .without_trace_evaluation(),
    ] {
        // Non-color traces are opaque numerator factors. Collecting them as
        // polynomial variables would distribute this power of a sum.
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((traces).as_atom_view().to_owned())
                .unwrap()
                .simplify_color(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            traces
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((&traces * &color).as_atom_view().to_owned())
                .unwrap()
                .simplify_color(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            Atom::num(24) * &traces
        );
    }
}

#[test]
fn color_collection_closes_generic_chains_after_metric_reduction() {
    test_initialize();
    // Opaque word factors still declare their structural channels.
    let _ = spenso::tensor_symbol!("spenso::opaque_matrix");
    let chains = parse!(
        "spenso::chain(spenso::mink(4,spenso::i),spenso::mink(4,spenso::j),spenso::g(spenso::in,spenso::out))*spenso::chain(spenso::mink(4,spenso::j),spenso::mink(4,spenso::i),spenso::opaque_matrix(spenso::in,spenso::out))"
    );
    let trace = parse!(
        "trace(mink(4), cyclic(opaque_matrix(in, out)))",
        default_namespace = "spenso"
    );
    let spectator = parse_lit!((opaque(x) + opaque(y)) ^ 3);
    // The identity chain first becomes a metric. Contracting that metric then
    // closes the other chain, so one pass of independent node rewrites is not enough.
    let expected = &spectator * trace;
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((spectator * chains).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expected
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((expected).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expected
    );
}

#[test]
fn color_collection_preserves_generic_chain_normalization() {
    test_initialize();
    // Opaque word factors still declare their structural channels.
    let _ = spenso::tensor_symbol!("spenso::spectator_tensor");
    let spectator = parse_lit!((opaque(x) + opaque(y)) ^ 5);
    let chain = parse!(
        "chain(mink(4, i), mink(4, i), spectator_tensor(in, out))",
        default_namespace = "spenso"
    );
    let expected = parse!(
        "trace(mink(4), cyclic(spectator_tensor(in, out)))",
        default_namespace = "spenso"
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((spectator.clone() * chain).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        spectator * expected
    );
    let identity = parse!(
        "spenso::chain(spenso::mink(4,spenso::i),spenso::mink(4,spenso::j),spenso::g(spenso::in,spenso::out))"
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((identity).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        parse_lit!(g(mink(4, i), mink(4, j)), default_namespace = "spenso")
    );
}

#[test]
fn color_collection_preserves_wrappers_and_mixed_tensor_slots() {
    test_initialize();
    // Opaque word factors still declare their structural channels.
    let _ = spenso::tensor_symbol!("spenso::mixed_tensor");
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    let color = parse_lit!(
        f(coad(8, a), coad(8, b), coad(8, c)) ^ 2,
        default_namespace = "spenso"
    );
    let spectator = parse_lit!((opaque(x) + opaque(y)) ^ 3);
    for wrapper in [
        SPENSO_TAG.pure_scalar,
        SPENSO_TAG.bracket,
        symbol!("color_wrapper"),
    ] {
        let input = function!(wrapper, &color * &spectator);
        let scalar = Atom::num(24) * &spectator;
        let expected = if wrapper == SPENSO_TAG.bracket {
            scalar
        } else {
            function!(wrapper, scalar)
        };
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                .unwrap()
                .simplify_color(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            expected
        );
    }
    let mixed = parse_lit!(
        mixed_tensor(cof(3, i), mink(4, mu)) * g(mink(4, mu), mink(4, nu)),
        default_namespace = "spenso"
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((mixed).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        parse_lit!(
            mixed_tensor(cof(3, i), mink(4, nu)),
            default_namespace = "spenso"
        )
    );
}

#[test]
fn color_payload_reduction_preserves_non_color_coefficients() {
    test_initialize();
    let _ = spenso::vector_symbol!("spenso::p");
    // Opaque word factors still declare their structural channels.
    let _ = spenso::tensor_symbol!("spenso::spectator_tensor");
    let _ = spenso::tensor_symbol!("spenso::mixed_tensor");
    let settings = ColorSimplifySettings {
        simplify_non_color: false,
        ..ColorSimplifySettings::default().with_cof_dimension_invariants()
    };
    let color = parse_lit!(
        f(coad(8, a), coad(8, b), coad(8, c)) ^ 2,
        default_namespace = "spenso"
    );
    let spectator = parse!(
        "(opaque(x) + opaque(y)) ^ 5 * g(mink(4, mu), mink(4, nu)) * p(mink(4, mu)) * trace(mink(4), cyclic(spectator_tensor(in, out)))",
        default_namespace = "spenso"
    );
    let input = &color * &spectator;
    let simplified = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
        .unwrap()
        .simplify_color(settings)
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    assert_eq!(simplified, Atom::num(24) * spectator);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((simplified).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        simplified
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((simplified).as_atom_view().to_owned())
            .unwrap()
            .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(ColorSimplifySettings::default().with_cof_dimension_invariants())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    );

    let mixed = parse_lit!(
        mixed_tensor(cof(3, i), mink(4, mu)) * g(mink(4, mu), mink(4, nu)),
        default_namespace = "spenso"
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((mixed).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        mixed
    );

    let fundamental = ColorFundamental {}.new_rep(3);
    let adjoint = ColorAdjoint {}.new_rep(8);
    let open = chain!(
        slot!(fundamental, i),
        slot!(fundamental.dual(), j),
        color_t!(slot!(adjoint, a)),
        color_t!(slot!(adjoint, a)),
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((open).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        parse_lit!(
            4 / 3 * g(cof(3, idenso::i), dind(cof(3, idenso::j))),
            default_namespace = "spenso"
        )
    );
}

#[test]
fn color_trace_projectors_preserve_adjoint_slots_and_algebra() {
    use spenso::network::parsing::{AtomStructureExt, NetworkParse};
    use spenso::structure::{OrderedStructure, TensorStructure, representation::LibraryRep};

    test_initialize();
    let fundamental = ColorFundamental {}.new_rep(3);
    let adjoint = ColorAdjoint {}.new_rep(8);
    let factors = [
        color_t!(slot!(adjoint, a)),
        color_t!(slot!(adjoint, b)),
        color_t!(slot!(adjoint, c)),
    ];
    let ordered = trace!(&fundamental; factors.clone());
    let symmetric = trace!(&fundamental, spenso::shadowing::sym(factors.clone()));
    let antisymmetric = trace!(&fundamental, spenso::shadowing::antisym(factors));
    let closing = color_f!(slot!(adjoint, a), slot!(adjoint, b), slot!(adjoint, c));
    let settings = ColorSimplifySettings {
        simplify_non_color: false,
        ..ColorSimplifySettings::default().with_cof_dimension_invariants()
    };
    for (expression, expected) in [
        (&ordered, Atom::i() * 6),
        (&symmetric, Atom::Zero),
        (&antisymmetric, Atom::i() * 6),
    ] {
        let fast = expression
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            .unwrap();
        let expanded = expression
            .parse_to_atom_net::<AbstractIndex>(&spenso::network::parsing::ParseSettings::default())
            .unwrap();
        let expanded = OrderedStructure::new(expanded.graph.dangling_indices());
        assert_eq!(fast.canonical().order(), 3);
        assert_eq!(fast.canonical(), expanded.canonical());
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((expression * &closing).as_atom_view().to_owned())
                .unwrap()
                .simplify_color(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            expected
        );
    }
    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            (&symmetric + &antisymmetric).as_atom_view().to_owned()
        )
        .unwrap()
        .simplify_color(settings)
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
        crate::tensor::SymbolicTensor::infer((ordered).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    );
}

#[test]
fn cof_dimension_simplification_resolves_new_color_invariants() {
    test_initialize();
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    let closed = parse_lit!(
        f(coad(8, a), coad(8, b), coad(8, c)) ^ 2,
        default_namespace = "spenso"
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((closed).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        Atom::num(24)
    );

    let open = parse_lit!(
        f(coad(8, a), coad(8, b), coad(8, c)) * f(coad(8, a), coad(8, b), coad(8, d)),
        default_namespace = "spenso"
    );
    let expected = parse_lit!(3 * g(coad(8, c), coad(8, d)), default_namespace = "spenso");
    let simplified = crate::tensor::SymbolicTensor::infer((open).as_atom_view().to_owned())
        .unwrap()
        .simplify_color(settings)
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    assert_eq!(simplified, expected);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((simplified).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        simplified
    );

    let fundamental = ColorFundamental {}.new_rep(3);
    let adjoint = ColorAdjoint {}.new_rep(8);
    let trace = trace!(
        &fundamental,
        color_t!(slot!(adjoint, a)),
        color_t!(slot!(adjoint, a)),
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((trace).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        Atom::num(4)
    );
}

#[test]
fn cof_dimension_simplification_preserves_unsupported_invariants() {
    test_initialize();
    let expression = parse_lit!(
        cas(3, cof(3)) * idx(3, cof(3)) * cas(2, coad(7)),
        default_namespace = "spenso"
    );
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((expression).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expression
    );
}

#[test]
fn cof_dimension_simplification_respects_disabled_trace_evaluation() {
    test_initialize();
    let fundamental = ColorFundamental {}.new_rep(3);
    let adjoint = ColorAdjoint {}.new_rep(8);
    let closed_chain = chain!(
        slot!(fundamental, i),
        slot!(fundamental.dual(), i),
        color_t!(slot!(adjoint, a)),
        color_t!(slot!(adjoint, b)),
    );
    let invariant = color_cas!(2, ColorAdjoint {}.to_symbolic([Atom::num(8)]));
    let expected = Atom::num(3)
        * trace!(
            &fundamental,
            color_t!(slot!(adjoint, a)),
            color_t!(slot!(adjoint, b)),
        );
    let settings = ColorSimplifySettings::default()
        .without_trace_evaluation()
        .with_cof_dimension_invariants();
    let simplified =
        crate::tensor::SymbolicTensor::infer((invariant * closed_chain).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
    assert_eq!(simplified, expected);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((simplified).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        simplified
    );
}

#[test]
fn cof_dimension_simplification_respects_disabled_fierz_expansion() {
    test_initialize();
    let fundamental = ColorFundamental {}.new_rep(3);
    let adjoint = ColorAdjoint {}.new_rep(8);
    let chains = chain!(
        slot!(fundamental, i),
        slot!(fundamental.dual(), j),
        color_t!(slot!(adjoint, a)),
    ) * chain!(
        slot!(fundamental, k),
        slot!(fundamental.dual(), l),
        color_t!(slot!(adjoint, a)),
    );
    let invariant = color_idx!(2, ColorFundamental {}.to_symbolic([Atom::num(3)]));
    let settings = ColorSimplifySettings::default()
        .without_cross_chain_fierz_expansion()
        .with_cof_dimension_invariants();
    let simplified = crate::tensor::SymbolicTensor::infer(
        (invariant * chains.clone()).as_atom_view().to_owned(),
    )
    .unwrap()
    .simplify_color(settings)
    .unwrap()
    .resolved()
    .unwrap()
    .into_expression();
    assert_eq!(simplified, chains / Atom::num(2));
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((simplified).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        simplified
    );
}

#[test]
fn color_structure_symbol_is_antisymmetric() {
    test_initialize();
    let r = TestReps::new();

    let a = slot!(r.coad_da, a);
    let b = slot!(r.coad_da, b);
    let c = slot!(r.coad_da, c);

    assert_snapshot!(
        color_f!(a, c, b).to_bare_ordered_string(),
        @"-1*f(coad(dA,a),coad(dA,b),coad(dA,c))"
    );
    assert!(color_f!(a, b, b).is_zero());
}

#[test]
fn permuted_structure_constant_square_simplifies_with_sign() {
    test_initialize();
    let atom = parse_lit!(
        f(coad(8, a), coad(8, b), coad(8, c)) * f(coad(8, a), coad(8, c), coad(8, b)),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((atom).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"-8*cas(2,coad(8))");
}

#[test]
fn three_loop_pole_part_color() {
    test_initialize();
    let input = parse_lit!(
        ((-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 56 / 3 * eps)
            * f(
                coad(8, hedge(11)),
                coad(8, hedge(13)),
                coad(8, vertex(2, 1))
            )
            * f(coad(8, hedge(11)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
            * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            * f(coad(8, hedge(4)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
            + (-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 56 / 3 * eps)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(3, 1))
                )
                * f(coad(8, hedge(11)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(6)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 88 / 3 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
            + (-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 88 / 3 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-16 / 3 * eps ^ 2 * 𝜋 ^ 2 + -32 + -52 * eps ^ 2 + 176 / 3 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
            + (-16 / 3 * eps ^ 2 * 𝜋 ^ 2 + -32 + -52 * eps ^ 2 + 48 * eps)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(2, 1))
                )
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(3, 1))
                )
                * f(coad(8, hedge(4)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(6)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-16 / 3 * eps ^ 2 * 𝜋 ^ 2 + -32 + -52 * eps ^ 2 + 48 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-88 / 3 * eps + 16 + 26 * eps ^ 2 + 8 / 3 * eps ^ 2 * 𝜋 ^ 2)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(2, 1))
                )
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(4)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
            + (-88 / 3 * eps + 16 + 26 * eps ^ 2 + 8 / 3 * eps ^ 2 * 𝜋 ^ 2)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(3, 1))
                )
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(6)), coad(8, hedge(9)), coad(8, vertex(3, 1))))
            * 1
            / 64
            * CA
            * dot(P(0, mink(4)), P(0, mink(4)))
            * g(coad(8, hedge(4)), coad(8, hedge(6)))
            * gs
            ^ 6 * eps
            ^ (-3),
        default_namespace = "spenso"
    );

    let color_zero_candidate =
        crate::tensor::SymbolicTensor::infer((input.cook_indices()).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .collect_reps([
                ColorAdjoint {}.into(),
                ColorFundamental {}.into(),
                ColorSextet {}.into(),
            ]);

    let common = parse_lit!(
        CA * eps ^ (-3) * gs ^ 6 * cas(2, coad(8)) ^ 2 * dot(P(0, mink(4)), P(0, mink(4))),
        default_namespace = "spenso"
    );
    let expected_coefficient = parse_lit!(
        -18 - 117 / 4 * eps ^ 2 - 3 * eps ^ 2 * 𝜋 ^ 2 + 29 * eps,
        default_namespace = "spenso"
    );
    let coefficient = color_zero_candidate.collect_symbol::<i16>(SPENSO_TAG.dot) / common;
    assert!((coefficient - expected_coefficient).expand().is_zero());

    let input = parse_lit!(
        ((8 * eps + 8 / 3) * 1 / 64
            * f(coad(8, hedge(1)), coad(8, hedge(11)), coad(8, hedge(15)))
            * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(5)))
            * f(coad(8, hedge(11)), coad(8, hedge(13)), coad(8, hedge(9)))
            * f(coad(8, hedge(13)), coad(8, hedge(5)), coad(8, hedge(7)))
            * f(coad(8, hedge(15)), coad(8, hedge(17)), coad(8, hedge(7)))
            + -1 / 16
                * f(coad(8, hedge(1)), coad(8, hedge(11)), coad(8, hedge(15)))
                * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(4)))
                * f(coad(8, hedge(11)), coad(8, hedge(13)), coad(8, hedge(9)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, hedge(7)))
                * f(coad(8, hedge(15)), coad(8, hedge(17)), coad(8, hedge(7)))
            + 1 / 16
                * f(coad(8, hedge(1)), coad(8, hedge(10)), coad(8, hedge(14)))
                * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(5)))
                * f(coad(8, hedge(10)), coad(8, hedge(13)), coad(8, hedge(9)))
                * f(coad(8, hedge(13)), coad(8, hedge(5)), coad(8, hedge(7)))
                * f(coad(8, hedge(14)), coad(8, hedge(17)), coad(8, hedge(7))))
            * dot(P(0, mink(4)), P(0, mink(4)))
            * f(coad(8, hedge(17)), coad(8, hedge(3)), coad(8, hedge(9)))
            * gs
            ^ 6 * eps
            ^ (-2),
        default_namespace = "spenso"
    );

    let color_zero_candidate =
        crate::tensor::SymbolicTensor::infer((input.cook_indices()).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .collect_reps([
                ColorAdjoint {}.into(),
                ColorFundamental {}.into(),
                ColorSextet {}.into(),
            ]);

    assert_snapshot!(
        &color_zero_candidate
            .collect_symbol::<i16>(SPENSO_TAG.dot)
            .to_bare_ordered_string(),
        @"0"
    )
}

#[test]
fn kaapo_gl34_color_input_simplifies_to_zero() {
    test_initialize();

    // Captured from feyngen GL34 after `to_param_color()` and before `simplify_color()`.
    let input = parse!(
        "
        -1𝑖 * UFO::G^6
        * (
            -Q(2, mink(4, hedge(11))) * g(mink(4, hedge(5)), mink(4, hedge(17)))
            + Q(2, mink(4, hedge(17))) * g(mink(4, hedge(5)), mink(4, hedge(11)))
            + Q(5, mink(4, hedge(5))) * g(mink(4, hedge(11)), mink(4, hedge(17)))
            - Q(5, mink(4, hedge(17))) * g(mink(4, hedge(5)), mink(4, hedge(11)))
            - Q(8, mink(4, hedge(5))) * g(mink(4, hedge(11)), mink(4, hedge(17)))
            + Q(8, mink(4, hedge(11))) * g(mink(4, hedge(5)), mink(4, hedge(17)))
        )
        * (
            -Q(6, mink(4, hedge(15))) * g(mink(4, hedge(13)), mink(4, hedge(16)))
            + Q(6, mink(4, hedge(16))) * g(mink(4, hedge(13)), mink(4, hedge(15)))
            + Q(7, mink(4, hedge(13))) * g(mink(4, hedge(15)), mink(4, hedge(16)))
            - Q(7, mink(4, hedge(16))) * g(mink(4, hedge(13)), mink(4, hedge(15)))
            + Q(8, mink(4, hedge(13))) * g(mink(4, hedge(15)), mink(4, hedge(16)))
            - Q(8, mink(4, hedge(15))) * g(mink(4, hedge(13)), mink(4, hedge(16)))
        )
        * Q(0, mink(4, edge(0, 1)))
        * Q(1, mink(4, edge(1, 1)))
        * Q(3, mink(4, edge(3, 1)))
        * Q(4, mink(4, edge(4, 1)))
        * g(mink(4, hedge(4)), mink(4, hedge(5)))
        * g(mink(4, hedge(10)), mink(4, hedge(11)))
        * g(mink(4, hedge(12)), mink(4, hedge(13)))
        * g(mink(4, hedge(14)), mink(4, hedge(15)))
        * g(mink(4, hedge(16)), mink(4, hedge(17)))
        * g(dind(cof(Nc, hedge(1))), cof(Nc, hedge(0)))
        * g(dind(cof(Nc, hedge(2))), cof(Nc, hedge(3)))
        * g(dind(cof(Nc, hedge(6))), cof(Nc, hedge(7)))
        * g(dind(cof(Nc, hedge(9))), cof(Nc, hedge(8)))
        * g(coad(-1 + Nc^2, hedge(4)), coad(-1 + Nc^2, hedge(5)))
        * g(coad(-1 + Nc^2, hedge(10)), coad(-1 + Nc^2, hedge(11)))
        * g(coad(-1 + Nc^2, hedge(12)), coad(-1 + Nc^2, hedge(13)))
        * g(coad(-1 + Nc^2, hedge(14)), coad(-1 + Nc^2, hedge(15)))
        * g(coad(-1 + Nc^2, hedge(16)), coad(-1 + Nc^2, hedge(17)))
        * spenso::gamma(bis(4, hedge(0)), bis(4, hedge(2)), mink(4, hedge(4)))
        * spenso::gamma(bis(4, hedge(3)), bis(4, hedge(9)), mink(4, hedge(14)))
        * spenso::gamma(bis(4, hedge(7)), bis(4, hedge(1)), mink(4, hedge(12)))
        * spenso::gamma(bis(4, hedge(8)), bis(4, hedge(6)), mink(4, hedge(10)))
        * spenso::gamma(bis(4, hedge(1)), bis(4, hedge(0)), mink(4, edge(0, 1)))
        * spenso::gamma(bis(4, hedge(2)), bis(4, hedge(3)), mink(4, edge(1, 1)))
        * spenso::gamma(bis(4, hedge(6)), bis(4, hedge(7)), mink(4, edge(3, 1)))
        * spenso::gamma(bis(4, hedge(9)), bis(4, hedge(8)), mink(4, edge(4, 1)))
        * t(coad(-1 + Nc^2, hedge(4)), cof(Nc, hedge(2)), dind(cof(Nc, hedge(0))))
        * t(coad(-1 + Nc^2, hedge(10)), cof(Nc, hedge(6)), dind(cof(Nc, hedge(8))))
        * t(coad(-1 + Nc^2, hedge(12)), cof(Nc, hedge(1)), dind(cof(Nc, hedge(7))))
        * t(coad(-1 + Nc^2, hedge(14)), cof(Nc, hedge(9)), dind(cof(Nc, hedge(3))))
        * f(coad(-1 + Nc^2, hedge(5)), coad(-1 + Nc^2, hedge(11)), coad(-1 + Nc^2, hedge(17)))
        * f(coad(-1 + Nc^2, hedge(13)), coad(-1 + Nc^2, hedge(15)), coad(-1 + Nc^2, hedge(16)))
        ",
        default_namespace = "spenso"
    );
    let dimension_value = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let dimension = Atom::var(symbol!("kaapo_gl34_adjoint_dimension"));
    let indexed = input
        .replace(dimension_value.to_pattern())
        .with(dimension.to_pattern());
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    let admitted = cooking.try_cook_indices(indexed.as_view()).unwrap();
    assert_eq!(
        cooking
            .uncook(admitted.as_view())
            .replace(dimension.to_pattern())
            .with(dimension_value.to_pattern()),
        input
    );
    let color_zero_candidate = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .collect_reps([
            ColorAdjoint {}.into(),
            ColorFundamental {}.into(),
            ColorSextet {}.into(),
        ]);

    assert!(
        color_zero_candidate.is_zero(),
        "Expected color to be zero, got {}",
        color_zero_candidate.printer(SpensoPrintSettings::compact().nice_symbolica())
    );
}

#[test]
fn color_casimir_basis_rewrites_dimensions() {
    test_initialize();
    let fundamental_dimension = Atom::var(s!(dF));
    let adjoint_dimension = Atom::var(s!(dA));
    let fundamental_rep = ColorFundamental {}.to_symbolic([fundamental_dimension.clone()]);
    let adjoint_rep = ColorAdjoint {}.to_symbolic([adjoint_dimension.clone()]);
    let fundamental_casimir = color_cas!(2, fundamental_rep.clone());
    let adjoint_casimir = color_cas!(2, adjoint_rep.clone());
    let fundamental_index = color_idx!(2, fundamental_rep.clone());

    let rewritten = (fundamental_dimension + adjoint_dimension)
        .to_color_casimir(fundamental_rep.as_view(), adjoint_rep.as_view());
    let expected =
        adjoint_casimir.clone() + adjoint_casimir * fundamental_casimir / fundamental_index;
    assert_eq!(rewritten.expand(), expected.expand());
    assert_eq!(
        rewritten.to_color_casimir(fundamental_rep.as_view(), adjoint_rep.as_view()),
        rewritten
    );
}

#[test]
fn color_casimir_basis_policies_control_su_n_conventions() {
    test_initialize();
    let fundamental_dimension = Atom::var(s!(dF));
    let adjoint_dimension = Atom::var(s!(dA));
    let fundamental_rep = ColorFundamental {}.to_symbolic([fundamental_dimension.clone()]);
    let adjoint_rep = ColorAdjoint {}.to_symbolic([adjoint_dimension.clone()]);
    let fundamental_casimir = color_cas!(2, fundamental_rep.clone());
    let adjoint_casimir = color_cas!(2, adjoint_rep.clone());
    let fundamental_index = color_idx!(2, fundamental_rep.clone());

    let general = adjoint_dimension.clone().to_color_casimir_with(
        fundamental_rep.as_view(),
        adjoint_rep.as_view(),
        ColorCasimirSettings::default().without_fundamental_dimension_rewrite(),
    );
    assert_eq!(
        general.expand(),
        (fundamental_dimension * fundamental_casimir.clone() / fundamental_index).expand()
    );

    let normalized = adjoint_dimension.to_color_casimir_with(
        fundamental_rep.as_view(),
        adjoint_rep.as_view(),
        ColorCasimirSettings::default().with_fundamental_index_normalization(),
    );
    assert_eq!(
        normalized.expand(),
        (Atom::num(2) * adjoint_casimir * fundamental_casimir).expand()
    );
}

#[test]
fn color_casimir_basis_preserves_representation_subtrees_and_concrete_dimensions() {
    test_initialize();
    let fundamental_dimension = Atom::var(s!(dF));
    let adjoint_dimension = Atom::var(s!(dA));
    let fundamental_rep = ColorFundamental {}.to_symbolic([fundamental_dimension.clone()]);
    let adjoint_rep = ColorAdjoint {}.to_symbolic([adjoint_dimension.clone()]);
    let probe = symbol!("probe");
    let expression = function!(
        probe,
        fundamental_rep.clone(),
        adjoint_rep.clone(),
        fundamental_dimension.clone(),
        adjoint_dimension.clone()
    );
    let expected = expression.clone();
    assert_eq!(
        expression.to_color_casimir(fundamental_rep.as_view(), adjoint_rep.as_view()),
        expected
    );

    let tensor = function!(
        SPENSO_TAG.tensor_symbol("casimir_metadata"),
        adjoint_dimension.clone(),
        adjoint_rep.clone()
    );
    let rewritten = (adjoint_dimension * tensor.clone())
        .to_color_casimir(fundamental_rep.as_view(), adjoint_rep.as_view());
    let expected = color_cas!(2, adjoint_rep.clone()) * color_cas!(2, fundamental_rep.clone())
        / color_idx!(2, fundamental_rep.clone())
        * tensor;
    assert_eq!(rewritten.expand(), expected.expand());

    let fundamental_rep = ColorFundamental {}.to_symbolic([Atom::num(3)]);
    let adjoint_rep = ColorAdjoint {}.to_symbolic([Atom::num(8)]);
    let numeric = (Atom::num(3) + Atom::num(8))
        .to_color_casimir(fundamental_rep.as_view(), adjoint_rep.as_view());
    assert_eq!(numeric, Atom::num(11));

    let adjoint_dimension = Atom::var(s!(dA));
    let adjoint_rep = ColorAdjoint {}.to_symbolic([adjoint_dimension.clone()]);
    let rewritten =
        adjoint_dimension.to_color_casimir(fundamental_rep.as_view(), adjoint_rep.as_view());
    let expected =
        Atom::num(3) * color_cas!(2, fundamental_rep.clone()) / color_idx!(2, fundamental_rep);
    assert_eq!(rewritten.expand(), expected.expand());
}

#[test]
fn antisymmetric_two_generator_trace_vanishes() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(color_t!(slot!(r.coad_da, a)), color_t!(slot!(r.coad_da, b)))
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"0");
}

#[test]
fn antisymmetric_three_generator_trace_reduces_to_structure_constant() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(
            color_t!(slot!(r.coad_da, a)),
            color_t!(slot!(r.coad_da, b)),
            color_t!(slot!(r.coad_da, c)),
        )
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"f(coad(dA,a),coad(dA,b),coad(dA,c))*idx(2,cof(Nc))*𝑖/2");
}

#[test]
fn antisymmetric_trace_commutator_reduces_before_terminal_trace() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(color_t!(slot!(r.coad_da, a)), color_t!(slot!(r.coad_da, b))),
        color_t!(slot!(r.coad_da, c)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"f(coad(dA,a),coad(dA,b),coad(dA,c))*idx(2,cof(Nc))*𝑖/2");
}

#[test]
fn antisymmetric_trace_commutator_preserves_projector_sign() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(color_t!(slot!(r.coad_da, b)), color_t!(slot!(r.coad_da, a))),
        color_t!(slot!(r.coad_da, c)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"-𝑖/2*f(coad(dA,a),coad(dA,b),coad(dA,c))*idx(2,cof(Nc))");
}

#[test]
fn antisymmetric_chain_commutator_reduces_to_structure_constant() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        antisym!(color_t!(slot!(r.coad_da, a)), color_t!(slot!(r.coad_da, b))),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(dA,x),in,out))*f(coad(dA,a),coad(dA,b),coad(dA,x))*𝑖/2");
}

#[test]
fn cyclic_trace_structure_product_scans_nonleading_pairs() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, c)),
        color_t!(slot!(r.coad_da, d)),
        color_t!(slot!(r.coad_da, b)),
    ) * color_f!(
        slot!(r.coad_da, a),
        slot!(r.coad_da, d),
        slot!(r.coad_da, c),
    );

    let simplified = crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
        .unwrap()
        .simplify_color(ColorSimplifySettings::default().with_cof_dimension_invariants())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    let simplified = simplified.to_bare_ordered_string();

    assert!(
        !simplified.contains("trace(") && !simplified.contains("t(") && !simplified.contains("f("),
        "expected closed color structure to reduce to scalar invariants, got {simplified}"
    );
}

#[test]
fn trace_structure_product_orientation_sign_cancels() {
    test_initialize();
    let r = TestReps::new();
    let expr = color_f!(
        slot!(r.coad_da, q),
        slot!(r.coad_da, a),
        slot!(r.coad_da, b),
    ) * (trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, c)),
    ) + trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, c)),
    ));

    assert_color_zero(expr);
}

#[test]
fn coefficient_list_keeps_closed_traces_on_color_side() {
    test_initialize();
    let r = TestReps::new();
    let trace_factor = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, c)),
    );
    let expr = Atom::var(s!(x))
        * trace_factor
        * color_f!(
            slot!(r.coad_da, a),
            slot!(r.coad_da, b),
            slot!(r.coad_da, c),
        );

    let expanded = SymbolicTensor::infer(expr)
        .unwrap()
        .coefficient_list(TensorCollectFilter::Reps([
            ColorAdjoint {}.into(),
            ColorFundamental {}.into(),
            ColorSextet {}.into(),
        ]))
        .unwrap();

    assert_eq!(expanded.len(), 1);
    let (color_factor, residual) = &expanded[0];
    let color_factor = color_factor.expression().to_bare_ordered_string();
    assert!(
        color_factor.contains("trace(") && color_factor.contains("f("),
        "expected closed trace and structure constant to stay in color factor, got {color_factor}"
    );
    let residual = residual.expression().to_bare_ordered_string();
    assert_eq!(
        residual, "x",
        "expected residual factor to be color-free, got {}",
        residual,
    );
}

#[test]
fn color_simplify_defaults_match_simplify_color() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        color_t!(slot!(r.coad_da, a)),
    ) * chain!(
        slot!(r.cof_nc, k),
        slot!(r.cof_nc.dual(), l),
        color_t!(slot!(r.coad_da, a)),
    );

    assert_eq!(
        crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    );
}

#[test]
fn color_trace_evaluation_can_be_disabled() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
    );
    let expected = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
    );

    assert_eq!(
        crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(ColorSimplifySettings::default().without_trace_evaluation())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expected,
    );
}

#[test]
fn color_cross_chain_fierz_can_be_disabled() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        color_t!(slot!(r.coad_da, a)),
    ) * chain!(
        slot!(r.cof_nc, k),
        slot!(r.cof_nc.dual(), l),
        color_t!(slot!(r.coad_da, a)),
    );

    assert_eq!(
        crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(ColorSimplifySettings::default().without_cross_chain_fierz_expansion())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expr,
    );
}

#[test]
fn color_cross_chain_fierz_handles_longer_open_chains() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
    ) * chain!(
        slot!(r.cof_nc, k),
        slot!(r.cof_nc.dual(), l),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, c)),
    );
    let index = color_idx!(2, ColorFundamental {}.to_symbolic([Atom::var(s!(Nc))]));
    let expected = index.clone()
        * chain!(
            slot!(r.cof_nc, i),
            slot!(r.cof_nc.dual(), l),
            color_t!(slot!(r.coad_da, c)),
        )
        * chain!(
            slot!(r.cof_nc, k),
            slot!(r.cof_nc.dual(), j),
            color_t!(slot!(r.coad_da, b)),
        )
        - index
            * chain!(
                slot!(r.cof_nc, i),
                slot!(r.cof_nc.dual(), j),
                color_t!(slot!(r.coad_da, b)),
            )
            * chain!(
                slot!(r.cof_nc, k),
                slot!(r.cof_nc.dual(), l),
                color_t!(slot!(r.coad_da, c)),
            )
            / s!(Nc);

    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            (crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
                .unwrap()
                .simplify_color(crate::color::ColorSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression())
            .as_atom_view()
            .to_owned()
        )
        .unwrap()
        .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
        crate::tensor::SymbolicTensor::infer((expected).as_atom_view().to_owned())
            .unwrap()
            .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
    );
}

#[test]
fn symmetric_trace_d33_partial_contraction() {
    test_initialize();
    let r = TestReps::new();
    let left = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_da, a)),
            color_t!(slot!(r.coad_da, b)),
            color_t!(slot!(r.coad_da, c)),
        )
    );
    let right = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_da, a)),
            color_t!(slot!(r.coad_da, b)),
            color_t!(slot!(r.coad_da, d)),
        )
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((left * right).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"dA^(-1)*g(coad(dA,c),coad(dA,d))*gram(3,cof(Nc),cof(Nc))");
}

#[test]
fn symmetric_trace_d44_partial_contraction() {
    test_initialize();
    let r = TestReps::new();
    let left = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_da, a)),
            color_t!(slot!(r.coad_da, b)),
            color_t!(slot!(r.coad_da, c)),
            color_t!(slot!(r.coad_da, d)),
        )
    );
    let right = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_da, a)),
            color_t!(slot!(r.coad_da, b)),
            color_t!(slot!(r.coad_da, c)),
            color_t!(slot!(r.coad_da, e)),
        )
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((left * right).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"dA^(-1)*g(coad(dA,d),coad(dA,e))*gram(4,cof(Nc),cof(Nc))");
}

#[test]
fn symmetric_trace_d44_full_contraction() {
    test_initialize();
    let r = TestReps::new();
    let left = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_da, a)),
            color_t!(slot!(r.coad_da, b)),
            color_t!(slot!(r.coad_da, c)),
            color_t!(slot!(r.coad_da, d)),
        )
    );
    let right = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_da, a)),
            color_t!(slot!(r.coad_da, b)),
            color_t!(slot!(r.coad_da, c)),
            color_t!(slot!(r.coad_da, d)),
        )
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((left * right).as_atom_view().to_owned()).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"gram(4,cof(Nc),cof(Nc))");
}

mod feyncalc_reference;
mod form_reference;

fn colored_matrix_element() -> (Atom, Atom) {
    (
        parse!(
            "-G^2
            *(
                -g(mink(D,5),mink(D,6))*Q(2,mink(D,7))
                +g(mink(D,5),mink(D,6))*Q(3,mink(D,7))
                +g(mink(D,5),mink(D,7))*Q(2,mink(D,6))
                +g(mink(D,5),mink(D,7))*Q(4,mink(D,6))
                -g(mink(D,6),mink(D,7))*Q(3,mink(D,5))
                -g(mink(D,6),mink(D,7))*Q(4,mink(D,5))
            )
            *g(mink(D,4),mink(D,7))
            *t(coad(Nc^2-1,6),cof(Nc,5),dind(cof(Nc,4)))
            *f(coad(Nc^2-1,7),coad(Nc^2-1,8),coad(Nc^2-1,9))
            *g(bis(D,0),bis(D,5))
            *g(bis(D,1),bis(D,4))
            *g(mink(D,2),mink(D,5))
            *g(mink(D,3),mink(D,6))
            *g(coad(Nc^2-1,2),coad(Nc^2-1,7))
            *g(coad(Nc^2-1,3),coad(Nc^2-1,8))
            *g(coad(Nc^2-1,6),coad(Nc^2-1,9))
            *g(cof(Nc,0),dind(cof(Nc,5)))
            *g(cof(Nc,4),dind(cof(Nc,1)))
            *spenso::gamma(bis(D,5),bis(D,4),mink(D,4))
            *vbar(1,bis(D,1))
            *u(0,bis(D,0))
            *ϵbar(2,mink(D,2))
            *ϵbar(3,mink(D,3))",
            default_namespace = "spenso"
        ),
        parse_lit!(
            -4 * TR
                ^ 2 * Nc * G
                ^ 4 * (Nc - 1) * (Nc + 1) * (D - 2)
                ^ -2 * (-2 * dot(Q(0), Q(1)) * dot(Q(2), Q(2)) + dot(Q(0), Q(1)) * dot(Q(2), Q(3))
                    - 3 * dot(Q(0), Q(1)) * dot(Q(2), Q(4))
                    - 2 * dot(Q(0), Q(1)) * dot(Q(3), Q(3))
                    - 3 * dot(Q(0), Q(1)) * dot(Q(3), Q(4))
                    - 3 * dot(Q(0), Q(1)) * dot(Q(4), Q(4))
                    + 2 * dot(Q(0), Q(2)) * dot(Q(1), Q(2))
                    - dot(Q(0), Q(2)) * dot(Q(1), Q(3))
                    + dot(Q(0), Q(2)) * dot(Q(1), Q(4))
                    - dot(Q(0), Q(3)) * dot(Q(1), Q(2))
                    + 2 * dot(Q(0), Q(3)) * dot(Q(1), Q(3))
                    + dot(Q(0), Q(3)) * dot(Q(1), Q(4))
                    + dot(Q(0), Q(4)) * dot(Q(1), Q(2))
                    + dot(Q(0), Q(4)) * dot(Q(1), Q(3))
                    + 2 * dot(Q(0), Q(4)) * dot(Q(1), Q(4))
                    + D * dot(Q(0), Q(1)) * dot(Q(2), Q(2))
                    - D * dot(Q(0), Q(1)) * dot(Q(2), Q(3))
                    + D * dot(Q(0), Q(1)) * dot(Q(2), Q(4))
                    + D * dot(Q(0), Q(1)) * dot(Q(3), Q(3))
                    + D * dot(Q(0), Q(1)) * dot(Q(3), Q(4))
                    + D * dot(Q(0), Q(1)) * dot(Q(4), Q(4))
                    - D * dot(Q(0), Q(2)) * dot(Q(1), Q(2))
                    + D * dot(Q(0), Q(2)) * dot(Q(1), Q(3))
                    + D * dot(Q(0), Q(3)) * dot(Q(1), Q(2))
                    - D * dot(Q(0), Q(3)) * dot(Q(1), Q(3))),
            default_namespace = "spenso"
        ),
    )
}

#[test]
fn compact_printing() {
    test_initialize();

    let (atom, _) = colored_matrix_element();
    println!(
        "{}",
        atom.printer(SpensoPrintSettings::compact().nice_symbolica())
    );

    println!(
        "{}",
        atom.printer(
            SpensoPrintSettings {
                parens: true,
                with_dim: false,
                commas: false,
                index_subscripts: false,
                symbol_scripts: false,
                array_components: false,
            }
            .nice_symbolica()
        )
    );

    println!(
        "{}",
        atom.printer(SpensoPrintSettings::typst().nice_symbolica())
    )
}

#[test]
fn t_structure() {
    test_initialize();
    println!("{}", CS.t_strct::<AbstractIndex>(3, 8));

    let _ = crate::tensor::SymbolicTensor::infer((Atom::Zero).as_atom_view().to_owned())
        .unwrap()
        .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
}

#[test]
fn test_val() {
    test_initialize();
    let expr = parse_lit!(
        (G ^ 3
            * (g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                - g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9))))
            * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
            * g(mink(4, l(0)), mink(4, l(6)))
            * g(mink(4, l(1)), mink(4, l(7)))
            * g(mink(4, l(4)), mink(4, l(8)))
            * g(mink(4, l(5)), mink(4, l(9)))
            * g(bis(4, l(2)), bis(4, l(5)))
            * g(bis(4, l(3)), bis(4, l(6)))
            * g(dind(cof(Nc, l(6))), cof(Nc, 3))
            * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
            * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
            * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
            * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
            * spenso::gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
            * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
            * f(
                coad(Nc ^ 2 - 1, l(8)),
                coad(Nc ^ 2 - 1, l(11)),
                coad(Nc ^ 2 - 1, l(12))
            )
            * f(
                coad(Nc ^ 2 - 1, l(9)),
                coad(Nc ^ 2 - 1, l(10)),
                coad(Nc ^ 2 - 1, l(12))
            )
            * ubar(2, bis(4, l(2)))
            * v(3, bis(4, l(3)))
            * ϵ(0, mink(4, l(0)))
            * ϵ(1, mink(4, l(1)))
            * ϵbar(4, mink(4, l(4)))
            + G
            ^ 3 * (g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                - g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
                * g(mink(4, l(0)), mink(4, l(6)))
                * g(mink(4, l(1)), mink(4, l(7)))
                * g(mink(4, l(4)), mink(4, l(8)))
                * g(mink(4, l(5)), mink(4, l(9)))
                * g(bis(4, l(2)), bis(4, l(5)))
                * g(bis(4, l(3)), bis(4, l(6)))
                * g(dind(cof(Nc, l(6))), cof(Nc, 3))
                * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
                * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
                * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
                * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
                * spenso::gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
                * f(
                    coad(Nc ^ 2 - 1, l(8)),
                    coad(Nc ^ 2 - 1, l(10)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * f(
                    coad(Nc ^ 2 - 1, l(9)),
                    coad(Nc ^ 2 - 1, l(11)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * ubar(2, bis(4, l(2)))
                * v(3, bis(4, l(3)))
                * ϵ(0, mink(4, l(0)))
                * ϵ(1, mink(4, l(1)))
                * ϵbar(4, mink(4, l(4)))
                + G
            ^ 3 * (g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9)))
                - g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
                * g(mink(4, l(0)), mink(4, l(6)))
                * g(mink(4, l(1)), mink(4, l(7)))
                * g(mink(4, l(4)), mink(4, l(8)))
                * g(mink(4, l(5)), mink(4, l(9)))
                * g(bis(4, l(2)), bis(4, l(5)))
                * g(bis(4, l(3)), bis(4, l(6)))
                * g(dind(cof(Nc, l(6))), cof(Nc, 3))
                * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
                * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
                * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
                * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
                * spenso::gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
                * f(
                    coad(Nc ^ 2 - 1, l(8)),
                    coad(Nc ^ 2 - 1, l(9)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * f(
                    coad(Nc ^ 2 - 1, l(10)),
                    coad(Nc ^ 2 - 1, l(11)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * ubar(2, bis(4, l(2)))
                * v(3, bis(4, l(3)))
                * ϵ(0, mink(4, l(0)))
                * ϵ(1, mink(4, l(1)))
                * ϵbar(4, mink(4, l(4))))
            * (-G
                ^ 3 * (g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                    - g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9))))
                    * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                    * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                    * spenso::gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                    * f(
                        coad(Nc ^ 2 - 1, r(8)),
                        coad(Nc ^ 2 - 1, r(11)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, r(9)),
                        coad(Nc ^ 2 - 1, r(10)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * u(2, bis(4, r(2)))
                    * vbar(3, bis(4, r(3)))
                    * ϵ(4, mink(4, r(4)))
                    * ϵbar(0, mink(4, r(0)))
                    * ϵbar(1, mink(4, r(1)))
                    - G
                ^ 3 * (g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                    - g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                    * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                    * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                    * spenso::gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                    * f(
                        coad(Nc ^ 2 - 1, r(8)),
                        coad(Nc ^ 2 - 1, r(10)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, r(9)),
                        coad(Nc ^ 2 - 1, r(11)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * u(2, bis(4, r(2)))
                    * vbar(3, bis(4, r(3)))
                    * ϵ(4, mink(4, r(4)))
                    * ϵbar(0, mink(4, r(0)))
                    * ϵbar(1, mink(4, r(1)))
                    - G
                ^ 3 * (g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9)))
                    - g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                    * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                    * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                    * spenso::gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                    * f(
                        coad(Nc ^ 2 - 1, r(8)),
                        coad(Nc ^ 2 - 1, r(9)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, r(10)),
                        coad(Nc ^ 2 - 1, r(11)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * u(2, bis(4, r(2)))
                    * vbar(3, bis(4, r(3)))
                    * ϵ(4, mink(4, r(4)))
                    * ϵbar(0, mink(4, r(0)))
                    * ϵbar(1, mink(4, r(1)))),
        default_namespace = "spenso"
    );

    println!("{expr}");
    println!("Simplify_color");
    let dimension_value = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let dimension = Atom::var(symbol!("color_test_val_adjoint_dimension"));
    let indexed = expr
        .replace(dimension_value.to_pattern())
        .with(dimension.to_pattern());
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    let admitted = cooking.try_cook_indices(indexed.as_view()).unwrap();
    assert_eq!(
        cooking
            .uncook(admitted.as_view())
            .replace(dimension.to_pattern())
            .with(dimension_value.to_pattern()),
        expr
    );
    let metric_result = crate::tensor::SymbolicTensor::infer(admitted.clone())
        .unwrap()
        .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    println!(
        "Simplify_metrics: {}",
        cooking
            .uncook(metric_result.as_view())
            .replace(dimension.to_pattern())
            .with(dimension_value.to_pattern())
    );
    let simplified = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    let contracted = crate::tensor::SymbolicTensor::infer(simplified)
        .unwrap()
        .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    println!(
        "{:>}",
        cooking
            .uncook(contracted.as_view())
            .replace(dimension.to_pattern())
            .with(dimension_value.to_pattern())
            .expand()
            .to_dots()
    );
}

#[test]
fn ratio_simplify() {
    test_initialize();
    let expr = parse_lit!(
        G ^ 4 * ee
            ^ 2 * f(
                coad(ohoho, dummy(0)),
                coad(ohoho, dummy(1)),
                coad(ohoho, dummy(2))
            ) * t(
                coad(ohoho, dummy(0)),
                cof(ahaha, dummy(3)),
                dind(cof(ahaha, dummy(4)))
            ) * t(
                coad(ohoho, dummy(1)),
                cof(ahaha, dummy(5)),
                dind(cof(ahaha, dummy(3)))
            ) * t(
                coad(ohoho, dummy(2)),
                cof(ahaha, dummy(4)),
                dind(cof(ahaha, dummy(5)))
            ),
        default_namespace = "spenso"
    );

    let simplified =
        crate::tensor::SymbolicTensor::infer((expr.cook_indices()).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();

    assert_snapshot!(simplified.collect_with_map(|a| {
        matches!(a, AtomView::Var(a) if a.get_symbol() == CS.tr || a.get_symbol() == CS.ca || a.get_symbol() == CS.cf)
            || matches!(a, AtomView::Fun(f) if f.get_symbol() == CS.cas || f.get_symbol() == CS.idx || f.get_symbol() == CS.gram)
    }).into_inner().unwrap_collect().collect_factors().to_bare_ordered_string(), @"-𝑖/2*G^4*cas(2,coad(ohoho))*ee^2*idx(2,cof(ahaha))*ohoho");
}

#[test]
fn structure_pair_with_closed_generator_chain_matches_form() {
    test_initialize();
    let expr = parse_lit!(
        f(
            coad(Nc ^ 2 - 1, c4),
            coad(Nc ^ 2 - 1, c0),
            coad(Nc ^ 2 - 1, c2)
        ) * f(
            coad(Nc ^ 2 - 1, c6),
            coad(Nc ^ 2 - 1, c2),
            coad(Nc ^ 2 - 1, c0)
        ) * t(coad(Nc ^ 2 - 1, c4), cof(Nc, i0), dind(cof(Nc, i2)))
            * t(coad(Nc ^ 2 - 1, c6), cof(Nc, i2), dind(cof(Nc, i0))),
        default_namespace = "spenso"
    );

    let dimension_value = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let dimension = Atom::var(symbol!(
        "structure_pair_with_closed_generator_chain_matches_form_dimension"
    ));
    let admitted = expr
        .replace(dimension_value.to_pattern())
        .with(dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(dimension.to_pattern())
            .with(dimension_value.to_pattern()),
        expr
    );
    assert_snapshot!(
        crate::tensor::SymbolicTensor::infer(admitted).unwrap().simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().replace(dimension.to_pattern()).with(dimension_value.to_pattern()).to_bare_ordered_string(),
        @"(-1+Nc^2)*-1*cas(2,coad(-1+Nc^2))*idx(2,cof(Nc))"
    );
}

#[test]
fn minus_sign() {
    test_initialize();
    let expr1 = parse_lit!(
        f(coad(8, hedge(5)), coad(8, hedge(9)), coad(8, hedge(13)))
            * g(coad(8, hedge(12)), coad(8, hedge(13)))
            * g(coad(8, hedge(4)), coad(8, hedge(5)))
            * g(coad(8, hedge(8)), coad(8, hedge(9)))
            * g(cof(3, hedge(10)), dind(cof(3, hedge(11))))
            * g(cof(3, hedge(11)), dind(cof(3, hedge(15))))
            * g(cof(3, hedge(15)), dind(cof(3, hedge(14))))
            * g(cof(3, hedge(16)), dind(cof(3, hedge(17))))
            * g(cof(3, hedge(17)), dind(cof(3, hedge(7))))
            * g(cof(3, hedge(2)), dind(cof(3, hedge(3))))
            * g(cof(3, hedge(7)), dind(cof(3, hedge(6))))
            * t(
                coad(8, hedge(12)),
                cof(3, hedge(14)),
                dind(cof(3, hedge(16)))
            )
            * t(coad(8, hedge(4)), cof(3, hedge(6)), dind(cof(3, hedge(2))))
            * t(coad(8, hedge(8)), cof(3, hedge(3)), dind(cof(3, hedge(10)))),
        default_namespace = "spenso"
    );

    let expr2 = parse_lit!(
        f(coad(8, hedge(5)), coad(8, hedge(9)), coad(8, hedge(13)))
            * g(coad(8, hedge(12)), coad(8, hedge(13)))
            * g(coad(8, hedge(4)), coad(8, hedge(5)))
            * g(coad(8, hedge(8)), coad(8, hedge(9)))
            * g(cof(3, hedge(11)), dind(cof(3, hedge(10))))
            * g(cof(3, hedge(14)), dind(cof(3, hedge(15))))
            * g(cof(3, hedge(15)), dind(cof(3, hedge(11))))
            * g(cof(3, hedge(17)), dind(cof(3, hedge(16))))
            * g(cof(3, hedge(3)), dind(cof(3, hedge(2))))
            * g(cof(3, hedge(6)), dind(cof(3, hedge(7))))
            * g(cof(3, hedge(7)), dind(cof(3, hedge(17))))
            * t(
                coad(8, hedge(12)),
                cof(3, hedge(16)),
                dind(cof(3, hedge(14)))
            )
            * t(coad(8, hedge(4)), cof(3, hedge(2)), dind(cof(3, hedge(6))))
            * t(coad(8, hedge(8)), cof(3, hedge(10)), dind(cof(3, hedge(3)))),
        default_namespace = "spenso"
    );

    // println!(
    //     "{}",
    //     (expr1
    //         .cook_indices()
    //         .canonize(AbstractIndex::Dummy)
    //         .expect("test expression should canonicalize")
    //         / expr2
    //             .cook_indices()
    //             .canonize(AbstractIndex::Dummy)
    //             .expect("test expression should canonicalize"))
    //     .cancel()
    // );
    println!(
        "{}\n",
        crate::tensor::SymbolicTensor::infer(expr1.cook_indices())
            .unwrap()
            .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
            .expect("test expression should canonicalize")
    );
    println!(
        "{}\n",
        crate::tensor::SymbolicTensor::infer(expr2.cook_indices())
            .unwrap()
            .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
            .expect("test expression should canonicalize")
    );
    let residual = crate::tensor::SymbolicTensor::infer(
        ((crate::tensor::SymbolicTensor::infer(expr1.cook_indices())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            + crate::tensor::SymbolicTensor::infer(expr2.cook_indices())
                .unwrap()
                .simplify_color(crate::color::ColorSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression())
        .expand())
        .as_atom_view()
        .to_owned(),
    )
    .unwrap()
    .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
    .unwrap()
    .resolved()
    .unwrap()
    .into_expression()
    .cook_indices();

    println!("{}", residual);

    let residual = crate::tensor::SymbolicTensor::infer(
        (residual
            .canonize(AbstractIndex::Dummy)
            .expect("test expression should canonicalize"))
        .as_atom_view()
        .to_owned(),
    )
    .unwrap()
    .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
    .unwrap()
    .resolved()
    .unwrap()
    .into_expression()
    .expand();
    assert!(residual.is_zero());
}

mod failing {
    use super::*;

    #[test]
    fn test_color_matrix_element() {
        use crate::shorthands::{UndoShorthands, chain::Chain};

        test_initialize();
        let spin_sum_rule = parse!(
            "
                g(coad(Nc^2-1, left(3)), coad(Nc^2-1, right(3)))
                    * g(coad(Nc^2-1, left(2)), coad(Nc^2-1, right(2)))
                    * g(cof(Nc, right(0)), dind(cof(Nc, left(0))))
                    * g(cof(Nc, left(1)), dind(cof(Nc, right(1))))",
            default_namespace = "spenso"
        );

        let amplitude_color = parse!(
            "
                t(coad(Nc^2-1, 6), cof(Nc, 5), dind(cof(Nc, 4)))
                    * f(coad(Nc^2-1, 7), coad(Nc^2-1, 8), coad(Nc^2-1, 9))
                    * g(coad(Nc^2-1, 2), coad(Nc^2-1, 7))
                    * g(coad(Nc^2-1, 3), coad(Nc^2-1, 8))
                    * g(coad(Nc^2-1, 6), coad(Nc^2-1, 9))
                    * g(cof(Nc, 0), dind(cof(Nc, 5)))
                    * g(cof(Nc, 4), dind(cof(Nc, 1)))",
            default_namespace = "spenso"
        );
        // Network dimensions are integers or symbols. Previously the compound
        // adjoint dimension was silently discarded as metadata by the parser.
        // Admit it explicitly for index tooling, then restore the original algebra.
        let adjoint_dimension = Atom::var(symbol!("matrix_element_adjoint_dimension"));
        let dimension_value = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
        let indexed_color = amplitude_color
            .replace(dimension_value.clone())
            .with(adjoint_dimension.clone());
        assert_eq!(
            indexed_color
                .replace(adjoint_dimension.clone())
                .with(dimension_value.clone()),
            amplitude_color
        );
        // Scope while dimensions are admitted: compound dimensions would leave
        // adjoint-color indices untouched and collide between the two amplitudes.
        let amplitude_color_left = indexed_color
            .wrap_indices(symbol!("spenso::left"))
            .replace(adjoint_dimension.clone())
            .with(dimension_value.clone());
        let amplitude_color_right = indexed_color
            .dirac_adjoint::<AbstractIndex>(false)
            .unwrap()
            .wrap_indices(symbol!("spenso::right"))
            .replace(adjoint_dimension.clone())
            .with(dimension_value.clone());
        println!("left{amplitude_color_left}");

        println!("right{amplitude_color_right}");
        let amp_squared_color = amplitude_color_left * spin_sum_rule * amplitude_color_right;
        let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
        let indexed = amp_squared_color
            .replace(dimension_value.to_pattern())
            .with(adjoint_dimension.to_pattern());
        let admitted = cooking.try_cook_indices(indexed.as_view()).unwrap();
        assert_eq!(cooking.uncook(admitted.as_view()), indexed);
        let simplified_color = crate::tensor::SymbolicTensor::infer(admitted)
            .unwrap()
            .contract(crate::tensor::ContractionSettings::default().without_rank_one_tensors())
            .unwrap()
            .resolved()
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .replace(adjoint_dimension.to_pattern())
            .with(dimension_value.to_pattern());
        println!("simplified_color={}", simplified_color);

        let spin_sum_rule_src = parse_lit!(
            spenso::vbar(1, bis(D, left(1)))
                * spenso::v(1, bis(D, right(1)))
                * spenso::u(0, bis(D, left(0)))
                * spenso::ubar(0, bis(D, right(0)))
                * spenso::ϵbar(2, mink(D, left(2)))
                * spenso::ϵ(2, mink(D, right(2)))
                * spenso::ϵbar(3, mink(D, left(3)))
                * spenso::ϵ(3, mink(D, right(3))),
            default_namespace = "spenso"
        );

        let spin_sum_rule_trg = parse!(
            "
                1/4*1/(D-2)^2*
                (
                    (-1) * spenso::gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))*Q(1,mink(D,1337))
                    * spenso::gamma(bis(D,right(0)),bis(D,left(0)),mink(D,1338))*Q(0,mink(D,1338))
                    * (-1) * g(mink(D,left(2)),mink(D,right(2)))
                    * (-1) * g(mink(D,left(3)),mink(D,right(3)))
                    )
                    * (
                        g(coad(Nc^2-1, left(3)), coad(Nc^2-1, right(3)))
                        * g(coad(Nc^2-1, left(2)), coad(Nc^2-1, right(2)))
                        * g(cof(Nc, right(0)), dind(cof(Nc, left(0))))
                        * g(cof(Nc, left(1)), dind(cof(Nc, right(1))))
                    )",
            default_namespace = "spenso"
        );

        let (amplitude, tgt) = colored_matrix_element();

        let indexed_amplitude = amplitude
            .replace(dimension_value.clone())
            .with(adjoint_dimension.clone());
        assert_eq!(
            indexed_amplitude
                .replace(adjoint_dimension.clone())
                .with(dimension_value.clone()),
            amplitude
        );
        let amplitude_left = indexed_amplitude
            .wrap_indices(symbol!("spenso::left"))
            .replace(adjoint_dimension.clone())
            .with(dimension_value.clone());

        println!("Amplitude left:\n{}", amplitude_left.collect_factors());

        println!(
            "Amplitude left cooked:\n{}",
            amplitude_left.collect_factors().cook_indices()
        );
        let amplitude_right = indexed_amplitude
            .wrap_indices(symbol!("spenso::right"))
            .replace(adjoint_dimension.clone())
            .with(dimension_value.clone());

        println!("Amplitude right:\n{}", amplitude_right.factor());

        // The full amplitude needs the same explicit dimension admission as
        // its color-only part above. Reuse the checked result for both consumers.
        let indexed_right = amplitude_right
            .replace(dimension_value.clone())
            .with(adjoint_dimension.clone());
        assert_eq!(
            indexed_right
                .replace(adjoint_dimension.clone())
                .with(dimension_value.clone()),
            amplitude_right
        );
        let amplitude_right_adjoint = indexed_right
            .dirac_adjoint::<AbstractIndex>(false)
            .unwrap()
            .replace(adjoint_dimension.to_pattern())
            .with(dimension_value.to_pattern());
        println!(
            "Amplitude right conj:\n{}",
            amplitude_right_adjoint.factor()
        );

        let mut amp_squared = amplitude_left * amplitude_right_adjoint;

        println!("Amplitude squared:\n{}", amp_squared.factor());

        // let dangling_atoms = amp_squared.list_dangling::<AbstractIndex>();

        // assert_eq!(dangling_atoms.len(), 8);

        amp_squared = amp_squared
            .expand()
            .replace(spin_sum_rule_src.to_pattern())
            .with(spin_sum_rule_trg.to_pattern());

        println!("Amplitude squared spin-summed:\n{}", amp_squared);

        let indexed = amp_squared
            .replace(dimension_value.to_pattern())
            .with(adjoint_dimension.to_pattern());
        let cooked = cooking.try_cook_indices(indexed.as_view()).unwrap();
        assert_eq!(
            cooking
                .uncook(cooked.as_view())
                .replace(adjoint_dimension.to_pattern())
                .with(dimension_value.to_pattern()),
            amp_squared
        );
        let admitted = cooked.chainify(crate::representations::Bispinor {}.into());
        assert_eq!(
            admitted.undo_chain::<AbstractIndex>().unwrap(),
            cooked.undo_chain::<AbstractIndex>().unwrap()
        );
        let color = crate::tensor::SymbolicTensor::infer(admitted)
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap();
        println!(
            "Color-simplified amplitude squared:\n{}",
            color.resolved().unwrap().expression
        );
        let gamma = color
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        let mut simplified_amp_squared = cooking
            .uncook(gamma.as_view())
            .replace(adjoint_dimension.to_pattern())
            .with(dimension_value.to_pattern());
        println!(
            "Gamma+color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.to_dots();

        assert_ne!(tgt, simplified_amp_squared.factor());
    }

    #[test]
    fn test_color_matrix_element_two() {
        use crate::shorthands::chain::Chain;

        test_initialize();

        let spin_sum_rule_src = parse_lit!(
            vbar(1, bis(D, left(1)))
                * v(1, bis(D, right(1)))
                * u(0, bis(D, left(0)))
                * ubar(0, bis(D, right(0)))
                * ϵbar(2, mink(D, left(2)))
                * ϵ(2, mink(D, right(2)))
                * ϵbar(3, mink(D, left(3)))
                * ϵ(3, mink(D, right(3))),
            default_namespace = "spenso"
        );

        let spin_sum_rule_trg = parse!(
            "
                1/4*1/(D-2)^2*
                (
                    (-1) * spenso::gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))*Q(1,mink(D,1337))
                    * spenso::gamma(bis(D,right(0)),bis(D,left(0)),mink(D,1338))*Q(0,mink(D,1338))
                    * (-1) * g(mink(D,left(2)),mink(D,right(2)))
                    * (-1) * g(mink(D,left(3)),mink(D,right(3)))
                    )
                    ",
            default_namespace = "spenso"
        );

        let (amplitude, tgt) = colored_matrix_element();

        // Keep the supported symbolic-dimension boundary explicit, without
        // losing the original Nc dependence or omitting adjoint ports.
        let adjoint_dimension = Atom::var(symbol!("matrix_element_two_adjoint_dimension"));
        let dimension_value = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
        let indexed_amplitude = amplitude
            .replace(dimension_value.clone())
            .with(adjoint_dimension.clone());
        assert_eq!(
            indexed_amplitude
                .replace(adjoint_dimension.clone())
                .with(dimension_value.clone()),
            amplitude
        );
        let amplitude_left = indexed_amplitude
            .wrap_dummies::<AbstractIndex>(symbol!("spenso::left"))
            .unwrap()
            .replace(adjoint_dimension.clone())
            .with(dimension_value.clone());

        println!("Amplitude left:\n{}", amplitude_left.collect_factors());

        let amplitude_right = indexed_amplitude
            .wrap_dummies::<AbstractIndex>(symbol!("spenso::right"))
            .unwrap()
            .replace(adjoint_dimension.clone())
            .with(dimension_value.clone());

        println!("Amplitude right:\n{}", amplitude_right.conj().factor());

        let amp_squared = amplitude_left * amplitude_right.conj();
        let spin_summed = amp_squared
            .replace(spin_sum_rule_src.to_pattern())
            .with(spin_sum_rule_trg.to_pattern());
        // Scalar conjugation does not implement the tensor adjoint or turn
        // conjugated spinors into the spin-sum rule's conjugate wavefunctions.
        assert_eq!(spin_summed, amp_squared);

        let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
        let indexed = spin_summed
            .replace(dimension_value.to_pattern())
            .with(adjoint_dimension.to_pattern());
        let cooked = cooking.try_cook_indices(indexed.as_view()).unwrap();
        assert_eq!(
            cooking
                .uncook(cooked.as_view())
                .replace(adjoint_dimension.to_pattern())
                .with(dimension_value.to_pattern()),
            spin_summed
        );
        let admitted = cooked.chainify(crate::representations::Bispinor {}.into());
        // conj(Q(slot)) is an opaque scalar function at the typed boundary.
        // Its surrounding metric summands therefore have different interfaces;
        // accepting them used to silently discard the conjugated vector ports.
        let error = crate::tensor::SymbolicTensor::infer(admitted).unwrap_err();
        assert!(
            matches!(error, crate::tensor::inference::TensorInferenceError::Invalid(ref reason)
            if reason.contains("compatible tensor interfaces"))
        );
        // Retain the original negative reference. The supported tensor-adjoint
        // route is exercised by test_color_matrix_element above.
        assert_ne!(tgt, spin_summed.factor());
    }
}

#[test]
fn color_trace_metric_closures_reduce_before_terminal_invariants() {
    test_initialize();
    let expr = parse!(
        "trace(cof(Nc),t(coad(Na,aa),in,out),t(coad(Na,bb),in,out),t(coad(Na,cc),in,out),t(coad(Na,dd),in,out))*g(coad(Na,aa),coad(Na,cc))*g(coad(Na,bb),coad(Na,dd))",
        default_namespace = "spenso"
    );
    let result = crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .to_cof_dimension_invariants()
        .replace(parse!("Na", default_namespace = "spenso"))
        .with(parse!("Nc^2-1", default_namespace = "spenso"))
        .to_cof_dimension_invariants();
    let expected = parse!("-(Nc^2-1)/(4*Nc)", default_namespace = "spenso");
    assert_eq!((result - expected).expand(), Atom::Zero);
}

#[test]
fn symmetric_color_trace_with_contracted_pair_reuses_casimir_rules() {
    test_initialize();
    let expr = parse!(
        "trace(cof(Nc),sym(t(coad(Na,aa),in,out),t(coad(Na,bb),in,out),t(coad(Na,cc),in,out),t(coad(Na,dd),in,out)))*g(coad(Na,aa),coad(Na,cc))*g(coad(Na,bb),coad(Na,dd))",
        default_namespace = "spenso"
    );
    let result = crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression()
        .to_cof_dimension_invariants()
        .replace(parse!("Na", default_namespace = "spenso"))
        .with(parse!("Nc^2-1", default_namespace = "spenso"))
        .to_cof_dimension_invariants();
    let expected = parse!("(Nc^2-1)*(2*Nc^2-3)/(12*Nc)", default_namespace = "spenso");
    assert_eq!((result - expected).expand(), Atom::Zero);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            (crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
                .unwrap()
                .simplify_color(crate::color::ColorSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression())
            .as_atom_view()
            .to_owned()
        )
        .unwrap()
        .simplify_color(crate::color::ColorSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
        crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    );
}

#[test]
fn cooked_dimensions_preserve_color_invariants_and_zero() {
    use crate::{CookMode, CookSettings};
    test_initialize();
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let n = Atom::var(CS.nc);
    let adjoint = n.clone().pow(Atom::num(2)) - Atom::num(1);
    let cof = ColorFundamental {}.to_symbolic([n.clone()]);
    let coad = ColorAdjoint {}.to_symbolic([adjoint.clone()]);
    let supported = [
        color_idx!(2, cof.clone()),
        color_cas!(2, cof.clone()),
        color_cas!(2, coad.clone()),
        color_gram!(3, cof.clone(), cof.clone()),
        color_gram!(4, cof.clone(), cof),
    ];
    for source in supported {
        let encoded = cooking.cook(source.as_view());
        let actual = cooking.uncook(encoded.to_cof_dimension_invariants().as_view());
        let expected = source.to_cof_dimension_invariants();
        assert_eq!(actual, expected);
        assert_eq!(actual.to_cof_dimension_invariants(), expected);
    }
    let source = color_cas!(2, coad) - n;
    let result = SymbolicTensor::infer(cooking.cook(source.as_view()))
        .unwrap()
        .simplify_color(ColorSimplifySettings::default().with_cof_dimension_invariants())
        .unwrap()
        .resolved()
        .unwrap();
    assert!(result.expression().is_zero());
    assert_eq!(cooking.uncook(result.expression().as_view()), Atom::Zero);
    for source in [
        parse_lit!(cas(2, coad(8)), default_namespace = "spenso"),
        parse_lit!(cas(2, cof(3)), default_namespace = "spenso"),
    ] {
        assert_eq!(cooking.cook(source.as_view()), source);
        assert_eq!(
            cooking.cook(source.as_view()).to_cof_dimension_invariants(),
            source.to_cof_dimension_invariants()
        );
    }
}

#[test]
fn cooked_compound_dimension_color_contraction_matches_raw_invariant_identity() {
    use crate::{CookMode, CookSettings};
    test_initialize();
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let source = parse!(
        "f(coad(Nc^2-1,a),coad(Nc^2-1,b),coad(Nc^2-1,c))^2",
        default_namespace = "spenso"
    );
    let expected = Atom::var(CS.nc) * (Atom::var(CS.nc).pow(Atom::num(2)) - Atom::num(1));
    let result = SymbolicTensor::infer(cooking.cook(source.as_view()))
        .unwrap()
        .simplify_color(ColorSimplifySettings::default().with_cof_dimension_invariants())
        .unwrap()
        .resolved()
        .unwrap();
    let decoded = cooking.uncook(result.expression().as_view());
    assert_eq!((decoded - expected).together().cancel(), Atom::Zero);
    assert_eq!(
        result
            .clone()
            .simplify_color(ColorSimplifySettings::default().with_cof_dimension_invariants())
            .unwrap()
            .resolved()
            .unwrap(),
        result
    );
}

#[test]
fn scalar_invariant_prefactor_preserves_closed_color_trace_collection() {
    test_initialize();
    let input = symbolica::atom::Atom::parse(
        "idx(2,cof(2))*(-trace(cof(2),cyclic(t(coad(3,a),in,out),t(coad(3,c),in,out)))^2/2+trace(cof(2),cyclic(t(coad(3,a),in,out),t(coad(3,a),in,out),t(coad(3,c),in,out),t(coad(3,c),in,out))))",
        "spenso", symbolica::parser::ParseSettings::symbolica()).unwrap();
    let spectator = (symbolica::atom::Atom::var(symbolica::symbol!("closed_color_x"))
        + symbolica::atom::Atom::var(symbolica::symbol!("closed_color_y")))
    .pow(3);
    for coefficient in [symbolica::atom::Atom::one(), spectator] {
        let tensor = crate::tensor::SymbolicTensor::infer(&input * &coefficient).unwrap();
        let result = tensor
            .simplify_color(ColorSimplifySettings::default().with_cof_dimension_invariants())
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(
            result.into_expression(),
            symbolica::atom::Atom::num((3, 8)) * coefficient
        );
    }
}
