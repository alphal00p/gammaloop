#![allow(dead_code)]

use idenso::{
    representations::initialize,
    schoonschip::{Schoonschip, SchoonschipContractionOrder, SchoonschipSettings},
    tensor::{SymbolicNet, SymbolicNetParse},
};
use spenso::{
    network::parsing::ParseSettings,
    network::tags::SPENSO_TAG as T,
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName, Representation},
        slot::IsAbstractSlot,
    },
};
use symbolica::{
    LicenseManager,
    atom::{Atom, AtomCore},
    function,
    id::Pattern,
    parse, symbol,
    transformer::Transformer,
};

pub fn nested_dot_expression() -> Atom {
    activate_symbolica_license();

    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let q = T.rank_one_tensor_symbol("Q");

    let p1 = function!(p, 1, mink.slot::<AbstractIndex, _>(1).to_atom());
    let p2 = function!(p, 2, mink.slot::<AbstractIndex, _>(1).to_atom());
    let p2_2 = function!(p, 2, mink.slot::<AbstractIndex, _>(2).to_atom());

    let q2 = function!(
        q,
        2,
        symbol!("bla"),
        mink.slot::<AbstractIndex, _>(1).to_atom()
    );
    let q2_2 = function!(
        q,
        2,
        symbol!("bla"),
        mink.slot::<AbstractIndex, _>(2).to_atom()
    );
    let q3_2 = function!(q, 3, mink.slot::<AbstractIndex, _>(2).to_atom());

    &p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2))
}

pub fn run_schoonschip(expr: Atom, settings: &SchoonschipSettings) -> Atom {
    expr.schoonschip_with_net::<false, true, AbstractIndex>(settings)
}

pub fn assert_benchmark_outputs_match() {
    let expr = nested_dot_expression();
    let expected = run_schoonschip(expr.clone(), &SchoonschipSettings::full());

    for (name, settings) in benchmark_settings() {
        let result = run_schoonschip(expr.clone(), &settings);
        assert_eq!(
            result, expected,
            "benchmark mode {name} produced a different output"
        );
    }
}

#[allow(dead_code)]
pub fn checked_nested_dot_expression() -> Atom {
    assert_benchmark_outputs_match();
    nested_dot_expression()
}

pub fn benchmark_settings() -> Vec<(&'static str, SchoonschipSettings)> {
    vec![
        ("depth_first_depth_1", SchoonschipSettings::partial()),
        (
            "breadth_first_depth_1",
            SchoonschipSettings::breadth_first(Some(1)),
        ),
        ("full", SchoonschipSettings::full()),
    ]
}

pub struct NetworkVertexFixture {
    pub input: Atom,
    pub gluon_rule: Pattern,
    pub vertex_count: usize,
}

pub struct BareVertexFixture {
    pub input: Atom,
    pub rhs_subs: Pattern,
    pub vertex_count: usize,
}

pub fn network_vertex_fixture_5() -> NetworkVertexFixture {
    network_vertex_fixture(5)
}

pub fn network_vertex_fixture_8() -> NetworkVertexFixture {
    network_vertex_fixture(8)
}

pub fn network_vertex_fixture_with_count(vertex_count: usize) -> NetworkVertexFixture {
    network_vertex_fixture(vertex_count)
}

fn network_vertex_fixture(vertex_count: usize) -> NetworkVertexFixture {
    activate_symbolica_license();
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);
    symbol!("k"; tags=["spenso::tensor","spenso::rank1"]);

    let mut vertices = vec![
        parse!("vx(1,-k(0), k(0)-k(1), k(1), k(10), spenso::mink(4,mu1), spenso::mink(4,mu8))"),
        parse!(
            "vx(2,-k(1), k(2), k(1)-k(2), spenso::mink(4,mu1), spenso::mink(4,mu2), spenso::mink(4,mu9))"
        ),
        parse!(
            "vx(3,-k(2), k(3), k(2)-k(3), spenso::mink(4,mu2), spenso::mink(4,mu3), spenso::mink(4,mu10))"
        ),
        parse!(
            "vx(4,-k(3), k(4), k(3)-k(4), spenso::mink(4,mu3), spenso::mink(4,mu4), spenso::mink(4,mu11))"
        ),
        parse!("vx(5,-k(4), k(0), k(4)-k(0), spenso::mink(4,mu4), k(20), spenso::mink(4,mu5))"),
    ];
    match vertex_count {
        5 => {}
        6..=8 => {
            vertices.extend([
                parse!(
                    "vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), spenso::mink(4,mu5), spenso::mink(4,mu11), spenso::mink(4,mu6))"
                ),
                parse!(
                    "vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), spenso::mink(4,mu6), spenso::mink(4,mu10), spenso::mink(4,mu7))"
                ),
                parse!(
                    "vx(8,-k(2)+k(0), -k(1)+k(2), k(1)-k(0), spenso::mink(4,mu7), spenso::mink(4,mu9), spenso::mink(4,mu8))"
                ),
            ]);
            vertices.truncate(vertex_count);
        }
        _ => panic!("unsupported network vertex count: {vertex_count}"),
    }

    let input = vertices
        .into_iter()
        .reduce(|acc, vertex| acc * vertex)
        .unwrap();
    let gluon_rule = parse!(
        "(- spenso::g(mu1_,mu3_) * spenso::g(k1_,mu2_)
                + spenso::g(mu1_,mu2_) * spenso::g(k1_,mu3_)
                + spenso::g(mu2_,mu3_) * spenso::g(k2_,mu1_)
                - spenso::g(mu1_,mu2_) * spenso::g(k2_,mu3_)
                - spenso::g(mu2_,mu3_) * spenso::g(k3_,mu1_)
                + spenso::g(mu1_,mu3_) * spenso::g(k3_,mu2_)
                )"
    )
    .to_pattern();

    NetworkVertexFixture {
        input,
        gluon_rule,
        vertex_count,
    }
}

pub fn bare_vertex_fixture_5() -> BareVertexFixture {
    bare_vertex_fixture(5)
}

pub fn bare_vertex_fixture_8() -> BareVertexFixture {
    bare_vertex_fixture(8)
}

fn bare_vertex_fixture(vertex_count: usize) -> BareVertexFixture {
    activate_symbolica_license();

    let _d = symbol!("d"; Symmetric, Linear);
    let _muw1 = symbol!("muw1_", tags = ["spenso::index"]);
    symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    let input = match vertex_count {
        5 => parse!(
            "vx(1,-k0, k0-k1, k1, k10, mu1, mu8)
            vx(2,-k1, k2, k1-k2, mu1, mu2, mu9)
            vx(3,-k2, k3, k2-k3, mu2, mu3, mu10)
            vx(4,-k3, k4, k3-k4, mu3, mu4, mu11)
            vx(5,-k4, k0, k4-k0, mu4, k20, mu5)"
        ),
        8 => parse!(
            "vx(1,-k0, k0-k1, k1, k10, mu1, mu8)
            vx(2,-k1, k2, k1-k2, mu1, mu2, mu9)
            vx(3,-k2, k3, k2-k3, mu2, mu3, mu10)
            vx(4,-k3, k4, k3-k4, mu3, mu4, mu11)
            vx(5,-k4, k0, k4-k0, mu4, k20, mu5)
            vx(6,-k4+k0, -k3+k4, k3-k0, mu5, mu11, mu6)
            vx(7,-k3+k0, -k2+k3, k2-k0, mu6, mu10, mu7)
            vx(8,-k2+k0, -k1+k2, k1-k0, mu7, mu9, mu8)"
        ),
        _ => panic!("unsupported bare symbolic vertex count: {vertex_count}"),
    };

    let gluon_rule = parse!(
        "(- d(mu1_, mu3_) * d(k1_, mu2_)
                + d(mu1_, mu2_) * d(k1_, mu3_)
                + d(mu2_, mu3_) * d(k2_, mu1_)
                - d(mu1_, mu2_) * d(k2_, mu3_)
                - d(mu2_, mu3_) * d(k3_, mu1_)
                + d(mu1_, mu3_) * d(k3_, mu2_)
                )"
    )
    .to_pattern();

    let rhs_subs = Pattern::Transformer(Box::new((
        Some(gluon_rule),
        vec![Transformer::Map(Box::new(move |input, _state, out| {
            *out = input
                .expand()
                .replace(parse!("d(muw1_,k1_)*d(muw1_,k2_)"))
                .level_range((0, Some(0)))
                .repeat()
                .with(parse!("d(k1_,k2_)"))
                .replace(parse!("d(muw1_,k1_)^2"))
                .level_range((0, Some(0)))
                .with(parse!("d(k1_,k1_)"))
                .replace(parse!("d(muw1_,muw1_)"))
                .level_range((0, Some(0)))
                .with(4);
            Ok(())
        }))],
    )));

    BareVertexFixture {
        input,
        rhs_subs,
        vertex_count,
    }
}

pub fn network_vertex_substitution(fixture: NetworkVertexFixture) -> Atom {
    let mut result = fixture.input;

    for vertex_id in 1..=fixture.vertex_count {
        let gluon_rule = fixture.gluon_rule.clone();
        result = result
            .replace(parse!(format!(
                "vx({vertex_id}, k1_, k2_, k3_, mu1_, mu2_, mu3_)"
            )))
            .level_range((0, Some(0)))
            .rhs_cache_size(1000)
            .with_map(move |matches| {
                gluon_rule
                    .replace_wildcards_with_matches(matches)
                    .normalize_dots()
            });
    }

    result
}

pub fn bare_vertex_substitution(mut fixture: BareVertexFixture) -> Atom {
    for vertex_id in 1..=fixture.vertex_count {
        fixture.input = fixture
            .input
            .replace(parse!(format!(
                "vx({vertex_id}, k1_, k2_, k3_, mu1_, mu2_, mu3_)"
            )))
            .level_range((0, Some(0)))
            .rhs_cache_size(1000)
            .with(&fixture.rhs_subs);
        fixture.input = fixture.input.expand();
        fixture.input = fixture
            .input
            .replace(parse!("d(muw1_, k1_)*vx_(x___,muw1_,y___)"))
            .level_range((0, Some(0)))
            .repeat()
            .with(parse!("vx(x___,k1_,y___)"));
    }

    fixture.input
}

pub fn network_substituted_5() -> Atom {
    network_vertex_substitution(network_vertex_fixture_5())
}

pub fn network_substituted_8() -> Atom {
    network_vertex_substitution(network_vertex_fixture_8())
}

pub fn network_normalize_substituted(expr: Atom) -> Atom {
    expr.normalize_dots()
}

pub fn network_parse_normalized(expr: Atom) -> SymbolicNet<AbstractIndex> {
    expr.parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
        depth_limit: Some(1),
        take_first_term_from_sum: false,
        parse_inner_products: true,
        parse_composite_scalars_as_tensors: true,
        ..Default::default()
    })
    .unwrap()
}

pub fn network_schoonschip_substituted(expr: Atom) -> Atom {
    network_schoonschip_substituted_with_order(expr, SchoonschipContractionOrder::default())
}

pub fn network_schoonschip_substituted_with_order(
    expr: Atom,
    order: SchoonschipContractionOrder,
) -> Atom {
    let mut result = expr.schoonschip_with_net::<false, false, AbstractIndex>(
        &SchoonschipSettings::partial()
            .with_expanded_contracted_sums()
            .with_contraction_order(order),
    );

    let needs_cleanup = matches!(
        order,
        SchoonschipContractionOrder::LargestDegree
            | SchoonschipContractionOrder::MinLargestOperandBytes
            | SchoonschipContractionOrder::MinProductTerms
            | SchoonschipContractionOrder::MinProductBytes
    );

    if needs_cleanup {
        let cleanup_settings = SchoonschipSettings::partial()
            .with_expanded_contracted_sums()
            .with_contraction_order(SchoonschipContractionOrder::SmallestDegree);
        for _ in 0..4 {
            let next =
                result.schoonschip_with_net::<false, false, AbstractIndex>(&cleanup_settings);
            if next == result {
                break;
            }
            result = next;
        }
    }

    result
}

pub fn network_full_algebra_5(fixture: NetworkVertexFixture) -> Atom {
    network_schoonschip_substituted(network_vertex_substitution(fixture))
}

pub fn network_full_algebra_8(fixture: NetworkVertexFixture) -> Atom {
    network_schoonschip_substituted(network_vertex_substitution(fixture))
}

pub fn assert_no_network_internal_indices(result: &Atom) {
    assert_no_network_internal_indices_with_count(result, 8);
}

pub fn assert_no_network_internal_indices_with_count(result: &Atom, vertex_count: usize) {
    let residual = residual_network_internal_indices(result, vertex_count);
    assert!(
        residual.is_empty(),
        "internal indices remained: {residual:?}"
    );
}

pub fn residual_network_internal_indices(result: &Atom, vertex_count: usize) -> Vec<&'static str> {
    let (mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10, mu11) = symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    let indices: &[(&'static str, symbolica::atom::Symbol)] = match vertex_count {
        5 => &[("mu1", mu1), ("mu2", mu2), ("mu3", mu3), ("mu4", mu4)],
        6 => &[
            ("mu1", mu1),
            ("mu2", mu2),
            ("mu3", mu3),
            ("mu4", mu4),
            ("mu5", mu5),
            ("mu11", mu11),
        ],
        7 => &[
            ("mu1", mu1),
            ("mu2", mu2),
            ("mu3", mu3),
            ("mu4", mu4),
            ("mu5", mu5),
            ("mu6", mu6),
            ("mu10", mu10),
            ("mu11", mu11),
        ],
        8 => &[
            ("mu1", mu1),
            ("mu2", mu2),
            ("mu3", mu3),
            ("mu4", mu4),
            ("mu5", mu5),
            ("mu6", mu6),
            ("mu7", mu7),
            ("mu8", mu8),
            ("mu9", mu9),
            ("mu10", mu10),
            ("mu11", mu11),
        ],
        _ => panic!("unsupported network vertex count: {vertex_count}"),
    };

    indices
        .iter()
        .filter_map(|(name, index)| {
            result
                .replace(*index)
                .match_iter()
                .next()
                .is_some()
                .then_some(*name)
        })
        .collect()
}

fn activate_symbolica_license() {
    if let Ok(key) = std::env::var("SYMBOLICA_LICENSE") {
        let _ = LicenseManager::set_license_key(&key);
    }
}
