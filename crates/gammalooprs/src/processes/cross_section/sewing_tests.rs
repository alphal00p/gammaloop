//! Graph/UFO assignment is production code; the matrix adjoint and spin sums
//! below are independent component oracles, with no tensor conjugation API.

use std::collections::BTreeMap;

use linnet::half_edge::involution::HedgePair;

use spenso::{
    iterators::IteratableTensor,
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::symbolic::{ExplicitKey, TensorLibrary},
        parsing::ParseSettings,
    },
    structure::{TensorStructure, representation::Minkowski, slot::IsAbstractSlot},
    tensors::parametric::{ParamOrConcrete, ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
};

use crate::{
    graph::Graph,
    initialisation::test_initialise,
    model::{Model, UFOSymbol},
    numerator::{ParsingNet, aind::Aind},
    utils::{F, FUN_LIB, GS, load_generic_model},
};

type Matrix = [[Atom; 4]; 4];

fn product(left: &Matrix, right: &Matrix) -> Matrix {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..4).fold(Atom::Zero, |sum, k| sum + &left[i][k] * &right[k][j]))
    })
}

fn bar(matrix: &Matrix) -> Matrix {
    // Explicit Weyl gamma0 swaps 0<->2 and 1<->3. This is an independent
    // numerical-index implementation of gamma0 X^dagger gamma0.
    std::array::from_fn(|i| std::array::from_fn(|j| matrix[(j + 2) % 4][(i + 2) % 4].conj()))
}

fn negative(matrix: &Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| -&matrix[i][j]))
}

fn spin_sum(left: &Matrix, right: &Matrix) -> Atom {
    // G+ -> nu_tau tau+, M=5, m_tau=2, real on-shell momenta.
    // p_nu=(21/10,0,0,21/10), q_tau=(29/10,0,0,-21/10).
    // Construct pslash and qslash-m explicitly, without using the library's
    // gamma matrices or its spin-sum implementation.
    let p: Matrix = std::array::from_fn(|i| {
        std::array::from_fn(|j| match (i, j) {
            (1, 3) | (2, 0) => Atom::num((21, 5)),
            _ => Atom::Zero,
        })
    });
    let q: Matrix = std::array::from_fn(|i| {
        std::array::from_fn(|j| match (i, j) {
            (0, 0) | (1, 1) | (2, 2) | (3, 3) => Atom::num(-2),
            (0, 2) | (3, 1) => Atom::num(5),
            (1, 3) | (2, 0) => Atom::num((4, 5)),
            _ => Atom::Zero,
        })
    });
    let sewn = product(&product(&product(&p, left), &q), right);
    (0..4).fold(Atom::Zero, |sum, i| sum + &sewn[i][i])
}

fn vertex_matrix(
    graph: &Graph,
    vertex_name: &str,
    model: &Model,
    vector_component: usize,
    tau_yukawa: Atom,
    parameter_values: &[(&str, Atom)],
) -> Matrix {
    let (node, _, vertex) = graph
        .underlying
        .iter_nodes()
        .find(|(_, _, vertex)| {
            vertex.name.value == vertex_name
                || vertex
                    .vertex_rule
                    .as_ref()
                    .is_some_and(|rule| rule.name == vertex_name)
        })
        .unwrap();
    let mut index_rules = Vec::new();
    for hedge in graph.underlying.iter_crown(node) {
        let order = usize::from(graph.underlying[hedge].ufo_order.value);
        for slot in graph.underlying[hedge]
            .num_indices
            .spin_indices
            .vertex_indices
            .external_structure_iter()
        {
            index_rules.push(Replacement::new(
                slot.to_atom().to_pattern(),
                slot.rep().slot::<Aind, _>(Aind::Normal(order)).to_atom(),
            ));
        }
        for slot in graph.underlying[hedge]
            .num_indices
            .color_indices
            .vertex_indices
            .external_structure_iter()
        {
            index_rules.push(Replacement::new(
                slot.to_atom().to_pattern(),
                slot.rep()
                    .slot::<Aind, _>(Aind::Normal(3 + order))
                    .to_atom(),
            ));
        }
    }
    let mut expression = model
        .apply_coupling_replacement_rules(&vertex.num.value)
        .replace_multiple(&index_rules);
    // Fix the two quark colors to the same component: delta_00=1. This
    // color-singlet metric otherwise remains symbolic in the tensor library.
    for identity in [
        symbolica::parse!("spenso::g(spenso::cof(3,4),spenso::dind(spenso::cof(3,3)))"),
        symbolica::parse!("spenso::g(spenso::cof(3,3),spenso::dind(spenso::cof(3,4)))"),
    ] {
        expression = expression.replace(identity).with(Atom::one());
    }
    let ckm = model.get_parameter("CKM1x3");
    expression = expression
        .replace(Atom::from(ckm.name).to_pattern())
        .with(ckm.expression.as_ref().unwrap().to_pattern());
    for (name, value) in parameter_values {
        expression = expression
            .replace(Atom::from(UFOSymbol::from(*name)).to_pattern())
            .with(value.clone());
    }
    for (name, value) in [
        ("ytau", tau_yukawa),
        ("ee", Atom::one()),
        ("sw", Atom::one()),
        // The model's CKM1x3 = AWS lamWS^3 (rhoWS-i etaWS) is
        // genuinely complex while every Wolfenstein input remains real.
        ("AWS", Atom::one()),
        ("lamWS", Atom::one()),
        ("rhoWS", Atom::one()),
        ("etaWS", Atom::num(-2)),
    ] {
        expression = expression
            .replace(Atom::from(UFOSymbol::from(name)).to_pattern())
            .with(value);
    }

    tensor_matrix(expression, vector_component)
}

fn tensor_matrix(expression: Atom, vector_component: usize) -> Matrix {
    // Build the generic metric with Atom entries so its spatial signs
    // remain exact alongside the parametric gamma matrices and couplings.
    let mut library = spenso_hep_lib::hep_lib_atom::<Aind, F<f64>>();
    library.insert_generic(
        TensorLibrary::<ParamTensor<ExplicitKey<Aind>>, Aind>::id(Minkowski {}.into()),
        |key| {
            ParamOrConcrete::Param(
                TensorLibrary::<ParamTensor<ExplicitKey<Aind>>, Aind>::diag_unimodular_metric(key),
            )
        },
    );
    let mut network =
        ParsingNet::try_from_view(expression.as_view(), &library, &ParseSettings::default())
            .unwrap();
    network
        .execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB)
        .unwrap();
    let ExecutionResult::Val(tensor) = network.result_tensor(&library).unwrap() else {
        panic!("expected nonzero model vertex");
    };
    let mut tensor = tensor.into_owned();
    tensor.to_param();
    let tensor = tensor.try_into_parametric().unwrap();
    let slots = tensor.external_structure();
    let mut result: Matrix = std::array::from_fn(|_| std::array::from_fn(|_| Atom::Zero));
    for (indices, value) in tensor.iter_expanded() {
        let by_slot: BTreeMap<_, _> = slots
            .iter()
            .zip(indices.iter())
            .map(|(slot, index)| (slot.aind().to_string(), *index))
            .collect();
        // Select the requested W component (mu=1 for the transverse sewing
        // checks), and one equal pair of quark colors. Scalar/leptonic vertices
        // do not contain these slots.
        if by_slot
            .get("2")
            .is_some_and(|component| *component != vector_component)
            || by_slot.get("3").is_some_and(|component| *component != 0)
            || by_slot.get("4").is_some_and(|component| *component != 0)
        {
            continue;
        }
        result[by_slot["0"]][by_slot["1"]] = value.to_owned();
    }
    result
}

#[test]
fn generated_charged_scalar_forward_vertex_is_already_the_hermitian_partner() {
    test_initialise().unwrap();
    let model = load_generic_model("sm");
    let graphs = Graph::from_string(
        r#"
        digraph direct {
            ext [style=invis];
            ext -> amplitude [particle="G+"];
            amplitude -> ext [particle="vt"];
            amplitude -> ext [particle="ta+"];
        }
        digraph forward {
            ext [style=invis];
            ext -> left [particle="G+", is_cut=0];
            right -> ext [particle="G+", is_cut=0];
            left -> right [particle="vt"];
            left -> right [particle="ta+"];
        }
    "#,
        &model,
    )
    .unwrap();
    let amplitude = vertex_matrix(&graphs[0], "V_112", &model, 1, Atom::one(), &[]);
    let left = vertex_matrix(&graphs[1], "V_112", &model, 1, Atom::one(), &[]);
    let raw_right = vertex_matrix(&graphs[1], "V_147", &model, 1, Atom::one(), &[]);
    assert_eq!(left, amplitude);
    assert_eq!(negative(&raw_right), bar(&amplitude));
    assert_eq!(spin_sum(&amplitude, &bar(&amplitude)), Atom::num(21));
    assert_eq!(spin_sum(&left, &negative(&raw_right)), Atom::num(21));
    assert_eq!(spin_sum(&left, &raw_right), Atom::num(-21));
    assert_eq!(spin_sum(&left, &bar(&raw_right)), Atom::Zero);
}

#[test]
fn generated_complex_charged_current_preserves_the_ckm_norm() {
    test_initialise().unwrap();
    let model = load_generic_model("sm");
    let graphs = Graph::from_string(
        r#"
        digraph direct {
            ext [style=invis];
            ext -> amplitude [particle="W-"];
            amplitude -> ext [particle="b"];
            amplitude -> ext [particle="u~"];
        }
        digraph forward {
            ext [style=invis];
            ext -> left [particle="W-", is_cut=0];
            right -> ext [particle="W-", is_cut=0];
            left -> right [particle="b"];
            left -> right [particle="u~"];
        }
    "#,
        &model,
    )
    .unwrap();
    let amplitude = vertex_matrix(&graphs[0], "V_125", &model, 1, Atom::one(), &[]);
    let left = vertex_matrix(&graphs[1], "V_125", &model, 1, Atom::one(), &[]);
    let raw_right = vertex_matrix(&graphs[1], "V_95", &model, 1, Atom::one(), &[]);
    assert_eq!(left, amplitude);
    assert_eq!(negative(&raw_right), bar(&amplitude));
    // ee=sw=1 leaves the physical vertex factor 1/sqrt(2).
    // The real spin densities here are an algebraic helicity-sum oracle;
    // they are not a separate statement about physical SM quark masses.
    assert_eq!(spin_sum(&amplitude, &bar(&amplitude)), Atom::num((105, 2)));
    assert_eq!(spin_sum(&left, &negative(&raw_right)), Atom::num((105, 2)));
    assert_eq!(
        spin_sum(&left, &bar(&raw_right)),
        Atom::num((-63, 2)) + Atom::num(42) * Atom::i()
    );
}

#[test]
fn generated_sm_charged_ward_and_ghost_momentum_follow_the_action() {
    test_initialise().unwrap();
    // The approved canonical model carries the action-consistent Lorentz rules;
    // loading it rebuilds the vertex rules and their shared Lorentz references.
    let model = load_generic_model("sm");
    let graphs = Graph::from_string(
        r#"
            digraph goldstone {
                ext [style=invis];
                ext -> amplitude [particle="G+"];
                amplitude -> ext [particle="vt"];
                amplitude -> ext [particle="ta+"];
            }
            digraph charged_vector {
                ext [style=invis];
                ext -> amplitude [particle="W+"];
                amplitude -> ext [particle="vt"];
                amplitude -> ext [particle="ta+"];
            }
            digraph charged_ghost {
                ext [style=invis];
                ext -> amplitude [particle="ghWp"];
                ext -> amplitude [particle="ghWp~"];
                ext -> amplitude [particle="a"];
            }
            "#,
        &model,
    )
    .unwrap();

    // W+ -> nu_L tau+, M_W=5, m_tau=2, g=1 and y_tau=sqrt(2)/5.
    // These explicit Weyl spinors use p_nu=(21/10,0,0,21/10) and
    // p_tau=(29/10,0,0,-21/10). Their fixed neutrino chirality distinguishes PL from PR,
    // unlike the spin-summed sewing norm. No internal gauge propagator enters.
    // The Ward identity is homogeneous in both external spinors. Divide their
    // normalized representatives by sqrt(21/5) and sqrt(5), respectively, and
    // write y_tau=(2/5)/sqrt(2) so the exact check needs no radical identities.
    let inverse_sqrt_two = Atom::num(2).pow(Atom::num((-1, 2)));
    let yukawa = Atom::num((2, 5)) * &inverse_sqrt_two;
    let goldstone = vertex_matrix(&graphs[0], "V_112", &model, 0, yukawa.clone(), &[]);
    let vector = vertex_matrix(&graphs[1], "V_115", &model, 0, yukawa.clone(), &[]);
    for row in 0..4 {
        for column in 0..4 {
            // Independent Weyl entries of y_tau PR and i gamma0 PL/sqrt(2)
            // separate model/slot failures from the external-spinor contraction.
            let expected_goldstone = if row == column && row >= 2 {
                yukawa.clone()
            } else {
                Atom::Zero
            };
            let expected_vector = if matches!((row, column), (2, 0) | (3, 1)) {
                Atom::i() * &inverse_sqrt_two
            } else {
                Atom::Zero
            };
            assert_eq!(
                goldstone[row][column], expected_goldstone,
                "G+[{row},{column}]"
            );
            assert_eq!(vector[row][column], expected_vector, "W+_0[{row},{column}]");
        }
    }
    let ubar = [Atom::Zero, Atom::Zero, Atom::Zero, Atom::one()];
    let antitau = [Atom::Zero, Atom::num((-2, 5)), Atom::Zero, Atom::one()];
    let mut ward_residual = Atom::Zero;
    let mut goldstone_amplitude = Atom::Zero;
    for row in 0..4 {
        for column in 0..4 {
            goldstone_amplitude += &ubar[row] * &goldstone[row][column] * &antitau[column];
            ward_residual += &ubar[row]
                * (Atom::num(5) * &vector[row][column]
                    + Atom::i() * Atom::num(5) * &goldstone[row][column])
                * &antitau[column];
        }
    }
    assert_eq!(ward_residual, Atom::Zero);
    assert_eq!(goldstone_amplitude, yukawa);

    let graph = &graphs[2];
    let (node, _, vertex) = graph
        .underlying
        .iter_nodes()
        .find(|(_, _, vertex)| {
            vertex
                .vertex_rule
                .as_ref()
                .is_some_and(|rule| rule.name == "V_26")
        })
        .unwrap();
    let mut expression = model.apply_coupling_replacement_rules(&vertex.num.value);
    let mut antighost_edge = None;
    for hedge in graph.underlying.iter_crown(node) {
        match graph.underlying[hedge].ufo_order.value {
            1 => antighost_edge = Some(graph.underlying[&hedge]),
            2 => {
                for slot in graph.underlying[hedge]
                    .num_indices
                    .spin_indices
                    .vertex_indices
                    .external_structure_iter()
                {
                    expression = expression
                        .replace(slot.to_atom().to_pattern())
                        .with(GS.cind(0));
                }
            }
            _ => {}
        }
    }
    // The antighost flows into this vertex. With e=1 and pbar^0=3,
    // -barc partial D c, D=partial-ieA, requires -3i. Slot assignment
    // and incoming-momentum conversion above are the actual Graph pipeline.
    let ghost_value = expression
        .replace(Atom::from(UFOSymbol::from("ee")).to_pattern())
        .with(Atom::one())
        .replace(GS.emr_mom(antighost_edge.unwrap(), GS.cind(0)).to_pattern())
        .with(Atom::num(3));
    assert_eq!(ghost_value, Atom::num(-3) * Atom::i());
}

#[test]
fn generated_sm_virtual_vector_and_goldstone_exchange_matches_unitary_current() {
    test_initialise().unwrap();
    let model = load_generic_model("sm");
    let inverse_sqrt_two = Atom::num(2).pow(Atom::num((-1, 2)));
    let yukawa = Atom::num((2, 3)) * &inverse_sqrt_two;
    let bilinear = |bra: &[Atom; 4], matrix: &Matrix, ket: &[Atom; 4]| {
        (0..4).fold(Atom::Zero, |sum, row| {
            sum + (0..4).fold(Atom::Zero, |sum, column| {
                sum + &bra[row] * &matrix[row][column] * &ket[column]
            })
        })
    };

    for charged in [true, false] {
        let (vector, goldstone, fermion, mass_name) = if charged {
            ("W+", "G+", "vt", "MW")
        } else {
            ("Z", "G0", "ta-", "MZ")
        };
        let graphs = Graph::from_string(
            format!(
                r#"
                digraph vector_exchange {{
                    ext [style=invis];
                    left [name="left"];
                    right [name="right"];
                    ext -> left [particle="{fermion}"];
                    ext -> left [particle="ta+"];
                    left -> right [particle="{vector}", name="exchange"];
                    right -> ext [particle="{fermion}"];
                    right -> ext [particle="ta+"];
                }}
                digraph goldstone_exchange {{
                    ext [style=invis];
                    left [name="left"];
                    right [name="right"];
                    ext -> left [particle="{fermion}"];
                    ext -> left [particle="ta+"];
                    left -> right [particle="{goldstone}", name="exchange"];
                    right -> ext [particle="{fermion}"];
                    right -> ext [particle="ta+"];
                }}
                "#
            ),
            &model,
        )
        .unwrap();
        // In both sectors v=6 and m_tau=2. Charged: g=1, MW=3.
        // Neutral: g/cw=1, MZ=3. Each has sw=3/5, cw=4/5.
        let parameters = [
            ("ee", Atom::num(if charged { (3, 5) } else { (12, 25) })),
            ("sw", Atom::num((3, 5))),
            ("cw", Atom::num((4, 5))),
        ];
        // All spinors are exact on-shell Weyl representatives. Both channels
        // carry q=(5,0,0,0), away from the mediator pole q^2-M^2=16.
        // Charged incoming (nu,tau)=(21/10,+21/10),(29/10,-21/10)
        // along z; the outgoing pair is backscattered. Neutral tau energies
        // are 5/2 with momenta +/-3/2 along z, unchanged by the scattering.
        let (left_bra, left_ket, right_bra, right_ket) = if charged {
            (
                [Atom::Zero, Atom::one(), Atom::Zero, Atom::num((-2, 5))],
                [Atom::Zero, Atom::one(), Atom::Zero, Atom::Zero],
                [Atom::Zero, Atom::Zero, Atom::one(), Atom::Zero],
                [Atom::num((-2, 5)), Atom::Zero, Atom::one(), Atom::Zero],
            )
        } else {
            (
                [Atom::num(-1), Atom::Zero, Atom::num(2), Atom::Zero],
                [Atom::one(), Atom::Zero, Atom::num(2), Atom::Zero],
                [Atom::num(2), Atom::Zero, Atom::one(), Atom::Zero],
                [Atom::num(2), Atom::Zero, Atom::num(-1), Atom::Zero],
            )
        };
        let left: [Atom; 4] = std::array::from_fn(|component| {
            bilinear(
                &left_bra,
                &vertex_matrix(
                    &graphs[0],
                    "left",
                    &model,
                    component,
                    yukawa.clone(),
                    &parameters,
                ),
                &left_ket,
            )
        });
        let right: [Atom; 4] = std::array::from_fn(|component| {
            bilinear(
                &right_bra,
                &vertex_matrix(
                    &graphs[0],
                    "right",
                    &model,
                    component,
                    yukawa.clone(),
                    &parameters,
                ),
                &right_ket,
            )
        });
        let scalar_left = bilinear(
            &left_bra,
            &vertex_matrix(&graphs[1], "left", &model, 0, yukawa.clone(), &parameters),
            &left_ket,
        );
        let scalar_right = bilinear(
            &right_bra,
            &vertex_matrix(&graphs[1], "right", &model, 0, yukawa.clone(), &parameters),
            &right_ket,
        );
        // Opposite momentum directions fix the two local Ward signs. Keep
        // the raw inverse-process UFO vertices, with no extra adjoint.
        assert_eq!(
            Atom::num(5) * &left[0] - Atom::i() * Atom::num(3) * &scalar_left,
            Atom::Zero
        );
        assert_eq!(
            Atom::num(5) * &right[0] + Atom::i() * Atom::num(3) * &scalar_right,
            Atom::Zero
        );

        let (pair, _, edge) = graphs[0]
            .underlying
            .iter_edges()
            .find(|(_, _, edge)| edge.data.name.value == "exchange")
            .unwrap();
        let HedgePair::Paired { source, sink } = pair else {
            panic!("the vector exchange must be internal");
        };
        let mut numerator = edge.data.num.value.clone();
        for (order, hedge) in [source, sink].into_iter().enumerate() {
            for slot in graphs[0].underlying[hedge]
                .num_indices
                .spin_indices
                .edge_indices
                .external_structure_iter()
            {
                numerator = numerator
                    .replace(slot.to_atom().to_pattern())
                    .with(slot.rep().slot::<Aind, _>(Aind::Normal(order)).to_atom());
            }
        }
        let propagator = tensor_matrix(numerator, 0);
        let mass = edge
            .data
            .particle
            .mass_atom()
            .replace(Atom::from(UFOSymbol::from(mass_name)).to_pattern())
            .with(Atom::num(3));
        assert_eq!(mass, Atom::num(3));
        let denominator = Atom::num(25) - mass.pow(2);
        assert_eq!(denominator, Atom::num(16));
        let (_, _, scalar_edge) = graphs[1]
            .underlying
            .iter_edges()
            .find(|(_, _, edge)| edge.data.name.value == "exchange")
            .unwrap();
        assert_eq!(scalar_edge.data.num.value, Atom::i());
        assert_eq!(
            scalar_edge.data.particle.mass_atom(),
            edge.data.particle.mass_atom()
        );

        let mut covariant_vector = Atom::Zero;
        let mut unitary_vector = Atom::Zero;
        for mu in 0..4 {
            for nu in 0..4 {
                let metric = if mu != nu {
                    0
                } else if mu == 0 {
                    1
                } else {
                    -1
                };
                assert_eq!(propagator[mu][nu], -Atom::i() * Atom::num(metric));
                covariant_vector += &left[mu] * &propagator[mu][nu] * &right[nu] / &denominator;
                // Independent complete unitary propagator, with the same
                // denominator and q=(5,0,0,0); not a replacement graph rule.
                let qq_over_mass_squared = if mu == 0 && nu == 0 {
                    Atom::num((25, 9))
                } else {
                    Atom::Zero
                };
                unitary_vector +=
                    &left[mu] * Atom::i() * (qq_over_mass_squared - Atom::num(metric)) * &right[nu]
                        / &denominator;
            }
        }
        let scalar_exchange =
            scalar_left * &scalar_edge.data.num.value * scalar_right / denominator;
        let expected_vector = if charged { (1, 100) } else { (63, 1250) };
        let expected_scalar = if charged { (-1, 72) } else { (-25, 144) };
        let expected_total = if charged { (-7, 1800) } else { (-11089, 90000) };
        assert_eq!(
            covariant_vector,
            Atom::num(expected_vector) * Atom::i(),
            "{vector}"
        );
        assert_eq!(
            scalar_exchange,
            Atom::num(expected_scalar) * Atom::i(),
            "{goldstone}"
        );
        assert_eq!(
            unitary_vector,
            Atom::num(expected_total) * Atom::i(),
            "{vector}"
        );
        assert_eq!(&covariant_vector + &scalar_exchange, unitary_vector);
        // The previous hybrid counted the longitudinal/Goldstone contribution
        // twice. Its complete exchange must fail this same absolute oracle.
        assert_ne!(&unitary_vector + scalar_exchange, unitary_vector);
    }
}
