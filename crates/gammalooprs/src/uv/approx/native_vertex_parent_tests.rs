//! Private proof boundary for actual logarithmic children inside OS parents.
//! Public mixed-prescription admission is deliberately unchanged.

use super::*;
use crate::{
    initialisation::test_initialise,
    numerator::AppliedFeynmanRule,
    uv::{
        approx::{integrated, integrated::native_mass_tests::QcdMassReference},
        hedge_poset::OwnedForestNode,
        orchestrator::validate_local_counterterm_schemes,
    },
};
use idenso::{
    color::{CS, ColorSimplifySettings},
    dirac::GammaSimplifySettings,
    representations::Bispinor,
};
use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
use spenso::structure::{
    TensorStructure,
    abstract_index::{AIND_SYMBOLS, AbstractIndex},
    partial::PartialStructure,
    representation::{LibraryRep, Minkowski, RepName},
    slot::{IsAbstractSlot, Slot, SlotMatcher},
};

#[test]
fn native_logarithmic_vertex_prefixes_preserve_parent_mass_projection() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_parent_mass_projection(["1", "0", "0", "0"])
}

#[test]
fn native_logarithmic_vertex_prefixes_preserve_moving_parent_mass_projection() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_parent_mass_projection(["5/4", "1/4", "1/2", "1/2"])
}

#[test]
fn native_parent_scalar_certificate_preserves_factored_cancellation() -> Result<()> {
    test_initialise()?;
    let dimension = Atom::var(GS.dim);
    let imaginary = Atom::i();
    // D01 coefficient row 57, with its five complete denominator functions
    // replaced by independent scalar aliases. Keep their inverses and the
    // original four complex factors; do not clear a common denominator.
    let inverses = (0..5).map(|index| {
        Atom::var(symbol!(
            &format!("native_parent_scalar_regression::s{index}"); Scalar
        ))
        .pow(-1)
    });
    let common = Atom::var(symbol!("UFO::GC_10")).pow(2)
        * Atom::var(symbol!("UFO::GC_11")).pow(2)
        * Atom::mul_many(inverses);
    let scalar = Atom::add_many([
        (Atom::num(-2) * &imaginary * &dimension + Atom::num(4) * &imaginary)
            * Atom::num(4)
            * &common,
        (-&imaginary * &dimension + Atom::num(2) * &imaginary) * Atom::num(8) * &common,
        (Atom::num(2) * &imaginary * &dimension - Atom::num(4) * &imaginary)
            * Atom::num(4)
            * &common,
        (&imaginary * &dimension - Atom::num(2) * &imaginary) * Atom::num(8) * &common,
    ]);
    let metric = function!(
        ETS.metric,
        idenso::bis!(
            4,
            Atom::var(symbol!("native_parent_scalar_regression::row"))
        ),
        idenso::bis!(
            4,
            Atom::var(symbol!("native_parent_scalar_regression::column"))
        )
    );
    QcdMassReference::certify_parent_tensor_zero(
        "captured scalar coefficient on an open metric",
        &AppliedFeynmanRule::new(&scalar * &metric)?.expr,
    )?;
    assert!(
        QcdMassReference::certify_parent_tensor_zero(
            "unit perturbation of the same scalar coefficient",
            &AppliedFeynmanRule::new((&scalar + Atom::one()) * &metric)?.expr,
        )
        .is_err(),
        "the private certificate must reject a nonzero metric contribution"
    );
    Ok(())
}

#[test]
fn native_parent_scalar_certificate_accepts_ranked_zero() -> Result<()> {
    test_initialise()?;
    let adjoint = idenso::representations::ColorAdjoint {}.new_rep(8);
    let [a, b, c, d] = ["a", "b", "c", "d"].map(|name| {
        adjoint.to_symbolic([Atom::var(symbol!(&format!(
            "native_parent_ranked_zero::{name}"
        )))])
    });
    let metric = function!(ETS.metric, &c, &d);
    let product = function!(CS.f, &a, &b, &c) * function!(CS.f, &a, &b, &d);
    // The final SU(3) color reduction gives a zero with two retained ports.
    // Its coefficient need not have a scalar interface to certify this zero.
    QcdMassReference::certify_parent_tensor_zero(
        "open adjoint Casimir identity",
        &AppliedFeynmanRule::new(&product - Atom::num(3) * &metric)?.expr,
    )?;
    assert!(
        QcdMassReference::certify_parent_tensor_zero(
            "nonzero adjoint metric remainder",
            &AppliedFeynmanRule::new(product - Atom::num(2) * metric)?.expr,
        )
        .is_err()
    );
    Ok(())
}

#[test]
fn native_d09_laurent_grid_preserves_factors_and_rejects_unsupported_domains() -> Result<()> {
    test_initialise()?;
    let variables =
        ["A", "H", "G", "t"].map(|name| symbol!(&format!("native_d09_grid_test::{name}")));
    let [a, h, coupling, t] = variables.map(Atom::var);
    let b = (&t - t.pow(-1)) / Atom::num(2);
    let c = (&t + t.pow(-1)) / Atom::num(2);
    let hyperbola = c.pow(2) - b.pow(2) - Atom::one();
    let spectator = (&a + &h).pow(3) * coupling.pow(4);
    let original = &spectator * &hyperbola;
    let retained = original.clone();
    assert_eq!(
        QcdMassReference::certify_d09_phi_laurent_grid(&original, variables)?,
        400
    );
    assert_eq!(original, retained);
    assert!(
        QcdMassReference::certify_d09_phi_laurent_grid(&(&original + Atom::one()), variables)
            .is_err()
    );
    for unsupported in [
        function!(symbol!("native_d09_grid_test::opaque"), &a),
        (&t + Atom::one()).pow(-1),
        a.pow(-1),
        a.pow((1, 2)),
        Atom::num(0.125f64),
    ] {
        assert!(
            QcdMassReference::certify_d09_phi_laurent_grid(&(&unsupported * &hyperbola), variables)
                .is_err(),
            "unsupported chart factor was accepted: {unsupported}"
        );
    }
    let too_large = (&a + &h).pow(16) * coupling.pow(4) * hyperbola;
    assert!(QcdMassReference::certify_d09_phi_laurent_grid(&too_large, variables).is_err());
    Ok(())
}

impl QcdMassReference {
    fn certify_parent_tensor_zero(
        label: &str,
        difference: &SymbolicTensor<IndexedStructure>,
    ) -> Result<()> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        debug_tags!(#uv, #integrated, #inspect;
            stage = "native_vertex_parent_exact_input", label,
            file.difference = %cooking.uncook(difference.expression().as_view()).to_canonical_string(),
            "Retained factorized input to the parent proof certificate"
        );
        let reduced = integrated::simplify(difference)
            .wrap_err_with(|| format!("{label}: parent proof algebra"))?;
        let value = reduced.simplify_algebra(&idenso::tensor::AlgebraSettings {
            // The source has native SU(3) ports and the independent
            // vertex reference uses CF=4/3, CA=3, as in the child test.
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            ..Default::default()
        })?;
        // The coefficient engine allows partial words inside the complete tensor.
        // Move the reduced carrier into that boundary, retaining layout and proofs.
        let value: SymbolicTensor<PartialStructure> = value.into();
        // Native collection keeps each complete tensor word opaque. Normalize
        // only its scalar coefficient, as in the independent bare-source proof;
        // complete denominator functions and powers retain their owners.
        for (row, (word, coefficient)) in value
            .coefficient_list(spenso::shadowing::TensorCollectFilter::<0>::Tensors)?
            .into_iter()
            .enumerate()
        {
            // A certified zero can retain logical tensor ports. Preserve the
            // native certificate's zero-factor rule before testing scalar rank.
            if word.expression().is_zero() || coefficient.expression().is_zero() {
                continue;
            }
            eyre::ensure!(
                coefficient.is_scalar(),
                "{label}: nonzero coefficient {row} retains tensor ports"
            );
            let normalized = coefficient.expression().collect_compact_factors();
            if !normalized.is_zero() {
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_vertex_parent_exact_coefficient", label, row,
                    file.coefficient = %coefficient.expression().to_canonical_string(),
                    file.normalized = %normalized.to_canonical_string(),
                    "Unresolved scalar coefficient retains its factorized source"
                );
            }
            eyre::ensure!(
                normalized.is_zero(),
                "{label}: exact scalar coefficient {row} is not zero"
            );
        }
        Ok(())
    }

    /// Reversible external-frame certificate only: no internal-d completeness.
    pub(in crate::uv) fn reduce_parent_external_frame(
        label: &str,
        difference: &SymbolicTensor<IndexedStructure>,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let source_tensor = difference;
        let source_expression = cooking.uncook(source_tensor.expression().as_view());
        let difference = &source_expression;
        // Only the certified component OS operators introduce this fixed frame.
        // Embed its four Gram relations into formal d, without assuming that
        // these four vectors span the internal Lorentz space.
        let d = Atom::var(GS.dim);
        let rep = Minkowski {}.to_symbolic([&d]);
        let frame = symbol!(
            "native_vertex_parent::frame",
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        );
        let names = (0..4).map(|i| function!(frame, i)).collect::<Vec<_>>();
        let vectors = (0..4)
            .map(|i| function!(frame, i, &rep))
            .collect::<Vec<_>>();
        let selectors = (0..4)
            .map(|i| function!(ETS.delta, function!(AIND_SYMBOLS.cind, i), &rep))
            .collect::<Vec<_>>();
        let mut kinematics = feynkit_kinematics::Kinematics::in_dimension(&d)?;
        for i in 0..4 {
            for j in 0..4 {
                let gram = Atom::num(if i != j {
                    0
                } else if i == 0 {
                    1
                } else {
                    -1
                });
                assert_eq!(
                    function!(SPENSO_TAG.dot, &selectors[i], &selectors[j]),
                    gram
                );
                kinematics = kinematics.with_scalar_product(&names[i], &names[j], gram)?;
            }
        }
        // Embed the raw class before any gamma-zero-sensitive reduction.
        // A literal inverse below certifies the whole factorized source;
        // no completeness relation for the internal d-space is assumed.
        let momentum_leaves = |value: &Atom| {
            let mut leaves = BTreeSet::new();
            value.visitor(&mut |view| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    assert!(momentum.get_nargs() >= 2);
                    let AtomView::Fun(port) = momentum.get(momentum.get_nargs() - 1) else {
                        panic!("the frame certificate does not convert fixed internal momentum components")
                    };
                    assert_eq!(port.get_symbol(), LibraryRep::from(Minkowski {}).symbol());
                    assert!(matches!(port.get_nargs(), 1 | 2));
                    assert_eq!(port.get(0), d.as_view());
                    leaves.insert(view.to_owned());
                }
                true
            });
            leaves
        };
        assert!(!difference.contains_symbol(frame));
        let momenta = momentum_leaves(difference);
        let embedded = difference.replace_map(|view, _, output| {
            if let AtomView::Fun(value) = view {
                if value.get_symbol() == AGS.gamma0 {
                    assert_eq!(value.get_nargs(), 2);
                    **output = function!(AGS.gamma, value.get(0), value.get(1), &vectors[0]);
                } else if value.get_symbol() == ETS.delta {
                    assert_eq!(value.get_nargs(), 2);
                    let AtomView::Fun(coordinate) = value.get(0) else {
                        panic!("the external frame requires a native fixed coordinate")
                    };
                    assert_eq!(coordinate.get_symbol(), AIND_SYMBOLS.cind);
                    assert_eq!(coordinate.get_nargs(), 1);
                    let i =
                        usize::try_from(coordinate.get(0)).expect("literal external coordinate");
                    assert!(i < 4);
                    let AtomView::Fun(port) = value.get(1) else {
                        panic!("the external frame requires an explicit Minkowski port")
                    };
                    assert_eq!(port.get_symbol(), LibraryRep::from(Minkowski {}).symbol());
                    assert!(matches!(port.get_nargs(), 1 | 2));
                    assert_eq!(port.get(0), d.as_view());
                    **output = function!(frame, i, value.get(1));
                }
            }
        });
        assert_eq!(momentum_leaves(&embedded), momenta);
        let reverse = |value: &Atom| {
            value.replace_map(|view, _, output| {
                if let AtomView::Fun(value) = view {
                    if value.get_symbol() == AGS.gamma
                        && value.get_nargs() == 3
                        && value.get(2) == vectors[0].as_view()
                    {
                        **output = function!(AGS.gamma0, value.get(0), value.get(1));
                    } else if value.get_symbol() == frame {
                        assert_eq!(value.get_nargs(), 2);
                        **output = function!(
                            ETS.delta,
                            function!(AIND_SYMBOLS.cind, value.get(0)),
                            value.get(1)
                        );
                    }
                }
            })
        };
        assert_eq!(
            reverse(&embedded),
            *difference,
            "literal source reconstruction precedes all formal-d Clifford algebra"
        );
        let embedded_tensor =
            source_tensor.with_checked_expression(cooking.try_cook(embedded.as_view())?)?;
        assert_eq!(
            source_tensor.structure().logical_slots(),
            embedded_tensor.structure().logical_slots()
        );
        // All public ports are explicit. Commutative rebuilding may change
        // their discovery layout, but never their named tensor interface.
        assert_eq!(
            source_tensor.structure().canonical().external_structure(),
            embedded_tensor.structure().canonical().external_structure()
        );
        let algebra = idenso::tensor::AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            epsilon: true,
            ..Default::default()
        };
        let ordinary_tensor = embedded_tensor.simplify_algebra(&algebra)?.to_dots()?;
        let ordinary = cooking.uncook(ordinary_tensor.expression().as_view());
        // Chain traversal may transpose the embedded temporal matrix. Only
        // E0 inherits the native gamma0 endpoint symmetry; no generic gamma
        // or internal momentum is made symmetric by this certificate.
        let mut temporal_transposes = 0;
        let ordered = ordinary.replace_map(|view, _, output| {
            if let AtomView::Fun(value) = view
                && value.get_symbol() == AGS.gamma
                && value.get_nargs() == 3
                && value.get(2) == vectors[0].as_view()
            {
                let native = function!(AGS.gamma0, value.get(0), value.get(1));
                let AtomView::Fun(temporal) = native.as_view() else {
                    unreachable!()
                };
                let normalized =
                    function!(AGS.gamma, temporal.get(0), temporal.get(1), &vectors[0]);
                if normalized.as_view() != view {
                    temporal_transposes += 1;
                    **output = normalized;
                }
            }
        });
        assert_eq!(reverse(&ordered), reverse(&ordinary));
        assert_eq!(momentum_leaves(&ordered), momentum_leaves(&ordinary));
        let completed_tensor = ordinary_tensor
            .with_checked_expression(cooking.try_cook(ordered.as_view())?)?
            .simplify_algebra(&algebra)?
            .to_dots()?;
        let completed = cooking.uncook(completed_tensor.expression().as_view());
        let constrained = kinematics.apply(&completed);
        debug_tags!(#uv, #integrated, #inspect;
            stage = "native_vertex_parent_external_frame_certificate", label,
            file.unembedded = %difference.to_canonical_string(),
            file.embedded = %embedded.to_canonical_string(),
            file.constrained = %constrained.to_canonical_string(),
            unchanged_momentum_leaves = momenta.len(), temporal_transposes,
            "Reversible external frame, native temporal symmetry and Gram relations retain internal d"
        );
        Ok(completed_tensor.with_checked_expression(cooking.try_cook(constrained.as_view())?)?)
    }

    /// Prove a factorized D09 chart identity by a full exact unisolvent grid.
    /// The caller retains its actual-row, energy-chart and unit witnesses.
    fn certify_d09_phi_laurent_grid(
        expression: &Atom,
        variables: [symbolica::atom::Symbol; 4],
    ) -> Result<usize> {
        use symbolica::coefficient::CoefficientView;

        #[derive(Clone, Copy)]
        struct Bounds {
            lower: [i32; 4],
            upper: [i32; 4],
        }
        fn bounds(
            atom: AtomView<'_>,
            variables: &[symbolica::atom::Symbol; 4],
            depth: usize,
            nodes: &mut usize,
        ) -> Result<Bounds> {
            *nodes += 1;
            eyre::ensure!(
                depth <= 64 && *nodes <= 16384,
                "D09 Laurent syntax budget exhausted"
            );
            let zero = Bounds {
                lower: [0; 4],
                upper: [0; 4],
            };
            let result = match atom {
                AtomView::Num(number) => {
                    eyre::ensure!(
                        matches!(
                            number.get_coeff_view(),
                            CoefficientView::Natural(..) | CoefficientView::Large(..)
                        ),
                        "D09 Laurent coefficients must be exact QQ(i)"
                    );
                    zero
                }
                AtomView::Var(variable) => {
                    let index = variables
                        .iter()
                        .position(|symbol| *symbol == variable.get_symbol())
                        .ok_or_else(|| eyre!("D09 Laurent chart contains an unbound variable"))?;
                    let mut result = zero;
                    result.lower[index] = 1;
                    result.upper[index] = 1;
                    result
                }
                AtomView::Fun(_) => return Err(eyre!("D09 Laurent chart refuses functions")),
                AtomView::Add(sum) => {
                    let mut result = None::<Bounds>;
                    for term in sum {
                        let next = bounds(term, variables, depth + 1, nodes)?;
                        result = Some(match result {
                            None => next,
                            Some(previous) => Bounds {
                                lower: std::array::from_fn(|i| {
                                    previous.lower[i].min(next.lower[i])
                                }),
                                upper: std::array::from_fn(|i| {
                                    previous.upper[i].max(next.upper[i])
                                }),
                            },
                        });
                    }
                    result.unwrap_or(zero)
                }
                AtomView::Mul(product) => {
                    let mut result = zero;
                    for factor in product {
                        let next = bounds(factor, variables, depth + 1, nodes)?;
                        for i in 0..4 {
                            result.lower[i] = result.lower[i]
                                .checked_add(next.lower[i])
                                .ok_or_else(|| eyre!("D09 Laurent exponent overflow"))?;
                            result.upper[i] = result.upper[i]
                                .checked_add(next.upper[i])
                                .ok_or_else(|| eyre!("D09 Laurent exponent overflow"))?;
                        }
                    }
                    result
                }
                AtomView::Pow(power) => {
                    let exponent = i64::try_from(power.get_exp())
                        .map_err(|_| eyre!("D09 Laurent chart refuses fractional powers"))?;
                    eyre::ensure!(
                        exponent.unsigned_abs() <= 32,
                        "D09 Laurent power budget exhausted"
                    );
                    if exponent < 0 {
                        eyre::ensure!(
                            matches!(power.get_base(), AtomView::Var(variable) if variable.get_symbol() == variables[3]),
                            "D09 Laurent inverse must be a monomial in t"
                        );
                        let mut result = zero;
                        result.lower[3] = exponent as i32;
                        result.upper[3] = exponent as i32;
                        result
                    } else {
                        let mut result = bounds(power.get_base(), variables, depth + 1, nodes)?;
                        for i in 0..4 {
                            result.lower[i] = result.lower[i]
                                .checked_mul(exponent as i32)
                                .ok_or_else(|| eyre!("D09 Laurent exponent overflow"))?;
                            result.upper[i] = result.upper[i]
                                .checked_mul(exponent as i32)
                                .ok_or_else(|| eyre!("D09 Laurent exponent overflow"))?;
                        }
                        result
                    }
                }
            };
            Ok(result)
        }
        fn value(
            atom: AtomView<'_>,
            variables: &[symbolica::atom::Symbol; 4],
            point: &[i64; 4],
        ) -> Result<Atom> {
            let result = match atom {
                AtomView::Num(_) => atom.to_owned(),
                AtomView::Var(variable) => {
                    let index = variables
                        .iter()
                        .position(|symbol| *symbol == variable.get_symbol())
                        .ok_or_else(|| eyre!("D09 exact grid variable is unbound"))?;
                    Atom::num(point[index])
                }
                AtomView::Add(sum) => Atom::add_many(
                    sum.iter()
                        .map(|term| value(term, variables, point))
                        .collect::<Result<Vec<_>>>()?,
                ),
                AtomView::Mul(product) => Atom::mul_many(
                    product
                        .iter()
                        .map(|factor| value(factor, variables, point))
                        .collect::<Result<Vec<_>>>()?,
                ),
                AtomView::Pow(power) => value(power.get_base(), variables, point)?.pow(
                    i64::try_from(power.get_exp())
                        .map_err(|_| eyre!("D09 exact grid exponent is not integral"))?,
                ),
                AtomView::Fun(_) => return Err(eyre!("D09 exact grid encountered a function")),
            };
            eyre::ensure!(
                matches!(result.as_view(), AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Natural(..) | CoefficientView::Large(..))),
                "D09 exact grid arithmetic left QQ(i)"
            );
            eyre::ensure!(
                result.as_view().get_byte_size() <= 65536,
                "D09 exact grid rational size exhausted"
            );
            Ok(result)
        }
        eyre::ensure!(
            variables.into_iter().collect::<BTreeSet<_>>().len() == 4,
            "D09 Laurent chart variables must be distinct"
        );
        eyre::ensure!(
            expression.as_view().get_byte_size() <= 32768,
            "D09 factorized Laurent input exceeds 32KiB"
        );
        let mut nodes = 0;
        let degree = bounds(expression.as_view(), &variables, 0, &mut nodes)?;
        eyre::ensure!(
            (0..3).all(|i| degree.lower[i] >= 0),
            "D09 Laurent chart has a non-t inverse"
        );
        eyre::ensure!(
            (0..4)
                .all(|i| degree.lower[i].unsigned_abs() <= 64
                    && degree.upper[i].unsigned_abs() <= 64),
            "D09 Laurent exponent range exceeds 64"
        );
        // Include exponent zero in EVERY axis, so this same full box also
        // certifies the independently evaluated original input plus one.
        let lower = degree.lower.map(|n| n.min(0));
        let upper = degree.upper.map(|n| n.max(0));
        let counts: [usize; 4] = std::array::from_fn(|i| (upper[i] - lower[i] + 1) as usize);
        let grid_nodes = counts
            .iter()
            .try_fold(1usize, |n, count| n.checked_mul(*count))
            .ok_or_else(|| eyre!("D09 exact grid size overflow"))?;
        eyre::ensure!(grid_nodes <= 4096, "D09 exact grid exceeds 4096 nodes");
        let wrong = expression + Atom::one();
        let wrong_degree = bounds(wrong.as_view(), &variables, 0, &mut 0)?;
        eyre::ensure!(
            (0..4).all(|i| wrong_degree.lower[i] >= lower[i] && wrong_degree.upper[i] <= upper[i]),
            "D09 original +1 control escaped its degree box"
        );
        let mut visited = 0;
        for a in 1..=counts[0] {
            for h in 1..=counts[1] {
                for coupling in 1..=counts[2] {
                    for t in 2..=counts[3] + 1 {
                        let point = [a as i64, h as i64, coupling as i64, t as i64];
                        eyre::ensure!(
                            value(expression.as_view(), &variables, &point)?.is_zero(),
                            "D09 factorized Laurent identity failed at exact point {point:?}"
                        );
                        eyre::ensure!(
                            value(wrong.as_view(), &variables, &point)?.is_one(),
                            "D09 original +1 control is not exactly one at {point:?}"
                        );
                        visited += 1;
                    }
                }
            }
        }
        eyre::ensure!(visited == grid_nodes);
        // Multiplying only conceptually by t^(-lower[3]) gives an ordinary
        // polynomial. Its full distinct-node tensor grid is unisolvent.
        debug_tags!(#uv, #local, #four_d, #cff, #summary;
            stage = "native_d09_phi_exact_laurent_grid",
            file.variables = ?variables, file.lower_exponents = ?degree.lower,
            file.upper_exponents = ?degree.upper, file.padded_lower = ?lower,
            file.padded_upper = ?upper, file.grid_counts = ?counts,
            grid_nodes, syntax_nodes = nodes, exact_QQ_i = true,
            original_zero_at_every_node = true, original_plus_one_exact_one = true,
            direct_phi_polynomial_conversion = false, factorized_input_preserved = true,
            "A complete exact grid proves the factorized chart identity"
        );
        Ok(grid_nodes)
    }

    /// Diagnostic copies of CFF scalar weights only. Phi bodies never enter
    /// the polynomial or rational owners used here.
    fn normalized_d09_phi_coefficient(coefficient: &Atom) -> Result<Atom> {
        use crate::numerator::symbolica_ext::ScalarRationalSize;
        use spenso::shadowing::Collectable;

        fn tiny_polynomial(atom: AtomView<'_>) -> Option<(usize, u32)> {
            let (terms, degree) = match atom {
                AtomView::Num(_) => (1, 0),
                AtomView::Var(_) | AtomView::Fun(_) => (1, 1),
                AtomView::Add(sum) => {
                    let mut terms = 0usize;
                    let mut degree = 0;
                    for term in sum {
                        let (count, power) = tiny_polynomial(term)?;
                        terms = terms.checked_add(count)?;
                        degree = degree.max(power);
                        if terms > 128 || degree > 4 {
                            return None;
                        }
                    }
                    (terms, degree)
                }
                AtomView::Mul(product) => {
                    let mut terms = 1usize;
                    let mut degree = 0u32;
                    for factor in product {
                        let (count, power) = tiny_polynomial(factor)?;
                        terms = terms.checked_mul(count)?;
                        degree = degree.checked_add(power)?;
                        if terms > 128 || degree > 4 {
                            return None;
                        }
                    }
                    (terms, degree)
                }
                AtomView::Pow(power) => {
                    let exponent = u32::try_from(i64::try_from(power.get_exp()).ok()?).ok()?;
                    if exponent > 4 {
                        return None;
                    }
                    let (terms, degree) = tiny_polynomial(power.get_base())?;
                    (terms.checked_pow(exponent)?, degree.checked_mul(exponent)?)
                }
            };
            (terms <= 128 && degree <= 4).then_some((terms, degree))
        }

        eyre::ensure!(
            coefficient.as_view().get_byte_size() <= 65536,
            "D09 scalar coefficient exceeds its 64KiB byte scope"
        );
        for head in [
            symbolica::atom::Symbol::IF,
            GS.theta,
            GS.orientation_delta,
            crate::cff::expression::OrientationID::symbol(),
            symbol!("gammalooprs::uv::numerator_family"),
        ] {
            eyre::ensure!(
                !coefficient.contains_symbol(head),
                "D09 coefficient certificate refuses guards or numerator families"
            );
        }
        let mut exact = true;
        coefficient.visitor(&mut |view| {
            if let AtomView::Num(number) = view {
                exact &= matches!(
                    number.get_coeff_view(),
                    symbolica::coefficient::CoefficientView::Natural(..)
                        | symbolica::coefficient::CoefficientView::Large(..)
                );
            }
            true
        });
        eyre::ensure!(exact, "D09 scalar coefficient is not exact");
        let mut refused = false;
        let mut scalar =
            GS.lower_esurfaces(coefficient)
                .replace_map_bottom_up(|view, _, output| {
                    if let AtomView::Fun(call) = view {
                        if call.get_symbol() == GS.on_shell_energy {
                            if call.get_nargs() != 2 {
                                refused = true;
                                return;
                            }
                            **output = call.get(1).to_owned().pow((1, 2));
                        } else if call.get_symbol() == GS.tree_denom_wrapper {
                            if call.get_nargs() != 1 {
                                refused = true;
                                return;
                            }
                            **output = call.get(0).to_owned();
                        }
                    }
                });
        eyre::ensure!(
            !refused,
            "D09 coefficient has a malformed evaluator function"
        );
        scalar = scalar.replace_map_bottom_up(|view, _, output| {
            if let AtomView::Pow(power) = view
                && i64::try_from(power.get_exp()).is_err()
            {
                let supported = (-15i64..16)
                    .step_by(2)
                    .any(|n| power.get_exp() == Atom::num((n, 2)).as_view());
                if !supported
                    || power.get_base().get_byte_size() > 8192
                    || tiny_polynomial(power.get_base()).is_none()
                {
                    refused = true;
                    return;
                }
                // Only an existing half-power's tiny complete base expands.
                // No graph numerator or Phi body is supplied to this method.
                let base = power.get_base().expand();
                if base.as_view().get_byte_size() > 32768 {
                    refused = true;
                    return;
                }
                **output = base.pow(power.get_exp());
            }
        });
        eyre::ensure!(
            !refused,
            "D09 coefficient half-power invariant exceeds its scope"
        );
        let mut numeric_fixed = false;
        for _ in 0..16 {
            let next = scalar.replace_map_bottom_up(|view, _, output| {
                if matches!(view, AtomView::Add(_) | AtomView::Mul(_)) {
                    **output = view.collect_numeric_content();
                }
            });
            if next == scalar {
                numeric_fixed = true;
                break;
            }
            scalar = next;
        }
        eyre::ensure!(
            numeric_fixed,
            "D09 numeric-content coefficient factoring did not settle"
        );
        let opaque = symbol!("native_d09_coefficient_certificate::opaque_power"; Scalar);
        for _ in 0..5 {
            eyre::ensure!(
                scalar.as_view().get_byte_size() <= 65536 && !scalar.contains_symbol(opaque),
                "D09 scalar coefficient scope exhausted"
            );
            // Existing complete inverse bases and fractional powers are exact
            // scalar leaves, as in the existing scalar-pole cancellation owner.
            let protected = scalar.replace_map(|view, _, output| {
                if let AtomView::Pow(power) = view {
                    if i64::try_from(power.get_exp()).is_err() {
                        **output = opaque.call_args([Atom::one(), view.to_owned()]);
                    } else if i64::try_from(power.get_exp()).is_ok_and(|n| n < 0)
                        && matches!(power.get_base(), AtomView::Add(_))
                    {
                        **output = opaque
                            .call_args([Atom::Zero, power.get_base().to_owned()])
                            .pow(power.get_exp());
                    }
                }
            });
            let restore = |atom: &Atom| {
                atom.replace_map_bottom_up(|view, _, output| {
                    if let AtomView::Fun(call) = view
                        && call.get_symbol() == opaque
                    {
                        **output = call.get(1).to_owned();
                    }
                })
            };
            eyre::ensure!(
                restore(&protected) == scalar,
                "D09 scalar coefficient opaque-power roundtrip changed its input"
            );
            eyre::ensure!(
                ScalarRationalSize::estimate(protected.as_view()).is_some(),
                "D09 scalar coefficient refused the unchanged rational preflight budget"
            );
            let next = restore(&protected.together());
            eyre::ensure!(
                next.as_view().get_byte_size() <= 65536,
                "D09 restored scalar coefficient exceeds its byte scope"
            );
            if next.is_zero() || next == scalar {
                return Ok(next);
            }
            scalar = next;
        }
        Err(eyre!(
            "D09 scalar coefficient did not settle within five rational rounds"
        ))
    }

    /// Actual D09 rational and factorized-CFF proof boundaries; no evaluator.
    pub(in crate::uv) fn certify_d09_rational_finalizer_mass_projection(
        &self,
        diagram: &feynkit_graph::FeynmanDiagram,
        ctx: &UVCtx<'_>,
        given: &OwnedForestNode,
        prefix: &Full4dCts,
        actual: &(OwnedForestNode, Local4dCts, IntegratedCts),
    ) -> Result<()> {
        use crate::uv::approx::final_integrand::{
            FinalIntegrandBuilder,
            projected_assembly_certificate::{ProjectedAssemblyCertificate, ProjectedAssemblySkip},
        };
        use std::sync::Mutex;

        let (parent, local, integrated) = actual;
        let graph = ctx.graph;
        eyre::ensure!(diagram.name() == "D09" && graph.name == "D09");
        eyre::ensure!(parent.renormalization_scheme() == ApproximationType::OS);
        eyre::ensure!(
            ctx.settings.on_shell_velocity()? == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]
        );
        eyre::ensure!(local.recursive_completion().is_empty());
        // Integrated::run borrows local immutably. Observe every original
        // signed analytic sector before any physical diagnostic copy below.
        // Keep serialization inside the existing filtered debug event.
        {
            let raw_cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            for (sector_ordinal, sector) in local.active_sectors().iter().enumerate() {
                let vacuum = sector.integration_atom();
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_d09_actual_vacuum_sector_input", source = %diagram.name(),
                    given = %given.subgraph().string_label(), sector_ordinal,
                    file.ownership = %serde_json::to_string(&serde_json::json!({
                        "active_components": &sector.active_components,
                        "frozen_lmbs": &sector.frozen_lmbs,
                        "inherited_loop_edges": &sector.inherited_loop_edges,
                        "component_provenance": &sector.provenance_components,
                        "branch_provenance": &sector.provenance_branches,
                        "has_soft_ancestry": sector.has_soft_ancestry,
                        "has_on_shell_ancestry": sector.has_on_shell_ancestry,
                        "explicit_scheme_integrand": sector.scheme_integrand.is_some(),
                        "active_loop_edges": sector.active_components.iter()
                            .flat_map(|(_, _, frame)| frame.loop_edges.iter().copied())
                            .collect::<BTreeSet<_>>(),
                        "logical_slots": vacuum.structure().logical_slots().iter()
                            .map(|slot| raw_cooking.uncook(slot.to_atom().as_view())
                                .to_canonical_string())
                            .collect::<Vec<_>>(),
                    }))?,
                    file.integration_atom = %raw_cooking.uncook(vacuum.expression().as_view())
                        .to_canonical_string(),
                    integration_atom_bytes = vacuum.expression().as_view().get_byte_size(),
                    atom_bytes = sector.atom.expression().as_view().get_byte_size(),
                    integration_atom_same_value_as_atom = vacuum.same_value(&sector.atom),
                    integration_atom_expression_equals_atom = vacuum.expression() == sector.atom.expression(),
                    integration_atom_literal_zero = vacuum.is_zero(),
                    formal_dimension_present = vacuum.expression().contains_symbol(GS.dim),
                    epsilon_present = vacuum.expression().contains_symbol(GS.dim_epsilon),
                    vacuum_mass_present = vacuum.expression().contains_symbol(GS.m_uv_vacuum),
                    expansion_mass_present = vacuum.expression().contains_symbol(GS.m_uv_expansion),
                    renormalization_scale_present = vacuum.expression().contains_symbol(GS.mu_r_sq),
                    consumed_loop_scale_present = vacuum.expression().contains_symbol(GS.integrated_loop_scale),
                    "Retained original signed W sector with complete ownership and symbolic scalings"
                );
            }
        }
        let reduced = graph
            .full_filter()
            .subtract(parent.subgraph())
            .subtract(&graph.initial_state_cut);
        let normalize = FinalIntegrandBuilder::final_normalizer(graph, &reduced);
        let normalize_scalar = FinalIntegrandBuilder::scalar_normalizer(graph, &reduced);
        let scalar_tail = |atom: &Atom| {
            normalize_scalar(atom)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.dim_epsilon)
                .with(0)
        };
        let boundary = graph
            .dummy_less_full_crown(parent.subgraph())
            .included_iter()
            .collect::<Vec<_>>();
        eyre::ensure!(boundary.len() == 2);
        let incoming = boundary
            .iter()
            .find(|hedge| graph.flow(**hedge) == Flow::Sink)
            .ok_or_else(|| eyre!("D09 incoming parent port missing"))?;
        let owner = graph[incoming];
        let mass = graph[owner].particle.mass_atom(&graph.model);
        let (projector, reference_p) = self.mass_projector(diagram, &Atom::var(GS.dim))?;
        let basis = diagram.loop_momentum_basis();
        let routed_p = basis.route_expression(&reference_p);
        eyre::ensure!(
            (basis.route_expression(&function!(GS.emr_mom, owner.0)) - &routed_p)
                .expand()
                .is_zero()
        );
        let mut signs = BTreeMap::new();
        for hedge in &boundary {
            let edge = graph[hedge];
            let routed = basis.route_expression(&function!(GS.emr_mom, edge.0));
            let sign = if (&routed - &routed_p).expand().is_zero() {
                1
            } else if (&routed + &routed_p).expand().is_zero() {
                -1
            } else {
                return Err(eyre!(
                    "D09 boundary is not the certified signed external momentum"
                ));
            };
            signs.insert(edge, sign);
        }
        let velocity = ctx.settings.on_shell_velocity()?;
        let at_point = |tensor: &SymbolicTensor<IndexedStructure>, lmb: &LoopMomentumBasis| {
            let routed = GS
                .erase_uv_momentum_provenance(tensor.expression())
                .replace_multiple(graph.uv_wrapped_replacement(parent.subgraph(), lmb, &[W_.x___]));
            let mapped = routed.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() > 0
                    && let Ok(edge) = usize::try_from(momentum.get(0))
                    && let Some(sign) = signs.get(&EdgeIndex(edge))
                {
                    let point = GS.on_shell_momentum(owner, &(Atom::num(*sign) * &mass), &velocity);
                    let indices = momentum
                        .iter()
                        .skip(1)
                        .map(|index| index.to_owned())
                        .collect::<Vec<_>>();
                    **output = GS.indexed_momentum(&point, &indices);
                }
            });
            tensor.with_checked_expression(
                GS.resolve_os_momenta(&mapped)
                    .replace(mass.to_pattern())
                    .with(1),
            )
        };
        let projector = projector
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
        let projector =
            projector.with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
        let projector = at_point(
            &normalize(&projector.try_into_indexed()?)?,
            &graph.loop_momentum_basis,
        )?;
        let phi = |tensor: &SymbolicTensor<IndexedStructure>,
                   lmb: &LoopMomentumBasis|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let projected = projector.clone().multiply(at_point(tensor, lmb)?)?;
            eyre::ensure!(
                projected.structure().canonical().is_scalar(),
                "D09 finalized mass probe leaves ports: projector={:?}; operand={:?}; product={:?}",
                projector.structure().logical_slots(),
                tensor.structure().logical_slots(),
                projected.structure().logical_slots(),
            );
            Ok(projected)
        };
        // Use the existing physical-shell owner to consume loop bookkeeping.
        // Final tensor bodies and scalar denominators retain their separate owners.
        let emit = |sector: &FourDSector| -> Result<SymbolicTensor<IndexedStructure>> {
            let emit_term = |term: FourDTerm| -> Result<SymbolicTensor<IndexedStructure>> {
                let numerator = normalize(&term.numerator)?;
                let denominators =
                    Atom::mul_many(term.denominators.into_iter().map(|denominator| {
                        GS.den(
                            usize::from(denominator.source_edge),
                            denominator.momentum,
                            denominator.mass_squared,
                            denominator.full_expr,
                        )
                    }));
                Ok(numerator
                    .with_checked_expression(numerator.expression() / scalar_tail(&denominators))?)
            };
            let mut terms = sector.physical_terms()?.into_iter();
            let first = terms
                .next()
                .ok_or_else(|| eyre!("D09 physical operand has no rational terms"))?;
            let first = emit_term(first)?;
            terms.try_fold(first, |sum, term| -> Result<_> {
                Ok(sum.try_add(&emit_term(term)?)?)
            })
        };
        // The same finalizer and point map must preserve an independently
        // known mass insertion. Use vertex color slots, dual to the removed
        // external propagators, exactly as in the native amputation controls.
        let spin = boundary
            .iter()
            .map(|hedge| {
                graph.underlying[*hedge]
                    .edge_spin_slots()
                    .next()
                    .map(|slot| slot.to_atom())
                    .ok_or_else(|| eyre!("D09 boundary spin slot missing"))
            })
            .collect::<Result<Vec<_>>>()?;
        let color = boundary
            .iter()
            .map(|hedge| {
                graph.underlying[*hedge]
                    .num_indices
                    .color_indices
                    .vertex_indices
                    .external_structure_iter()
                    .next()
                    .map(|slot| slot.to_atom())
                    .ok_or_else(|| eyre!("D09 boundary vertex color slot missing"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mass_control = AppliedFeynmanRule::new(
            &mass
                * function!(ETS.metric, &spin[0], &spin[1])
                * function!(ETS.metric, &color[0], &color[1]),
        )?
        .expr;
        eyre::ensure!(
            mass_control.structure().canonical().external_structure()
                == local.atom().structure().canonical().external_structure(),
            "the authored mass probe must use the actual parent's indexed interface"
        );
        let projected_mass_control = phi(&normalize(&mass_control)?, &graph.loop_momentum_basis)?;
        let calibration = crate::uv::overlap_numerators::NumeratorJetCache::default()
            .canonicalize(&on_shell_tests::scalar_component(
                projected_mass_control.expression(),
            )?)?;
        eyre::ensure!(
            !calibration.is_zero
                && calibration.key.is_one()
                && calibration.scalar_prefactor.is_one(),
            "D09 finalized mass projector must give exactly one on the native mass identity"
        );
        // A single explicit coordinate witness proves that these reconstructed
        // scalar kernels are not identically zero. No D-dimensional integral
        // or actual Projected4dSector transport is inferred from this check.
        let kernel_at_origin = |kernel: &Atom, loops: usize| -> Result<()> {
            eyre::ensure!(!kernel.is_zero());
            let origin = kernel
                .replace(function!(GS.localizing_integrand, W_.a_))
                .with(W_.a_)
                .replace(function!(GS.emr_mom, W_.a_, W_.b_))
                .with(0)
                .replace(GS.renormalization_localization_scale)
                .with(1);
            let inverse_origin =
                Atom::var(symbolica::atom::Symbol::PI).pow(i64::try_from(2 * loops)?);
            eyre::ensure!(
                (&origin * &inverse_origin)
                    .collect_compact_factors()
                    .is_one()
            );
            Ok(())
        };
        let co_graph = grow_factor(ctx, parent, given)?;
        let mut domains = std::collections::HashSet::new();
        let mut matched = BTreeSet::new();
        let mut frozen_domains = 0;
        for input in prefix.sectors() {
            eyre::ensure!(
                !input.atom.is_zero(),
                "D09 cannot match a pruned input by position"
            );
            let single = Full4dCts(FourDSectors::new(vec![input.clone()], Vec::new())?);
            let coordinate_lmb = compatible_lmb(ctx, parent, given, &single)?;
            let retained = input
                .active_components
                .iter()
                .flat_map(|(_, _, lmb)| lmb.loop_edges.iter())
                .chain(
                    input
                        .frozen_lmbs
                        .iter()
                        .flat_map(|lmb| lmb.loop_edges.iter()),
                )
                .copied()
                .collect::<BTreeSet<_>>();
            let mut quotient = coordinate_lmb.clone();
            if !given.subgraph().is_empty() {
                for index in (0..quotient.loop_edges.len()).rev() {
                    let index = LoopIndex(index);
                    if retained.contains(&quotient.loop_edges[index]) {
                        quotient.put_loop_to_ext(index);
                    }
                }
                graph.canonicalize_lmb_external_order(&mut quotient);
            }
            eyre::ensure!(
                quotient.loop_edges.len()
                    == graph.n_loops(parent.subgraph()) - graph.n_loops(given.subgraph())
            );
            let mut components = input.active_components.clone();
            components.push((
                parent.reduced_subgraph(given),
                parent.subgraph().clone(),
                quotient,
            ));
            let key = (components, input.frozen_lmbs.clone());
            eyre::ensure!(
                domains.insert(key.clone()),
                "D09 input domains must be distinct"
            );
            let outputs = local
                .active_sectors()
                .iter()
                .enumerate()
                .filter(|(_, output)| {
                    output.active_components == key.0 && output.frozen_lmbs == key.1
                })
                .collect::<Vec<_>>();
            eyre::ensure!(
                outputs.len() == 1,
                "D09 source lacks one exact full-domain output"
            );
            let (index, output) = outputs[0];
            eyre::ensure!(matched.insert(index), "D09 output domain was reused");
            let source = FourDSector {
                atom: finalize(input.atom.clone().multiply(co_graph.clone())?)?,
                ..input.clone()
            };
            let source = emit(&source)?;
            let local = emit(output)?;
            let frozen = scalar_tail(&Atom::mul_many(
                input
                    .frozen_lmbs
                    .iter()
                    .map(|lmb| GS.localizing_integrand(lmb)),
            ));
            let frozen_tensor = SymbolicTensor::infer_indexed(frozen.clone())?;
            eyre::ensure!(at_point(&frozen_tensor, &coordinate_lmb)?.expression() == &frozen);
            kernel_at_origin(
                &frozen,
                input
                    .frozen_lmbs
                    .iter()
                    .map(|lmb| lmb.loop_edges.len())
                    .sum(),
            )?;
            if !input.frozen_lmbs.is_empty() {
                frozen_domains += 1;
                eyre::ensure!(input.active_components.is_empty());
                eyre::ensure!(input.frozen_lmbs.as_slice() == std::slice::from_ref(given.lmb()));
                eyre::ensure!(
                    !normalize(&input.atom)?.is_zero() && !source.is_zero(),
                    "actual finite-child carrier must survive finalization"
                );
            }
            let source = phi(&source, &coordinate_lmb)?;
            let local = phi(&local, &coordinate_lmb)?;
            let difference = source.try_add(&local)?;
            on_shell_tests::assert_rational_zero(
                "D09 separately finalized source and signed parent local output",
                graph,
                parent,
                &(difference.expression() * &frozen),
                on_shell_tests::SpinReduction::FiniteComponents,
            )?;
            // An omitted/misnormalized mass subtraction cannot hide behind
            // the matched zero: the independently finalized wrong local
            // operand must leave the unit times the same nonzero kernel.
            let perturbed = FourDSector {
                atom: output.atom.try_add(&mass_control)?,
                ..output.clone()
            };
            let perturbed = phi(&emit(&perturbed)?, &coordinate_lmb)?;
            let wrong = source.try_add(&perturbed)?;
            let wrong = wrong.with_checked_expression(wrong.expression() * &frozen)?;
            let expected = projected_mass_control
                .with_checked_expression(projected_mass_control.expression() * &frozen)?;
            let controlled_residual = wrong.try_sub(&expected)?;
            on_shell_tests::assert_rational_zero(
                "D09 finalized mass perturbation leaves the certified nonzero frozen kernel",
                graph,
                parent,
                controlled_residual.expression(),
                on_shell_tests::SpinReduction::FiniteComponents,
            )?;
        }
        eyre::ensure!(matched.len() == local.active_sectors().len());
        eyre::ensure!(domains.len() == if given.subgraph().is_empty() { 1 } else { 2 });
        eyre::ensure!(frozen_domains == usize::from(!given.subgraph().is_empty()));
        let physical_finite = integrated.physical_finite_counterterm()?;
        let finite = normalize(&physical_finite)?;
        if given.subgraph().is_empty() {
            eyre::ensure!(
                !finite.is_zero(),
                "actual nonzero finite field correction must survive finalization"
            );
        }
        let parent_kernel = scalar_tail(&GS.localizing_integrand(parent.lmb()));
        kernel_at_origin(&parent_kernel, parent.lmb().loop_edges.len())?;
        let kernel_tensor = SymbolicTensor::infer_indexed(parent_kernel.clone())?;
        eyre::ensure!(
            at_point(&kernel_tensor, &graph.loop_momentum_basis)?.expression() == &parent_kernel
        );
        // Exercise the actual recursive finite-carrier owner relative to the
        // actual Integrated result, including consumed loop bookkeeping.
        // A zero result has no surviving recursive payload; it cannot be zipped.
        let completion = Full4dCts::recursion_input(
            local,
            integrated,
            ApproximationType::OS,
            false,
            parent.lmb(),
        )?;
        let completions = &completion.0.recursive_completion;
        if physical_finite.is_zero() {
            eyre::ensure!(completions.is_empty());
        } else {
            eyre::ensure!(completions.len() == 1);
            let completion = &completions[0];
            eyre::ensure!(completion.active_components.is_empty());
            eyre::ensure!(completion.frozen_lmbs.as_slice() == std::slice::from_ref(parent.lmb()));
            let transported = emit(completion)?;
            Self::certify_parent_tensor_zero(
                "D09 signed recursive finite carrier matches the actual physical Integrated output",
                &transported.try_sub(&finite)?,
            )?;
            let kernel = scalar_tail(&Atom::mul_many(
                completion
                    .frozen_lmbs
                    .iter()
                    .map(|lmb| GS.localizing_integrand(lmb)),
            ));
            eyre::ensure!(kernel == parent_kernel);
        }
        let finite_projection = phi(&finite, &graph.loop_momentum_basis)?;
        on_shell_tests::assert_rational_zero(
            "D09 actual finalized finite parent correction",
            graph,
            parent,
            &(finite_projection.expression() * &parent_kernel),
            on_shell_tests::SpinReduction::FiniteComponents,
        )?;
        let perturbed_finite = normalize(&physical_finite.try_add(&mass_control)?)?;
        let wrong = phi(&perturbed_finite, &graph.loop_momentum_basis)?;
        let wrong = wrong.with_checked_expression(wrong.expression() * &parent_kernel)?;
        let expected = projected_mass_control
            .with_checked_expression(projected_mass_control.expression() * &parent_kernel)?;
        on_shell_tests::assert_rational_zero(
            "D09 finalized finite mass perturbation leaves the certified nonzero parent kernel",
            graph,
            parent,
            wrong.try_sub(&expected)?.expression(),
            on_shell_tests::SpinReduction::FiniteComponents,
        )?;
        // Construct both actual emitted lanes before the complete-monomial
        // mass certificate consumes their original factorized payloads.
        {
            use crate::{
                graph::feynman_graph::FeynmanGraph,
                settings::global::OrientationPattern,
                uv::{
                    Integrands,
                    approx::{
                        CutStructure, OrientationProjection,
                        local_3d::Localizer,
                        projected_4d::{Local4dProjectionContext, Projected4dApproximation},
                    },
                    settings::FinalIntegrandDimension,
                },
            };
            use symbolica::id::Replacement;

            let retained_local = local.clone();
            let retained_finite = physical_finite.clone();
            let mut projection_graph = graph.clone();
            let mut settings = ctx.settings.clone();
            settings.local_uv_cts_from_expanded_4d_integrands = true;
            settings.final_integrand = FinalIntegrandDimension::ThreeD;
            eyre::ensure!(settings.generate_integrated && !settings.add_marker);
            let pattern = OrientationPattern::default();
            let numerator = projection_graph
                .numerator(&graph.full_filter(), &graph.empty_subgraph())?
                .state
                .expr;
            let options = projection_graph.denominator_only_cff_3d_expression_options();
            let surfaces = projection_graph.get_esurface_canonization(&graph.loop_momentum_basis);
            let contracted_edges =
                graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
            let production = projection_graph.generate_3d_expression_for_integrand(
                &contracted_edges,
                &surfaces,
                &options,
                Some(numerator.expression()),
            )?;
            let orientation =
                OrientationProjection::exact_expression(&production, &options, &pattern, true);
            eyre::ensure!(!orientation.orientation_ids()?.is_empty());
            let cuts = CutStructure::empty(&projection_graph);
            let localizer = Localizer::new(&cuts.cuts[0], orientation);
            let allowed_keys = cuts.cuts[0].residue_selector.generate_allowed_keys();
            eyre::ensure!(!allowed_keys.is_empty());

            let family = symbol!("gammalooprs::uv::numerator_family");
            let calls = |integrands: &Integrands| {
                let mut calls = BTreeSet::new();
                for (_, root) in integrands.iter() {
                    root.visitor(&mut |view| {
                        if let AtomView::Fun(call) = view
                            && call.get_symbol() == family
                        {
                            calls.insert(view.to_owned());
                        }
                        true
                    });
                }
                calls
            };
            let check_store = |integrands: &Integrands| -> Result<()> {
                // Reuse the store's flatness, binding and collision admission.
                let checked = integrands
                    .clone()
                    .with_numerators(integrands.numerators().iter().cloned())?;
                eyre::ensure!(&checked == integrands);
                for call in calls(integrands) {
                    let AtomView::Fun(call) = call.as_view() else {
                        unreachable!()
                    };
                    let matches = integrands
                        .numerators()
                        .iter()
                        .filter(|entry| {
                            call.get_nargs() == entry.tags.len() + entry.args.len()
                                && entry
                                    .tags
                                    .iter()
                                    .enumerate()
                                    .all(|(index, tag)| call.get(index) == tag.as_view())
                        })
                        .count();
                    eyre::ensure!(
                        matches == 1,
                        "D09 CFF call lacks one complete retained binding"
                    );
                }
                for (_, root) in integrands.iter() {
                    eyre::ensure!(
                        ![AGS.gamma, AGS.gamma0, GS.u, GS.ubar]
                            .into_iter()
                            .any(|head| root.contains_symbol(head)),
                        "D09 CFF dispatch root contains an explicit spin tensor"
                    );
                }
                Ok(())
            };

            // The producer emits one sector per already grouped full domain,
            // in this stable order. Never zip the ungrouped Taylor sectors.
            let groups = local.projection_sectors()?;
            eyre::ensure!(!groups.is_empty());
            let mut full_domains = std::collections::HashSet::new();
            for group in &groups {
                eyre::ensure!(
                    full_domains
                        .insert((group.active_components.clone(), group.frozen_lmbs.clone())),
                    "D09 projection retained duplicate full domains"
                );
            }
            let mut projection_context = Local4dProjectionContext::default();
            eyre::ensure!(projection_context.source_preparation_builds == 0);
            eyre::ensure!(projection_context.denominator_reconstruction_certificates == 0);
            projection_context.certify_denominator_reconstruction = true;
            let projected =
                Projected4dApproximation::new(localizer, &mut projection_graph, &settings)
                    .project_local_4d(local, &mut projection_context)?;
            eyre::ensure!(projected.sectors().len() == groups.len());
            let certificates = projection_context.denominator_reconstruction_certificates;
            eyre::ensure!(certificates > 0);
            eyre::ensure!(certificates == projection_context.source_preparation_builds);
            eyre::ensure!(projection_context.projected_component_requests >= certificates);
            let mut projected_families = 0;
            for (index, group) in groups.iter().enumerate() {
                let sector = &projected.sectors()[index];
                let expected_kernel = group.frozen_lmbs.iter().fold(Atom::one(), |product, lmb| {
                    product * GS.localizing_integrand(lmb)
                });
                eyre::ensure!(sector.frozen_factor == expected_kernel);
                // The preceding rational gate supplies the nonzero origin
                // witness for these exact same frozen LMB kernels.
                eyre::ensure!(!sector.frozen_factor.is_zero());
                kernel_at_origin(
                    &scalar_tail(&sector.frozen_factor),
                    group
                        .frozen_lmbs
                        .iter()
                        .map(|lmb| lmb.loop_edges.len())
                        .sum(),
                )?;
                check_store(&sector.coefficient)?;
                projected_families += sector.coefficient.numerators().len();
            }
            eyre::ensure!(projected_families > 0);

            // Certify the actual finite localizer's ordered carrier transport
            // relative to this very Integrated output, before finalization/Phi.
            eyre::ensure!(
                !physical_finite
                    .expression()
                    .contains_symbol(GS.integrated_loop_scale)
                    && !physical_finite
                        .expression()
                        .contains_symbol(GS.localizing_integrand)
            );
            let localized = localizer.localize(&physical_finite, &mut projection_graph, parent)?;
            let mut finite_calls = 0;
            if !physical_finite.is_zero() {
                let kernel = GS.localizing_integrand(parent.lmb());
                eyre::ensure!(!kernel.is_zero());
                eyre::ensure!(localized.frozen_integrands.iter().count() == allowed_keys.len());
                for (key, actual_kernel) in localized.frozen_integrands.iter() {
                    eyre::ensure!(allowed_keys.contains(key) && actual_kernel == &kernel);
                    kernel_at_origin(&scalar_tail(actual_kernel), parent.lmb().loop_edges.len())?;
                }
                for (host, source_map, integrands) in localized.active.iter_orientations() {
                    check_store(integrands)?;
                    eyre::ensure!(integrands.numerators().len() == 1);
                    let entry = &integrands.numerators()[0];
                    let expected = orientation.map_numerator(
                        &projection_graph,
                        host,
                        source_map,
                        physical_finite.expression(),
                    )?;
                    let expected = physical_finite.with_checked_expression(expected)?;
                    for call in calls(integrands) {
                        let AtomView::Fun(call) = call.as_view() else {
                            unreachable!()
                        };
                        let replacements = entry.args.iter().enumerate().map(|(index, formal)| {
                            Replacement::new(
                                Atom::from(formal.clone()).to_pattern(),
                                call.get(entry.tags.len() + index).to_owned(),
                            )
                        });
                        // Specialize one retained typed body as a proof view;
                        // never substitute the body into its dispatch roots.
                        let transported = entry.rhs.with_checked_expression(
                            entry.rhs.expression().replace_multiple(replacements),
                        )?;
                        Self::certify_parent_tensor_zero(
                            "D09 actual CFF finite carrier and ordered branch arguments",
                            &transported.try_sub(&expected)?,
                        )?;
                        finite_calls += 1;
                    }
                }
                eyre::ensure!(finite_calls > 0);
            }
            // This is the actual final builder with both retained operands.
            // into_integrands keeps every numerator body separate from roots.
            let finalized = FinalIntegrandBuilder::new(localizer, &settings)
                .build_projected(&mut projection_graph, parent, &projected, integrated)?
                .into_integrands();
            check_store(&finalized)?;
            eyre::ensure!(finalized.scalar_definitions().is_empty());
            eyre::ensure!(!finalized.numerators().is_empty());
            eyre::ensure!(finalized.iter().any(|(_, root)| !root.is_zero()));
            eyre::ensure!(finalized.iter().count() == allowed_keys.len());
            eyre::ensure!(finalized.iter().all(|(key, _)| allowed_keys.contains(key)));
            eyre::ensure!(local == &retained_local);
            eyre::ensure!(integrated.physical_finite_counterterm()? == retained_finite);
            debug_tags!(#uv, #local, #four_d, #cff, #integrated, #summary;
                stage = "native_d09_parent_cff_reconstruction",
                source = %graph.name, given = %given.subgraph().string_label(),
                full_domains = groups.len(), certificates,
                file.grouped_domains = ?groups.iter().map(|group| (&group.active_components, &group.frozen_lmbs)).collect::<Vec<_>>(),
                component_requests = projection_context.projected_component_requests,
                projected_families, finite_calls,
                final_families = finalized.numerators().len(),
                final_root_bytes = finalized.iter().map(|(_, atom)| atom.as_view().get_byte_size()).sum::<usize>(),
                final_body_bytes = finalized.numerators().iter().map(|entry| entry.rhs.expression().as_view().get_byte_size()).sum::<usize>(),
                actual_parent_cff_constructed = true, finite_carrier_transport_checked = true,
                emitted_mass_cancellation_proved = false, dimensional_R2_condition_proved = false,
                "Actual D09 parent CFF retains its certified factorized operands"
            );
            // Construct the actual source lanes next. Keep the external matrix
            // open in this separate projection view; no mass projector or
            // matched-point specialization enters source/assignment selection.
            {
                use crate::numerator::symbolica_ext::NumeratorAtomExt;
                use crate::uv::approx::direct_3d::{Direct3dCts, DirectResidueBranches};
                use std::{collections::HashMap, sync::Arc};

                let retained_prefix = prefix.clone();
                let retained_global_projector = projection_graph.global_prefactor.projector.clone();
                let mut source_graph = projection_graph.clone();
                source_graph.global_prefactor.projector =
                    SymbolicTensor::infer_indexed(Atom::one())?;
                eyre::ensure!(
                    source_graph.global_prefactor.num == projection_graph.global_prefactor.num
                );
                eyre::ensure!(source_graph.overall_factor == projection_graph.overall_factor);
                eyre::ensure!(
                    source_graph.loop_momentum_basis == projection_graph.loop_momentum_basis
                );
                // The production expression above uses this full graph's
                // numerator without its external global projector. Reusing it
                // retains the actual stored energy bounds/maps for this view.
                eyre::ensure!(
                    !retained_global_projector
                        .expression()
                        .contains_symbol(GS.loop_mom)
                );
                let external_owners = boundary
                    .iter()
                    .map(|hedge| graph[hedge].0)
                    .collect::<BTreeSet<_>>();
                let mut unexpected_external = Vec::new();
                retained_global_projector.expression().visitor(&mut |view| {
                    if let AtomView::Fun(factor) = view
                        && [GS.u, GS.ubar, GS.emr_mom].contains(&factor.get_symbol())
                        && (factor.get_nargs() == 0
                            || usize::try_from(factor.get(0))
                                .ok()
                                .is_none_or(|edge| !external_owners.contains(&edge)))
                    {
                        unexpected_external.push(view.to_owned());
                    }
                    true
                });
                eyre::ensure!(
                    unexpected_external.is_empty(),
                    "the external projection view cannot remove a hard-loop-dependent factor: {unexpected_external:?}"
                );
                // Reuse the already certified parent coefficients in the SAME
                // open external graph view as the source. This is an actual
                // final-builder call with the original Integrated result.
                let parent_assembly = Mutex::new(ProjectedAssemblyCertificate::default());
                let open_parent = FinalIntegrandBuilder::new(localizer, &settings)
                    .with_projected_assembly_certificate(&parent_assembly)
                    .build_projected(&mut source_graph, parent, &projected, integrated)?
                    .into_integrands();
                let parent_assembly = parent_assembly
                    .into_inner()
                    .map_err(|_| eyre!("D09 parent assembly certificate was poisoned"))?;
                parent_assembly.certify(
                    &source_graph,
                    parent,
                    &projected,
                    localizer,
                    &settings,
                    &open_parent,
                )?;
                let mut wrong_key = parent_assembly.clone();
                wrong_key
                    .emissions
                    .first_mut()
                    .ok_or_else(|| eyre!("D09 parent emitted no residue key"))?
                    .key
                    .selector_host
                    .0 ^= 1;
                eyre::ensure!(
                    wrong_key
                        .certify(
                            &source_graph,
                            parent,
                            &projected,
                            localizer,
                            &settings,
                            &open_parent
                        )
                        .is_err(),
                    "D09 assembly certificate accepted a changed full residue key"
                );
                check_store(&open_parent)?;
                eyre::ensure!(open_parent.scalar_definitions().is_empty());
                eyre::ensure!(!open_parent.numerators().is_empty());
                eyre::ensure!(open_parent.iter().count() == allowed_keys.len());
                eyre::ensure!(
                    open_parent
                        .iter()
                        .all(|(key, _)| allowed_keys.contains(key))
                );
                let source_result = if given.subgraph().is_empty() {
                    let direct = Direct3dCts::root(&source_graph, localizer)?;
                    FinalIntegrandBuilder::new(localizer, &settings)
                        .build_direct(&mut source_graph, given, &direct, &IntegratedCts::root())?
                        .into_integrands()
                } else {
                    let sectors = prefix.sectors().collect::<Vec<_>>();
                    let active = sectors
                        .iter()
                        .filter(|sector| sector.frozen_lmbs().is_empty())
                        .map(|sector| (*sector).clone())
                        .collect::<Vec<_>>();
                    let finite = sectors
                        .iter()
                        .filter(|sector| !sector.frozen_lmbs().is_empty())
                        .copied()
                        .collect::<Vec<_>>();
                    eyre::ensure!(active.len() == 1 && finite.len() == 1);
                    eyre::ensure!(finite[0].active_components.is_empty());
                    eyre::ensure!(finite[0].frozen_lmbs() == std::slice::from_ref(given.lmb()));
                    eyre::ensure!(graph.n_loops(given.subgraph()) == 1);
                    let active_source = Local4dCts(FourDSectors::new(active, Vec::new())?);
                    Self::certify_parent_tensor_zero(
                        "D09 actual source prefix active/finite recomposition",
                        &prefix
                            .atom()
                            .try_sub(&active_source.atom().try_add(&finite[0].atom)?)?,
                    )?;
                    let mut physical_child_finite =
                        finite[0].atom.with_checked_expression(Atom::Zero)?;
                    let terms = finite[0].physical_terms()?;
                    eyre::ensure!(!terms.is_empty());
                    for term in terms {
                        eyre::ensure!(
                            term.denominators.is_empty(),
                            "a retained finite child cannot own physical propagators"
                        );
                        physical_child_finite = physical_child_finite.try_add(&term.numerator)?;
                    }
                    eyre::ensure!(!physical_child_finite.is_zero());
                    eyre::ensure!(
                        !physical_child_finite
                            .expression()
                            .contains_symbol(GS.localizing_integrand)
                    );
                    let bookkeeping = physical_child_finite.with_checked_expression(
                        physical_child_finite.expression()
                            * Atom::var(GS.integrated_loop_scale).pow(4),
                    )?;
                    Self::certify_parent_tensor_zero(
                        "D09 actual finite child retains exactly its consumed one-loop measure",
                        &finite[0].atom.try_sub(&bookkeeping)?,
                    )?;
                    let active_groups = active_source.projection_sectors()?;
                    let mut child_context = Local4dProjectionContext {
                        certify_denominator_reconstruction: true,
                        ..Default::default()
                    };
                    eyre::ensure!(child_context.source_preparation_builds == 0);
                    eyre::ensure!(child_context.denominator_reconstruction_certificates == 0);
                    let child_projected =
                        Projected4dApproximation::new(localizer, &mut source_graph, &settings)
                            .project_local_4d(&active_source, &mut child_context)?;
                    eyre::ensure!(child_projected.sectors().len() == active_groups.len());
                    eyre::ensure!(child_context.source_preparation_builds > 0);
                    eyre::ensure!(
                        child_context.denominator_reconstruction_certificates
                            == child_context.source_preparation_builds
                    );
                    for (ordinal, group) in active_groups.iter().enumerate() {
                        let expected = group
                            .frozen_lmbs()
                            .iter()
                            .fold(Atom::one(), |product, lmb| {
                                product * GS.localizing_integrand(lmb)
                            });
                        eyre::ensure!(child_projected.sectors()[ordinal].frozen_factor == expected);
                    }
                    // This zero means no ADDITIONAL finite child term here.
                    // The actual retained finite completion is localized below.
                    let source_assembly = Mutex::new(ProjectedAssemblyCertificate::default());
                    let active_result = FinalIntegrandBuilder::new(localizer, &settings)
                        .with_projected_assembly_certificate(&source_assembly)
                        .build_projected(
                            &mut source_graph,
                            given,
                            &child_projected,
                            &IntegratedCts::root(),
                        )?
                        .into_integrands();
                    let source_assembly = source_assembly.into_inner().map_err(|_| {
                        eyre!("D09 active source assembly certificate was poisoned")
                    })?;
                    source_assembly.certify(
                        &source_graph,
                        given,
                        &child_projected,
                        localizer,
                        &settings,
                        &active_result,
                    )?;
                    eyre::ensure!(
                        source_assembly.integrated_is_zero == Some(true)
                            && source_assembly.integrated_key_domain.len() == 1
                            && source_assembly.integrated_key_domain[0]
                                .source_edge_energy_map()
                                .is_none()
                            && Some(&source_assembly.integrated_key_domain[0].selector_host)
                                == orientation.orientation_ids()?.first()
                            && source_assembly
                                .skips
                                .iter()
                                .all(|skip| matches!(skip, ProjectedAssemblySkip::EmptyIntegrated)),
                        "active source must retain the actual integrated-zero policy and its first production fallback key"
                    );
                    // Ordinals belong to the exact active_source projection groups
                    // checked above, before the retained finite completion is added.
                    check_store(&active_result)?;
                    let finite_localized =
                        localizer.localize(&physical_child_finite, &mut source_graph, given)?;
                    let actual_kernel = GS.localizing_integrand(given.lmb());
                    eyre::ensure!(!actual_kernel.is_zero());
                    eyre::ensure!(
                        finite_localized.frozen_integrands.iter().count() == allowed_keys.len()
                    );
                    for (key, kernel) in finite_localized.frozen_integrands.iter() {
                        eyre::ensure!(allowed_keys.contains(key) && kernel == &actual_kernel);
                    }
                    let finite_branches =
                        DirectResidueBranches::from_transient(&finite_localized.combine()?)?;
                    let mut finite_entries = HashMap::new();
                    for (key, integrands) in finite_branches.iter_keys() {
                        check_store(integrands)?;
                        eyre::ensure!(integrands.numerators().len() == 1);
                        eyre::ensure!(
                            finite_entries
                                .insert(key.clone(), Arc::clone(&integrands.numerators()[0]))
                                .is_none()
                        );
                    }
                    let reduced_source = graph
                        .full_filter()
                        .subtract(given.subgraph())
                        .subtract(&graph.initial_state_cut);
                    let outside = source_graph
                        .numerator(&reduced_source, given.subgraph())?
                        .state
                        .expr
                        .simplify_algebra(&idenso::tensor::AlgebraSettings {
                            color: Some(ColorSimplifySettings::default()),
                            ..Default::default()
                        })?
                        .multiply(source_graph.global_numerator()?)?;
                    eyre::ensure!(
                        !outside.is_scalar(),
                        "the actual D09 finite child has a tensor cograph"
                    );
                    let complete_source =
                        physical_child_finite.clone().multiply(outside.clone())?;
                    // Compare the actual raw-grown input before any CFF map or
                    // Phi. The F/G witness below only certifies transport of
                    // complete_source and cannot supply this separate join.
                    let raw_reduced = parent.reduced_subgraph(given);
                    eyre::ensure!(
                        graph.paired_edges(&raw_reduced)
                            == source_graph.paired_edges(&reduced_source),
                        "D09 raw and source cographs retain the same physical edge occurrences"
                    );
                    let global = source_graph.global_numerator()?;
                    let raw_outside = graph
                        .numerator(&raw_reduced, given.subgraph())?
                        .to_d_dim(GS.dim)?
                        .color_simplify()?
                        .state
                        .expr
                        .multiply(global.clone())?;
                    Self::certify_parent_tensor_zero(
                        "D09 raw and source cographs have the same complete physical-four numerator and one global",
                        &normalize(&raw_outside)?.try_sub(&normalize(&outside)?)?,
                    )?;
                    let shared_denominator = graph.denominator(&raw_reduced, ctx.model, |_| 1);
                    let source_denominator =
                        source_graph.denominator(&reduced_source, &source_graph.model, |_| 1);
                    eyre::ensure!(
                        shared_denominator == source_denominator,
                        "D09 cograph momenta, masses, powers and full prescribed polynomials differ"
                    );
                    let source_tree = source_graph
                        .tree_edges
                        .subtract(&source_graph.initial_state_cut);
                    eyre::ensure!(reduced_source.includes(&source_tree));
                    eyre::ensure!(
                        source_denominator
                            == source_graph.denominator(
                                &reduced_source.subtract(&source_tree),
                                &source_graph.model,
                                |_| 1,
                            ) * source_graph.denominator(
                                &source_tree,
                                &source_graph.model,
                                |_| 1,
                            ),
                        "D09 source loop and production-tree denominators must partition the raw cograph once"
                    );
                    let raw_source = FourDSector {
                        atom: finalize(finite[0].atom.clone().multiply(co_graph.clone())?)?,
                        ..finite[0].clone()
                    };
                    let raw_terms = raw_source.physical_terms()?;
                    eyre::ensure!(
                        raw_terms.len() == 1,
                        "D09 raw finite input must retain its one factorized cograph denominator"
                    );
                    let literal_denominator =
                        Atom::mul_many(raw_terms[0].denominators.iter().map(|denominator| {
                            GS.den(
                                usize::from(denominator.source_edge),
                                denominator.momentum.clone(),
                                denominator.mass_squared.clone(),
                                denominator.full_expr.clone(),
                            )
                        }));
                    eyre::ensure!(literal_denominator == shared_denominator);
                    let raw_emitted = emit(&raw_source)?.multiply(global)?;
                    // Cancel only this independently matched literal product;
                    // no denominator is cleared across a Phi or forest sum.
                    let raw_numerator = raw_emitted.with_checked_expression(
                        raw_emitted.expression() * scalar_tail(&literal_denominator),
                    )?;
                    eyre::ensure!(!raw_numerator.expression().contains_symbol(GS.den));
                    Self::certify_parent_tensor_zero(
                        "D09 raw frozen-child input equals the actual unmapped source before Phi",
                        &normalize(&raw_numerator)?.try_sub(&normalize(&complete_source)?)?,
                    )?;
                    eyre::ensure!(
                        raw_source.frozen_lmbs() == finite[0].frozen_lmbs()
                            && actual_kernel == GS.localizing_integrand(given.lmb())
                    );
                    eyre::ensure!(
                        complete_source.structure().canonical().external_structure()
                            == numerator.structure().canonical().external_structure()
                    );
                    let completed_branches = finite_branches.multiply_key_mapped(
                        orientation,
                        &source_graph,
                        &outside,
                        DirectResidueBranches::numerator_scope(),
                    )?;
                    let mut witnessed = 0usize;
                    let mut source_keys = 0usize;
                    for (key, integrands) in completed_branches.iter_keys() {
                        source_keys += 1;
                        check_store(integrands)?;
                        let finite_entry = finite_entries
                            .get(key)
                            .ok_or_else(|| eyre!("finite child lost its exact residue key"))?;
                        eyre::ensure!(
                            integrands
                                .numerators()
                                .iter()
                                .any(|entry| Arc::ptr_eq(entry, finite_entry))
                        );
                        let outside_entries = integrands
                            .numerators()
                            .iter()
                            .filter(|entry| !Arc::ptr_eq(entry, finite_entry))
                            .collect::<Vec<_>>();
                        eyre::ensure!(
                            outside_entries.len() == 1 && integrands.numerators().len() == 2
                        );
                        let outside_entry = outside_entries[0];
                        let actual_calls = calls(integrands).into_iter().collect::<Vec<_>>();
                        for (cut, root) in integrands.iter() {
                            let Some(rows) = root.coefficient_list_exact(&actual_calls) else {
                                debug_tags!(#uv, #local, #cff, #summary;
                                    stage = "native_d09_source_monomial_refused", file.residue_key = ?key,
                                    file.cut = ?cut, file.root = %root.to_canonical_string(),
                                    calls = actual_calls.len(), root_bytes = root.as_view().get_byte_size(),
                                    "Existing bounded complete-family-call polynomial owner refused; no fallback"
                                );
                                return Err(eyre!(
                                    "D09 source complete-call polynomial extraction refused"
                                ));
                            };
                            let reconstructed = rows
                                .iter()
                                .map(|(monomial, coefficient)| monomial * coefficient)
                                .sum::<Atom>();
                            let residual = (root - &reconstructed)
                                .coefficient_list_exact(&actual_calls)
                                .ok_or_else(|| {
                                    eyre!("D09 source complete-call reconstruction refused")
                                })?;
                            eyre::ensure!(
                                residual
                                    .iter()
                                    .all(|(_, coefficient)| coefficient.is_zero())
                            );
                            for (monomial, coefficient) in rows {
                                eyre::ensure!(
                                    !coefficient.contains_symbol(family),
                                    "D09 complete-call extraction left an unbound family in its coefficient"
                                );
                                if coefficient.is_zero() {
                                    continue;
                                }
                                let AtomView::Mul(product) = monomial.as_view() else {
                                    return Err(eyre!(
                                        "D09 actual finite source did not retain its complete two-factor call monomial"
                                    ));
                                };
                                let factors = product
                                    .iter()
                                    .map(|factor| factor.to_owned())
                                    .collect::<Vec<_>>();
                                eyre::ensure!(factors.len() == 2);
                                let matches = |call: &Atom, entry: &crate::integrands::process::param_builder::NumeratorDefinition| {
                                    matches!(call.as_view(), AtomView::Fun(call) if call.get_symbol() == family
                                        && call.get_nargs() == entry.tags.len() + entry.args.len()
                                        && entry.tags.iter().enumerate().all(|(index, tag)| call.get(index) == tag.as_view()))
                                };
                                let child_call = factors
                                    .iter()
                                    .filter(|call| matches(call, finite_entry))
                                    .collect::<Vec<_>>();
                                let outside_call = factors
                                    .iter()
                                    .filter(|call| matches(call, outside_entry))
                                    .collect::<Vec<_>>();
                                eyre::ensure!(
                                    child_call.len() == 1
                                        && outside_call.len() == 1
                                        && child_call[0] != outside_call[0]
                                );
                                // Every actual nonzero-coefficient monomial,
                                // under every actual key/cut, gets the typed
                                // transport proof before materialization.
                                {
                                    let bind = |call: &Atom, entry: &crate::integrands::process::param_builder::NumeratorDefinition| -> Result<_> {
                                        let AtomView::Fun(function) = call.as_view() else { unreachable!() };
                                        let replacements = entry.args.iter().enumerate().map(|(index, formal)| Replacement::new(
                                            Atom::from(formal.clone()).to_pattern(), function.get(entry.tags.len() + index).to_owned(),
                                        )).collect::<Vec<_>>();
                                        eyre::ensure!(entry.lhs.replace_multiple(&replacements) == *call);
                                        Ok(entry.rhs.with_checked_expression(entry.rhs.expression().replace_multiple(&replacements))?)
                                    };
                                    let child_body = bind(child_call[0], finite_entry)?;
                                    let outside_body = bind(outside_call[0], outside_entry)?;
                                    for (label, body, original) in [
                                        (
                                            "actual finite child call",
                                            &child_body,
                                            &physical_child_finite,
                                        ),
                                        (
                                            "actual finite child cograph call",
                                            &outside_body,
                                            &outside,
                                        ),
                                    ] {
                                        let expected = original.with_checked_expression(
                                            key.map_numerator(
                                                orientation,
                                                &source_graph,
                                                original.expression(),
                                            )?,
                                        )?;
                                        Self::certify_parent_tensor_zero(
                                            label,
                                            &body.try_sub(&expected)?,
                                        )?;
                                    }
                                    let complete_body = child_body.multiply(outside_body)?;
                                    eyre::ensure!(
                                        complete_body.structure().canonical().external_structure()
                                            == complete_source
                                                .structure()
                                                .canonical()
                                                .external_structure()
                                    );
                                    let expected = complete_source.with_checked_expression(
                                        key.map_numerator(
                                            orientation,
                                            &source_graph,
                                            complete_source.expression(),
                                        )?,
                                    )?;
                                    Self::certify_parent_tensor_zero(
                                        "actual same-key finite child/cograph complete monomial",
                                        &complete_body.try_sub(&expected)?,
                                    )?;
                                    debug_tags!(#uv, #local, #four_d, #cff, #summary;
                                        stage = "native_d09_source_complete_monomial", file.residue_key = ?key,
                                        file.cut = ?cut, file.monomial = %monomial.to_canonical_string(),
                                        file.child_call = %child_call[0].to_canonical_string(),
                                        file.cograph_call = %outside_call[0].to_canonical_string(),
                                        file.frozen_domain = ?finite[0].frozen_lmbs(),
                                        file.active_domains = ?active_groups.iter().map(|group| (&group.active_components, group.frozen_lmbs())).collect::<Vec<_>>(),
                                        file.kernel = %actual_kernel.to_canonical_string(),
                                        whole_parent_rank = complete_body.rank(), actual_call_pair_bound = true,
                                        mass_projector_applied = false,
                                        "Two actual partial families form one native open parent matrix under the same exact source key"
                                    );
                                    witnessed += 1;
                                }
                            }
                        }
                    }
                    eyre::ensure!(source_keys == finite_entries.len() && witnessed > 0);
                    let normalize_source =
                        FinalIntegrandBuilder::final_normalizer(&source_graph, &reduced_source);
                    let normalize_source_scalar =
                        FinalIntegrandBuilder::scalar_normalizer(&source_graph, &reduced_source);
                    let finite_result = completed_branches
                        .resolved_scalars()?
                        .fallible_map(|_, atom| {
                            Ok(normalize_source_scalar(atom)
                                .replace(GS.m_uv_expansion)
                                .with(GS.m_uv_vacuum)
                                .replace(GS.dim_epsilon)
                                .with(0))
                        })?
                        .map_numerators(|entry| normalize_source(&entry.rhs))?
                        .materialize(false)?;
                    check_store(&finite_result)?;
                    active_result.zip_add([finite_result])?
                };
                check_store(&source_result)?;
                eyre::ensure!(source_result.scalar_definitions().is_empty());
                eyre::ensure!(!source_result.numerators().is_empty());
                eyre::ensure!(source_result.iter().count() == allowed_keys.len());
                eyre::ensure!(
                    source_result
                        .iter()
                        .all(|(key, _)| allowed_keys.contains(key))
                );
                // Admit a common spatial chart without changing any component's
                // independent hard or frozen coordinates. A nontrivial chart
                // needs the actual component owner, not a low-rank map guess.
                let mut component_coordinates = BTreeSet::new();
                let mut frozen_kernels = BTreeSet::from([
                    GS.localizing_integrand(parent.lmb()),
                    GS.localizing_integrand(given.lmb()),
                ]);
                for sector in prefix.sectors().chain(local.active_sectors().iter()) {
                    frozen_kernels.extend(
                        sector
                            .frozen_lmbs
                            .iter()
                            .map(|lmb| GS.localizing_integrand(lmb)),
                    );
                    for lmb in sector
                        .active_components
                        .iter()
                        .map(|(_, _, lmb)| lmb)
                        .chain(&sector.frozen_lmbs)
                    {
                        component_coordinates.extend(lmb.loop_edges.iter().copied());
                    }
                }
                component_coordinates.extend(parent.lmb().loop_edges.iter().copied());
                component_coordinates.extend(given.lmb().loop_edges.iter().copied());
                for edge in &component_coordinates {
                    eyre::ensure!(
                        graph.loop_momentum_basis.loop_edges.contains(edge),
                        "D09 emitted Phi requires the actual component chart for non-global loop coordinate {edge:?}"
                    );
                    let routed =
                        graph
                            .loop_momentum_basis
                            .loop_atom::<Atom>(*edge, GS.emr_mom, &[], true)
                            + graph.loop_momentum_basis.ext_atom::<Atom>(
                                *edge,
                                GS.emr_mom,
                                &[],
                                true,
                            );
                    eyre::ensure!(
                        (routed - function!(GS.emr_mom, edge.0)).expand().is_zero(),
                        "D09 emitted Phi cannot change an independent component coordinate"
                    );
                }
                let spatial_routing = graph.uv_spatial_wrapped_replacement(
                    &graph.full_filter(),
                    &graph.loop_momentum_basis,
                    &[W_.x___],
                );
                let tree_routing = graph.uv_wrapped_replacement(
                    &graph.tree_edges,
                    &graph.loop_momentum_basis,
                    &[W_.x___],
                );
                let physical_energies = graph
                    .iter_edges()
                    .filter(|(pair, _, data)| pair.is_paired() && !data.data.is_dummy)
                    .map(|(_, edge, _)| {
                        Replacement::new(
                            GS.ose(edge).to_pattern(),
                            GS.ose_full(
                                edge,
                                EdgeIndex(0),
                                graph[edge].particle.mass_atom(&graph.model),
                                None,
                            ),
                        )
                    })
                    .collect::<Vec<_>>();
                let coordinates = |atom: &Atom| {
                    GS.resolve_os_momenta(atom)
                        .replace_multiple(&physical_energies)
                        .replace_multiple(&spatial_routing)
                        .replace_multiple(&tree_routing)
                };
                let emitted_point = |atom: &Atom| -> Result<Atom> {
                    eyre::ensure!(
                        GS.erase_uv_momentum_provenance(atom) == *atom,
                        "D09 emitted Phi requires every local momentum owner to be consumed before final observation"
                    );
                    let routed = coordinates(atom);
                    eyre::ensure!(
                        coordinates(&routed) == routed,
                        "D09 emitted Phi chart binding is not idempotent"
                    );
                    eyre::ensure!(
                        !routed.contains_symbol(GS.ose) && !routed.contains_symbol(GS.loop_mom),
                        "D09 emitted Phi has an unbound physical energy or loop-coordinate owner"
                    );
                    let mut unsupported = Vec::new();
                    let pointed = routed.replace_map(|view, _, output| {
                        let AtomView::Fun(momentum) = view else {
                            return;
                        };
                        if ![GS.emr_mom, GS.emr_vec].contains(&momentum.get_symbol()) {
                            return;
                        }
                        if momentum.get_nargs() != 2 {
                            unsupported.push(view.to_owned());
                            return;
                        }
                        let Ok(edge) = usize::try_from(momentum.get(0)).map(EdgeIndex) else {
                            unsupported.push(view.to_owned());
                            return;
                        };
                        let index = momentum.get(1).to_owned();
                        if let Some(sign) = signs.get(&edge) {
                            let point =
                                GS.on_shell_momentum(owner, &(Atom::num(*sign) * &mass), &velocity);
                            let vector = GS.resolve_os_momenta(
                                &GS.indexed_momentum(&point, std::slice::from_ref(&index)),
                            );
                            **output = if momentum.get_symbol() == GS.emr_vec {
                                vector
                                    - GS.resolve_os_momenta(
                                        &GS.indexed_momentum(&point, &[GS.cind(0)]),
                                    ) * GS.energy_delta(&index)
                            } else {
                                vector
                            };
                        } else if momentum.get_symbol() == GS.emr_mom
                            && !(1..=3).any(|component| index == GS.cind(component))
                        {
                            // A full vector could introduce an unassigned time
                            // component during the later finite contraction.
                            unsupported.push(view.to_owned());
                        }
                    });
                    eyre::ensure!(
                        unsupported.is_empty(),
                        "D09 emitted Phi has unresolved momentum dependencies: {unsupported:?}"
                    );
                    let pointed = pointed.replace(mass.to_pattern()).with(1);
                    eyre::ensure!(
                        coordinates(&pointed) == pointed,
                        "D09 emitted Phi changed the established final chart"
                    );
                    Ok(pointed)
                };
                for kernel in &frozen_kernels {
                    eyre::ensure!(
                        emitted_point(kernel)? == *kernel,
                        "D09 emitted Phi changed a frozen localization kernel"
                    );
                }
                let cooking = CookSettings::indices()
                    .with_mode(CookMode::ReversibleEncoding)
                    .with_representation_payloads(true, true);
                let explicit_momenta = |tensor: &SymbolicTensor<IndexedStructure>| -> Result<_> {
                    use idenso::tensor::{SymbolicNetExt, SymbolicNetParse};
                    use spenso::network::parsing::{
                        ParseSettings, SchoonschipExpansionMode, ShorthandParsing,
                    };
                    // Proper Schoonschip vector slots are materialized through
                    // their existing owner on one complete factorized matrix.
                    let explicit = cooking
                        .uncook(tensor.expression().as_view())
                        .parse_to_symbolic_net::<Aind>(&ParseSettings {
                            shorthand_parsing: ShorthandParsing::Expand {
                                schoonschip: SchoonschipExpansionMode {
                                    inner_products: false,
                                    expand_inside_chains: true,
                                    expand_schoonship: true,
                                },
                                trace: false,
                                chain: false,
                            },
                            ..Default::default()
                        })?
                        .simple_execute::<()>()?;
                    Ok(tensor.with_checked_expression(cooking.try_cook(explicit.as_view())?)?)
                };
                let probe = explicit_momenta(&normalize(&mass_control)?)?;
                let probe = probe.with_checked_expression(emitted_point(probe.expression())?)?;
                let probe_projection = projector.clone().multiply(probe.clone())?;
                Self::certify_parent_tensor_zero(
                    "D09 emitted Phi native mass identity keeps the independent unit calibration",
                    &probe_projection.try_sub(&projected_mass_control)?,
                )?;
                let project_complete_monomials = |label: &str, input: &Integrands| -> Result<_> {
                    // A captured emission can vanish under the actual finalizer.
                    // Its exact linear contribution is zero, with no Phi probe.
                    if input.iter().all(|(_, root)| root.is_zero()) {
                        return Ok(Vec::new());
                    }
                    let mut observations = Vec::new();
                    let mut controls = 0usize;
                    for (cut, root) in input.iter() {
                        let mut actual_calls = BTreeSet::new();
                        root.visitor(&mut |view| {
                            if let AtomView::Fun(call) = view
                                && call.get_symbol() == family
                            {
                                actual_calls.insert(view.to_owned());
                            }
                            true
                        });
                        let actual_calls = actual_calls.into_iter().collect::<Vec<_>>();
                        let rows = root.coefficient_list_exact(&actual_calls).ok_or_else(|| {
                            debug_tags!(#uv, #local, #four_d, #cff, #summary;
                                stage = "native_d09_emitted_phi_refused", label, file.cut = ?cut,
                                root_bytes = root.as_view().get_byte_size(), calls = actual_calls.len(),
                                file.root = %root.to_canonical_string(), file.actual_calls = ?actual_calls,
                                "Existing bounded complete-call polynomial owner refused; no fallback"
                            );
                            eyre!("D09 {label} complete-monomial Phi extraction refused")
                        })?;
                        let reconstructed = rows
                            .iter()
                            .map(|(monomial, coefficient)| monomial * coefficient)
                            .sum::<Atom>();
                        let residual = (root - &reconstructed)
                            .coefficient_list_exact(&actual_calls)
                            .ok_or_else(|| {
                                eyre!("D09 {label} complete-monomial reconstruction refused")
                            })?;
                        for (monomial, coefficient) in residual {
                            let coefficient = coefficient.collect_compact_factors();
                            eyre::ensure!(!coefficient.contains_symbol(family));
                            if coefficient.is_zero() {
                                continue;
                            }
                            // Verify distributivity only in the scalar CFF
                            // coefficient. Complete functions and powers stay
                            // opaque; no denominator or graph numerator opens.
                            let leaf = symbol!("native_d09_recomposition::leaf"; Scalar);
                            eyre::ensure!(!coefficient.contains_symbol(leaf));
                            let mut leaves = BTreeMap::new();
                            let mut exact = true;
                            let protected = coefficient.replace_map(|view, _, output| {
                                if let AtomView::Num(number) = view {
                                    exact &= matches!(
                                        number.get_coeff_view(),
                                        symbolica::coefficient::CoefficientView::Natural(..)
                                            | symbolica::coefficient::CoefficientView::Large(..)
                                    );
                                } else if matches!(
                                    view,
                                    AtomView::Var(_) | AtomView::Fun(_) | AtomView::Pow(_)
                                ) {
                                    let next = leaves.len();
                                    **output = leaves
                                        .entry(view.to_owned())
                                        .or_insert_with(|| {
                                            leaf.call_args([Atom::num(next), view.to_owned()])
                                        })
                                        .clone();
                                }
                            });
                            eyre::ensure!(
                                exact,
                                "D09 {label} reconstruction has an inexact arithmetic coefficient"
                            );
                            let restored = protected.replace_map(|view, _, output| {
                                if let AtomView::Fun(call) = view
                                    && call.get_symbol() == leaf
                                {
                                    output.set_from_view(&call.get(1));
                                }
                            });
                            eyre::ensure!(restored == coefficient);
                            let keys = leaves.into_values().collect::<Vec<_>>();
                            let proof = protected.coefficient_list_exact(&keys);
                            let proved = proof.as_ref().is_some_and(|rows| {
                                rows.iter().all(|(_, coefficient)| coefficient.is_zero())
                            });
                            if !proved {
                                debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                                    stage = "native_d09_scalar_recomposition_refused", label,
                                    file.monomial = %monomial.to_canonical_string(),
                                    file.coefficient = %coefficient.to_canonical_string(),
                                    file.protected = %protected.to_canonical_string(),
                                    leaves = keys.len(), collector_admitted = proof.is_some(),
                                    "Scalar-only exact ring certificate did not prove recomposition"
                                );
                            }
                            eyre::ensure!(
                                proved,
                                "D09 {label} exact scalar recomposition not established"
                            );
                        }
                        for (monomial, coefficient) in rows {
                            eyre::ensure!(!coefficient.contains_symbol(family));
                            if coefficient.is_zero() {
                                continue;
                            }
                            let factors = match monomial.as_view() {
                                AtomView::Fun(call) if call.get_symbol() == family => {
                                    vec![monomial.clone()]
                                }
                                AtomView::Mul(product) => product
                                    .iter()
                                    .map(|factor| factor.to_owned())
                                    .collect::<Vec<_>>(),
                                _ => {
                                    return Err(eyre!(
                                        "D09 {label} Phi refuses a scalar-only, inverse or powered family monomial"
                                    ));
                                }
                            };
                            eyre::ensure!(
                                (1..=2).contains(&factors.len()),
                                "D09 {label} Phi requires a proven complete one- or two-factor matrix"
                            );
                            let mut bodies = Vec::new();
                            let mut ordered_arguments = Vec::new();
                            for factor in &factors {
                                let AtomView::Fun(call) = factor.as_view() else {
                                    return Err(eyre!(
                                        "D09 {label} Phi found a powered or non-call factor"
                                    ));
                                };
                                eyre::ensure!(call.get_symbol() == family);
                                let entries = input
                                    .numerators()
                                    .iter()
                                    .filter(|entry| {
                                        call.get_nargs() == entry.tags.len() + entry.args.len()
                                            && entry.tags.iter().enumerate().all(|(index, tag)| {
                                                call.get(index) == tag.as_view()
                                            })
                                    })
                                    .collect::<Vec<_>>();
                                eyre::ensure!(
                                    entries.len() == 1,
                                    "D09 {label} Phi lacks one exact complete call binding"
                                );
                                let entry = entries[0];
                                let arguments = call
                                    .iter()
                                    .skip(entry.tags.len())
                                    .map(|argument| argument.to_owned())
                                    .collect::<Vec<_>>();
                                let replacements = entry
                                    .args
                                    .iter()
                                    .zip(&arguments)
                                    .map(|(formal, actual)| {
                                        Replacement::new(
                                            Atom::from(formal.clone()).to_pattern(),
                                            actual.clone(),
                                        )
                                    })
                                    .collect::<Vec<_>>();
                                eyre::ensure!(entry.lhs.replace_multiple(&replacements) == *factor);
                                bodies.push(entry.rhs.with_checked_expression(
                                    entry.rhs.expression().replace_multiple(&replacements),
                                )?);
                                ordered_arguments.push(
                                    arguments
                                        .iter()
                                        .map(&emitted_point)
                                        .collect::<Result<Vec<_>>>()?,
                                );
                            }
                            let mut bodies = bodies.into_iter();
                            let mut complete = bodies
                                .next()
                                .ok_or_else(|| eyre!("D09 Phi has no native matrix"))?;
                            for body in bodies {
                                complete = complete.multiply(body)?;
                            }
                            eyre::ensure!(
                                complete.structure().canonical().external_structure()
                                    == numerator.structure().canonical().external_structure(),
                                "D09 {label} complete family monomial has the wrong native parent interface: expected={:?}; actual={:?}",
                                numerator.structure().logical_slots(),
                                complete.structure().logical_slots(),
                            );
                            // Native multiplication connects written spin/color
                            // endpoints. It does not multiply scalar traces of
                            // the separate partial F/G families.
                            let complete = explicit_momenta(&normalize(&complete)?)?;
                            let complete = complete
                                .with_checked_expression(emitted_point(complete.expression())?)?;
                            let projected_raw = projector.clone().multiply(complete.clone())?;
                            eyre::ensure!(
                                projected_raw.structure().canonical().is_scalar(),
                                "D09 {label} whole-monomial Phi left tensor ports: projector={:?}; operand={:?}; result={:?}",
                                projector.structure().logical_slots(),
                                complete.structure().logical_slots(),
                                projected_raw.structure().logical_slots(),
                            );
                            // Close the complete matrix with the shared native
                            // Dirac/color owner before finite materialization.
                            // This preserves the graph's factorization instead
                            // of multiplying explicit gamma components first.
                            let projected = integrated::simplify(&projected_raw)?
                                .simplify_algebra(&idenso::tensor::AlgebraSettings {
                                    color: Some(
                                        ColorSimplifySettings::default()
                                            .with_cof_dimension_invariants(),
                                    ),
                                    ..Default::default()
                                })?;
                            let scalar = emitted_point(&on_shell_tests::scalar_component(
                                projected.expression(),
                            )?)?;
                            eyre::ensure!(!scalar.contains_symbol(family));
                            if controls == 0 {
                                let wrong =
                                    projector.clone().multiply(complete.try_add(&probe)?)?;
                                Self::certify_parent_tensor_zero(
                                    "D09 emitted whole-matrix mass shift leaves the independently calibrated unit",
                                    &wrong.try_sub(&projected_raw)?.try_sub(&probe_projection)?,
                                )?;
                                controls += 1;
                            }
                            let coefficient = emitted_point(&coefficient)?;
                            debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                                stage = "native_d09_emitted_complete_monomial_phi", label, file.cut = ?cut,
                                file.original_monomial = %monomial.to_canonical_string(),
                                file.ordered_pointed_arguments = ?ordered_arguments,
                                file.scalar_coefficient = %coefficient.to_canonical_string(),
                                file.phi_numerator = %scalar.to_canonical_string(),
                                native_input_bytes = projected_raw.expression().as_view().get_byte_size(),
                                native_reduced_bytes = projected.expression().as_view().get_byte_size(),
                                finite_scalar_bytes = scalar.as_view().get_byte_size(),
                                native_projector_once_per_complete_matrix = true,
                                mass_cancellation_proved = false, dimensional_R2_condition_proved = false,
                                "Native Phi acts once after actual ordered-call binding and complete matrix transport"
                            );
                            observations.push((
                                *cut,
                                monomial,
                                coefficient,
                                ordered_arguments,
                                scalar,
                            ));
                        }
                    }
                    eyre::ensure!(
                        !observations.is_empty() && controls == 1,
                        "D09 {label} Phi has no actual nonzero-coefficient matrix witness"
                    );
                    Ok(observations)
                };
                let source_phi = project_complete_monomials("actual source", &source_result)?;
                eyre::ensure!(
                    !source_phi.is_empty(),
                    "D09 source requires an actual nonzero-coefficient Phi witness"
                );
                let reduced_parent = source_graph
                    .full_filter()
                    .subtract(parent.subgraph())
                    .subtract(&source_graph.initial_state_cut);
                let parent_views = parent_assembly.finalized_views(
                    &source_graph,
                    &reduced_parent,
                    &allowed_keys,
                    &open_parent,
                )?;
                eyre::ensure!(parent_views.len() == parent_assembly.emissions.len());
                let mut parent_phi = Vec::new();
                for (ordinal, (emission, view)) in parent_assembly
                    .emissions
                    .iter()
                    .zip(&parent_views)
                    .enumerate()
                {
                    let first_row = parent_phi.len();
                    let rows = project_complete_monomials(
                        &format!("actual open parent emitted chunk {ordinal}"),
                        view,
                    )?;
                    let row_count = rows.len();
                    parent_phi.extend(rows);
                    debug_tags!(#uv, #local, #four_d, #cff, #summary;
                        stage = "native_d09_emitted_phi_partition", ordinal,
                        file.residue_key = ?emission.key,
                        first_row, row_count,
                        exact_zero_chunk = view.iter().all(|(_, root)| root.is_zero()),
                        exact_factorized_recomposition = true,
                        unchanged_coefficient_extraction_bound = true,
                        denominator_clearing = false,
                        "Linear Phi observes one actual finalized chunk and retains its complete key"
                    );
                }
                eyre::ensure!(
                    !parent_phi.is_empty(),
                    "D09 parent requires an actual nonzero-coefficient Phi witness across its emitted chunks"
                );
                // Diagnostic proof copies only. Keep every original Phi row,
                // call, ordered argument tuple and denominator coefficient.
                let normalize_phi_body = |scalar: &Atom| -> Result<Atom> {
                    // This is the evaluator's defining function map, applied
                    // only after every operator and point/chart substitution.
                    // Keep the actual tuples in their original compact form.
                    let lower_energy = |atom: &Atom| {
                        atom.replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
                            .with(Atom::var(W_.prop_).pow((1, 2)))
                    };
                    let normalized = lower_energy(&scalar_tail(scalar)).collect_compact_factors();
                    eyre::ensure!(
                        lower_energy(&scalar_tail(&normalized)).collect_compact_factors()
                            == normalized,
                        "D09 finite Phi scalar/even-energy normalization is not idempotent"
                    );
                    Ok(normalized)
                };
                let source_phi_bodies = source_phi
                    .iter()
                    .map(|row| normalize_phi_body(&row.4))
                    .collect::<Result<Vec<_>>>()?;
                let parent_phi_bodies = parent_phi
                    .iter()
                    .map(|row| normalize_phi_body(&row.4))
                    .collect::<Result<Vec<_>>>()?;
                if given.subgraph().is_empty() {
                    // Each actual row contributes once. Pairwise alias witnesses
                    // are evidence, never extra summands in a complete bipartite class.
                    let mut coefficient_classes = BTreeMap::<
                        (crate::cff::CutCFFIndex, Atom),
                        Vec<(bool, usize, i64, Atom)>,
                    >::new();
                    let mut covered_rows = BTreeSet::new();
                    let mut explicit_zero_rows = 0usize;
                    for (is_source, rows, bodies) in [
                        (true, source_phi.as_slice(), source_phi_bodies.as_slice()),
                        (false, parent_phi.as_slice(), parent_phi_bodies.as_slice()),
                    ] {
                        eyre::ensure!(rows.len() == bodies.len());
                        for (index, (row, body)) in rows.iter().zip(bodies).enumerate() {
                            eyre::ensure!(
                                covered_rows.insert((is_source, index)),
                                "D09 Phi coefficient class reused an actual row"
                            );
                            if body.is_zero() {
                                explicit_zero_rows += 1;
                                continue;
                            }
                            let negative = normalize_phi_body(&(-body))?;
                            let (representative, sign) = if body <= &negative {
                                (body.clone(), 1i64)
                            } else {
                                (negative, -1i64)
                            };
                            eyre::ensure!(
                                normalize_phi_body(&(body - Atom::num(sign) * &representative))?
                                    .is_zero(),
                                "D09 Phi row lacks its exact signed opaque-body binding"
                            );
                            let coefficient = Self::normalized_d09_phi_coefficient(&row.2)?;
                            coefficient_classes
                                .entry((row.0, representative))
                                .or_default()
                                .push((is_source, index, sign, coefficient));
                        }
                    }
                    let class_members = coefficient_classes.values().map(Vec::len).sum::<usize>();
                    eyre::ensure!(
                        !coefficient_classes.is_empty()
                            && covered_rows.len() == source_phi.len() + parent_phi.len()
                            && explicit_zero_rows + class_members == covered_rows.len(),
                        "D09 Phi coefficient classes do not cover every actual row exactly once"
                    );
                    let class_count = coefficient_classes.len();
                    for (class, ((cut, _body), members)) in
                        coefficient_classes.into_iter().enumerate()
                    {
                        // The retained body is only a key; all algebra below acts
                        // on its scalar CFF coefficients, never on that Phi body.
                        let sum = Atom::add_many(
                            members
                                .iter()
                                .map(|(_, _, sign, coefficient)| Atom::num(*sign) * coefficient),
                        );
                        let residual = Self::normalized_d09_phi_coefficient(&sum)?;
                        eyre::ensure!(
                            residual.is_zero(),
                            "D09 signed Phi coefficient class {class} did not cancel"
                        );
                        let wrong = Self::normalized_d09_phi_coefficient(&(sum + Atom::one()))?;
                        eyre::ensure!(
                            wrong.is_one(),
                            "D09 scalar class +1 perturbation must leave exactly one"
                        );
                        debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                            stage = "native_d09_emitted_mass_coefficient_class", class,
                            given = %given.subgraph().string_label(), file.cut = ?cut,
                            file.members = ?members.iter().map(|(source, index, sign, _)|
                                (*source, *index, *sign)).collect::<Vec<_>>(),
                            scalar_coefficient_zero = true, plus_one_residual_exactly_one = true,
                            phi_body_entered_scalar_algebra = false,
                            "Each actual signed opaque-body class has an exact zero scalar weight"
                        );
                    }
                    debug_tags!(#uv, #local, #four_d, #cff, #summary;
                        stage = "native_d09_emitted_mass_coefficient_identity",
                        given = %given.subgraph().string_label(),
                        actual_rows = covered_rows.len(), explicit_zero_rows, class_members, class_count,
                        every_actual_row_once = true, signed_class_coefficients_zero = true,
                        unchanged_rational_preflight_budget = true, scalar_round_limit = 5,
                        plus_one_nonzero_controls = true, phi_body_entered_scalar_algebra = false,
                        emitted_mass_cancellation_proved = true,
                        dimensional_R2_condition_proved = false,
                        "Actual emitted CFF Phi rows cancel by bounded scalar coefficient identities"
                    );
                } else {
                    use spenso::shadowing::Collectable;
                    struct LinearGroup {
                        cut: crate::cff::CutCFFIndex,
                        coefficient: Atom,
                        body: Atom,
                        members: Vec<(bool, usize)>,
                    }
                    let mut groups = BTreeMap::<_, (Vec<(bool, usize)>, Vec<Atom>)>::new();
                    let mut covered = BTreeSet::new();
                    let mut closed = BTreeSet::new();
                    let mut literal_zero_rows = 0usize;
                    let mut coefficient_zero_rows = 0usize;
                    let mut radicands = BTreeSet::new();
                    for (source, rows, bodies) in [
                        (true, source_phi.as_slice(), source_phi_bodies.as_slice()),
                        (false, parent_phi.as_slice(), parent_phi_bodies.as_slice()),
                    ] {
                        eyre::ensure!(rows.len() == bodies.len());
                        for (index, (row, body)) in rows.iter().zip(bodies).enumerate() {
                            eyre::ensure!(
                                covered.insert((source, index)),
                                "D09 linear proof reused an actual row"
                            );
                            if source {
                                row.2.visitor(&mut |part| {
                                    match part {
                                        AtomView::Fun(call)
                                            if call.get_symbol() == GS.on_shell_energy
                                                && call.get_nargs() == 2 =>
                                        {
                                            radicands.insert(call.get(1).to_owned());
                                        }
                                        AtomView::Pow(power)
                                            if i64::try_from(power.get_exp()).is_err() =>
                                        {
                                            radicands.insert(power.get_base().to_owned());
                                        }
                                        _ => {}
                                    }
                                    true
                                });
                            }
                            if body.is_zero() {
                                literal_zero_rows += 1;
                                eyre::ensure!(closed.insert((source, index)));
                                continue;
                            }
                            let coefficient = Self::normalized_d09_phi_coefficient(&row.2)?;
                            if coefficient.is_zero() {
                                coefficient_zero_rows += 1;
                                eyre::ensure!(closed.insert((source, index)));
                                continue;
                            }
                            let numeric = coefficient.as_view().collect_numeric_content();
                            let number = match numeric.as_view() {
                                AtomView::Mul(product) => product
                                    .get_coefficient()
                                    .map(|n| n.to_owned())
                                    .unwrap_or_else(Atom::one),
                                AtomView::Num(_) => numeric.clone(),
                                _ => Atom::one(),
                            };
                            let key = &numeric / &number;
                            eyre::ensure!(
                                Self::normalized_d09_phi_coefficient(
                                    &(&coefficient - &number * &key)
                                )?
                                .is_zero()
                            );
                            let (members, weighted_bodies) =
                                groups.entry((row.0, key)).or_default();
                            members.push((source, index));
                            weighted_bodies.push(number * body);
                        }
                    }
                    eyre::ensure!(covered.len() == source_phi.len() + parent_phi.len());
                    // The original component inventory already binds every global
                    // coordinate to its actual source/parent domain and full keys.
                    eyre::ensure!(
                        component_coordinates.len() == 2,
                        "D09 linear proof requires its two certified component coordinates"
                    );
                    let spatial = |edge: EdgeIndex| {
                        (1..=3)
                            .map(|i| function!(GS.emr_mom, edge.0, GS.cind(i)))
                            .collect::<Vec<_>>()
                    };
                    let mut charts = Vec::new();
                    for hard in &component_coordinates {
                        for soft in &component_coordinates {
                            if hard == soft
                                || !given.lmb().loop_edges.contains(hard)
                                || !parent.lmb().loop_edges.contains(soft)
                                || given.lmb().loop_edges.contains(soft)
                            {
                                continue;
                            }
                            let k = spatial(*hard);
                            let l = spatial(*soft);
                            let hard_base = Atom::add_many(k.iter().map(|x| x.pow(2)))
                                + Atom::var(GS.m_uv_vacuum).pow(2);
                            let soft_base = Atom::add_many(l.iter().map(|x| x.pow(2)));
                            if [&hard_base, &soft_base, &(&soft_base + Atom::one())]
                                .into_iter()
                                .all(|base| radicands.contains(base))
                            {
                                charts.push((*hard, *soft, k, l, hard_base, soft_base));
                            }
                        }
                    }
                    eyre::ensure!(
                        charts.len() == 1,
                        "D09 coefficient radicands do not select one actual component chart: candidates={}",
                        charts.len()
                    );
                    let (hard, soft, k, l, hard_base, soft_base) = charts.pop().unwrap();
                    let dot = Atom::add_many(k.iter().zip(&l).map(|(x, y)| x * y));
                    let names = ["A", "B", "C", "H", "K", "L", "t"];
                    let symbols =
                        names.map(|name| symbol!(&format!("native_d09_phi_chart::{name}"); Scalar));
                    for head in symbols {
                        eyre::ensure!(
                            source_phi
                                .iter()
                                .chain(&parent_phi)
                                .all(|row| !row.2.contains_symbol(head)
                                    && !row.4.contains_symbol(head)),
                            "D09 chart symbol collides with an actual row"
                        );
                    }
                    let [a, b, c, h, big_k, big_l, t] = symbols.map(Atom::var);
                    let actual_b = soft_base.pow((1, 2));
                    let actual_c = (&soft_base + Atom::one()).pow((1, 2));
                    eyre::ensure!(
                        (actual_c.pow(2) - actual_b.pow(2))
                            .collect_compact_factors()
                            .is_one(),
                        "D09 physical energy roots do not have unit squared difference"
                    );
                    eyre::ensure!(
                        (&actual_c * actual_c.pow(-1)).is_one(),
                        "D09 physical massive energy is not a certified unit"
                    );
                    let definitions = [
                        (hard_base.pow((1, 2)), a.clone()),
                        (actual_b.clone(), b.clone()),
                        (actual_c.clone(), c.clone()),
                        (hard_base.clone(), big_k.clone()),
                        (soft_base.clone(), big_l.clone()),
                        (dot, h.clone()),
                    ];
                    let restore = |atom: &Atom| {
                        let mut restored = atom.clone();
                        for (original, alias) in definitions.iter().rev() {
                            restored = restored
                                .replace(alias.to_pattern())
                                .with(original.to_pattern());
                        }
                        restored
                    };
                    let chart_view = |atom: &Atom| -> Result<Atom> {
                        eyre::ensure!(
                            atom.as_view().get_byte_size() <= 32768,
                            "D09 factorized chart input exceeds 32KiB"
                        );
                        let mut view = atom.clone();
                        for (original, alias) in &definitions {
                            view = view.replace(original.to_pattern()).with(alias.to_pattern());
                        }
                        view = view
                            .replace(Atom::var(GS.m_uv_vacuum).pow(2).to_pattern())
                            .with(
                                (&big_k + Atom::add_many(k.iter().map(|x| -x.pow(2)))).to_pattern(),
                            )
                            .replace((&k[2] * &l[2]).to_pattern())
                            .with((&h - &k[0] * &l[0] - &k[1] * &l[1]).to_pattern())
                            .replace(l[2].pow(2).to_pattern())
                            .with((&big_l - l[0].pow(2) - l[1].pow(2)).to_pattern());
                        eyre::ensure!(
                            restore(&view) == *atom,
                            "D09 affine chart does not literally restore its complete input"
                        );
                        view = view
                            .replace(symbols[4])
                            .with(a.pow(2))
                            .replace(symbols[5])
                            .with(c.pow(2) - Atom::one());
                        eyre::ensure!(
                            restore(&view) == *atom,
                            "D09 root chart does not literally restore its complete input"
                        );
                        Ok(view)
                    };
                    let compact_chart = |atom: &Atom| -> Result<Atom> {
                        let mut current = atom.clone();
                        for _ in 0..8 {
                            let next = current
                                .replace_map_bottom_up(|part, context, output| {
                                    if context.function_level == 0
                                        && matches!(part, AtomView::Add(_) | AtomView::Mul(_))
                                    {
                                        **output = part.collect_compact_factors();
                                    }
                                })
                                .collect_compact_factors();
                            eyre::ensure!(
                                next.as_view().get_byte_size() <= 32768,
                                "D09 protected compact chart exceeded its byte scope"
                            );
                            if next == current
                                || next.as_view().get_byte_size()
                                    >= current.as_view().get_byte_size()
                            {
                                return Ok(current);
                            }
                            current = next;
                        }
                        Ok(current)
                    };
                    let mut unresolved = Vec::new();
                    let mut direct_groups = 0usize;
                    let mut invariant_groups = 0usize;
                    for ((cut, coefficient), (members, bodies)) in groups {
                        let original_body = Atom::add_many(bodies);
                        let body = original_body.collect_compact_factors();
                        let direct = body.is_zero();
                        let invariant = if direct {
                            false
                        } else {
                            compact_chart(&chart_view(&body)?)?.is_zero()
                        };
                        if direct || invariant {
                            let original_control = &original_body + Atom::one();
                            let mut original_control_one = original_control
                                .collect_compact_factors()
                                .is_one()
                                || (invariant
                                    && compact_chart(&chart_view(&original_control)?)?.is_one());
                            let mut coupling_control_rows = 0usize;
                            if !original_control_one {
                                // Only the actual registered coupling is a column.
                                // Every coupling-free complete graph factor stays
                                // an opaque AtomField coefficient in this owner.
                                let mut coupling_keys = graph
                                    .model
                                    .couplings()
                                    .iter()
                                    .map(|coupling| {
                                        Atom::var(symbol!(&format!("UFO::{}", coupling.name)))
                                    })
                                    .filter(|coupling| {
                                        original_control.contains(coupling.as_view())
                                    })
                                    .collect::<Vec<_>>();
                                eyre::ensure!(
                                    coupling_keys.len() == 1,
                                    "D09 original+1 requires one actual registered coupling"
                                );
                                let coupling = coupling_keys.pop().unwrap();
                                let keys = std::slice::from_ref(&coupling);
                                let rows = original_control.coefficient_list_exact(keys)
                                    .ok_or_else(|| eyre!("D09 original+1 coupling collection refused its unchanged budget"))?;
                                let recomposed = Atom::add_many(
                                    rows.iter()
                                        .map(|(monomial, coefficient)| monomial * coefficient),
                                );
                                eyre::ensure!(
                                    recomposed.coefficient_list_exact(keys).as_ref() == Some(&rows),
                                    "D09 original+1 coupling rows do not recompose in their exact owner"
                                );
                                let difference = (&recomposed - &original_control).coefficient_list_exact(keys)
                                    .ok_or_else(|| eyre!("D09 original+1 coupling reconstruction refused its unchanged budget"))?;
                                eyre::ensure!(
                                    difference.iter().all(|(_, coefficient)| coefficient
                                        .collect_compact_factors()
                                        .is_zero()),
                                    "D09 original+1 coupling recomposition changed an original complete coefficient"
                                );
                                let mut constants = 0usize;
                                for (monomial, coefficient) in &rows {
                                    eyre::ensure!(
                                        !coefficient.contains(coupling.as_view()),
                                        "D09 coupling coefficient retained its collection key"
                                    );
                                    if monomial.is_one() {
                                        eyre::ensure!(
                                            coefficient.is_one(),
                                            "D09 original+1 constant coefficient is not literally one"
                                        );
                                        constants += 1;
                                    } else {
                                        eyre::ensure!(
                                            *monomial == coupling.pow(4),
                                            "D09 original+1 has an unsupported coupling degree"
                                        );
                                        eyre::ensure!(
                                            coefficient.collect_compact_factors().is_zero()
                                                || (invariant
                                                    && compact_chart(&chart_view(coefficient)?)?
                                                        .is_zero()),
                                            "D09 original+1 retains a nonzero complete coupling coefficient"
                                        );
                                    }
                                    debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                                        stage = "native_d09_phi_original_unit_coupling_control",
                                        file.members = ?members, file.actual_coupling = %coupling.to_canonical_string(),
                                        file.monomial = %monomial.to_canonical_string(), file.complete_coefficient = %coefficient.to_canonical_string(),
                                        constant_coefficient_literal_one = monomial.is_one(),
                                        nonconstant_complete_coefficient_zero = !monomial.is_one(),
                                        exact_original_control_recomposition = true,
                                        only_actual_coupling_is_polynomial_column = true,
                                        "The original summed Phi expression plus one has one constant coefficient and zero complete coupling coefficients"
                                    );
                                }
                                eyre::ensure!(
                                    constants == 1,
                                    "D09 original+1 did not retain one unit constant row"
                                );
                                coupling_control_rows = rows.len();
                                original_control_one = true;
                            }
                            eyre::ensure!(
                                original_control_one,
                                "D09 closed linear group lost its original+1 control"
                            );
                            for member in &members {
                                eyre::ensure!(
                                    closed.insert(*member),
                                    "D09 linear group closed a row twice"
                                );
                            }
                            direct_groups += usize::from(direct);
                            invariant_groups += usize::from(invariant);
                            debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                                stage = "native_d09_phi_linear_group", file.cut = ?cut,
                                file.coefficient = %coefficient.to_canonical_string(), file.members = ?members,
                                protected_compact_zero = direct, reversible_invariant_zero = invariant,
                                original_plus_one_exact_one = true, coupling_control_rows,
                                "One actual scalar-coefficient group closes by an exact factorized body identity"
                            );
                        } else {
                            unresolved.push(LinearGroup {
                                cut,
                                coefficient,
                                body,
                                members,
                            });
                        }
                    }
                    eyre::ensure!(
                        unresolved.len() == 6,
                        "D09 supported completed chart requires six unresolved scalar groups; actual={}",
                        unresolved.len()
                    );
                    struct Triple {
                        reference: usize,
                        sign: i64,
                        slots: [usize; 3],
                        weights: [Atom; 3],
                    }
                    let mut candidates = Vec::new();
                    for reference in 0..unresolved.len() {
                        for sign in [1i64, -1] {
                            let weights = [
                                actual_c.clone(),
                                Atom::one() + Atom::num(sign) * &actual_b,
                                Atom::one()
                                    + Atom::num(sign) * &actual_b
                                    + Atom::num(sign) * &actual_c,
                            ];
                            let mut matches = [Vec::new(), Vec::new(), Vec::new()];
                            for (slot, weight) in weights.iter().enumerate() {
                                for (index, group) in unresolved.iter().enumerate() {
                                    if group.cut != unresolved[reference].cut {
                                        continue;
                                    }
                                    // Only actual scalar weights enter this unchanged
                                    // rational certificate; no division by c_ref.
                                    let residual = Self::normalized_d09_phi_coefficient(
                                        &(&actual_c * &group.coefficient
                                            - weight * &unresolved[reference].coefficient),
                                    )?;
                                    if residual.is_zero() {
                                        matches[slot].push(index);
                                    }
                                }
                            }
                            eyre::ensure!(
                                matches.iter().all(|slot| slot.len() <= 1),
                                "D09 coefficient pattern has ambiguous group membership"
                            );
                            if matches.iter().all(|slot| slot.len() == 1) {
                                let slots = std::array::from_fn(|i| matches[i][0]);
                                if slots.into_iter().collect::<BTreeSet<_>>().len() == 3 {
                                    candidates.push(Triple {
                                        reference,
                                        sign,
                                        slots,
                                        weights,
                                    });
                                }
                            }
                        }
                    }
                    let mut covers = Vec::new();
                    for first in 0..candidates.len() {
                        for second in first + 1..candidates.len() {
                            let mut members = BTreeSet::new();
                            if candidates[first]
                                .slots
                                .iter()
                                .chain(&candidates[second].slots)
                                .all(|slot| members.insert(*slot))
                                && members.len() == unresolved.len()
                            {
                                covers.push([first, second]);
                            }
                        }
                    }
                    eyre::ensure!(
                        covers.len() == 1,
                        "D09 scalar coefficients do not give one exact two-triple cover: {}",
                        covers.len()
                    );
                    let tb = (&t - t.pow(-1)) / Atom::num(2);
                    let tc = (&t + t.pow(-1)) / Atom::num(2);
                    eyre::ensure!((&tc + &tb - &t).collect_compact_factors().is_zero());
                    eyre::ensure!((&tc - &tb - t.pow(-1)).collect_compact_factors().is_zero());
                    eyre::ensure!((&t * t.pow(-1)).is_one());
                    eyre::ensure!(
                        ((((&c + &b) - (&c - &b)) / Atom::num(2)) - &b)
                            .collect_compact_factors()
                            .is_zero()
                    );
                    eyre::ensure!(
                        ((((&c + &b) + (&c - &b)) / Atom::num(2)) - &c)
                            .collect_compact_factors()
                            .is_zero()
                    );
                    let mut grid_nodes = 0usize;
                    let mut triple_members = 0usize;
                    for candidate in covers.pop().unwrap() {
                        let triple = &candidates[candidate];
                        let mut terms = Vec::new();
                        for (slot, weight) in triple.slots.iter().zip(&triple.weights) {
                            let group = &unresolved[*slot];
                            eyre::ensure!(
                                Self::normalized_d09_phi_coefficient(
                                    &(&actual_c * &group.coefficient
                                        - weight * &unresolved[triple.reference].coefficient)
                                )?
                                .is_zero()
                            );
                            // Bind each of the three complete weighted bodies before
                            // compaction; never distribute a product or power.
                            terms.push(compact_chart(&chart_view(&(weight * &group.body))?)?);
                            debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                                stage = "native_d09_phi_coefficient_linear_relation", file.cut = ?group.cut,
                                file.members = ?group.members, file.coefficient = %group.coefficient.to_canonical_string(),
                                file.reference_coefficient = %unresolved[triple.reference].coefficient.to_canonical_string(),
                                file.weight = %weight.to_canonical_string(), file.unit_energy = %actual_c.to_canonical_string(),
                                scalar_cross_residual_zero = true, no_phi_in_scalar_certificate = true,
                                "The original CFF weight is bound to one known-energy-unit linear relation"
                            );
                        }
                        let original_chart = Atom::add_many(terms);
                        let laurent = original_chart
                            .replace(symbols[1])
                            .with(tb.clone())
                            .replace(symbols[2])
                            .with(tc.clone());
                        let mut free = BTreeSet::new();
                        laurent.visitor(&mut |part| {
                            if let AtomView::Var(variable) = part {
                                free.insert(variable.get_symbol());
                            }
                            true
                        });
                        for head in [symbols[0], symbols[3], symbols[6]] {
                            free.remove(&head);
                        }
                        eyre::ensure!(
                            free.len() == 1,
                            "D09 exact chart requires one actual remaining coupling variable; actual={free:?}"
                        );
                        let coupling = *free.iter().next().unwrap();
                        grid_nodes += Self::certify_d09_phi_laurent_grid(
                            &laurent,
                            [symbols[0], symbols[3], coupling, symbols[6]],
                        )?;
                        for slot in triple.slots {
                            for member in &unresolved[slot].members {
                                eyre::ensure!(
                                    closed.insert(*member),
                                    "D09 exact triple reused an actual row"
                                );
                                triple_members += 1;
                            }
                        }
                        debug_tags!(#uv, #local, #four_d, #cff, #summary;
                            stage = "native_d09_phi_exact_linear_triple", sign = triple.sign,
                            file.group_indices = ?triple.slots, file.hard_coordinate = ?hard, file.soft_coordinate = ?soft,
                            file.hard_radicand = %hard_base.to_canonical_string(), file.soft_radicand = %soft_base.to_canonical_string(),
                            physical_mass_unit = true, actual_unit_energy_positive = true,
                            actual_given_hard_coordinate = true, actual_parent_retained_coordinate = true,
                            reversible_whole_term_bindings = true, actual_C_squared_minus_B_squared_one = true,
                            rational_chart_t_positive = true, emitted_mass_linear_identity_proved = true,
                            dimensional_R2_condition_proved = false,
                            "The actual completed D09 coefficient triple closes on a full exact Laurent grid"
                        );
                    }
                    eyre::ensure!(
                        closed == covered,
                        "D09 linear proof left an actual emitted row uncovered"
                    );
                    debug_tags!(#uv, #local, #four_d, #cff, #summary;
                        stage = "native_d09_emitted_mass_linear_identity", given = %given.subgraph().string_label(),
                        actual_rows = covered.len(), literal_zero_rows, coefficient_zero_rows,
                        direct_groups, invariant_groups, exact_triples = 2, triple_members, grid_nodes,
                        every_actual_row_once = true, scalar_coefficient_relations_exact = true,
                        unchanged_rational_preflight_budget = true, original_chart_plus_one_controls = true,
                        direct_phi_polynomial_conversion = false, graph_numerator_factorization_preserved = true,
                        emitted_mass_cancellation_proved = true, dimensional_R2_condition_proved = false,
                        "Every actual completed D09 row is covered by an exact factorized linear proof"
                    );
                }
                let mut source_body_aliases =
                    BTreeMap::<(crate::cff::CutCFFIndex, Atom), Vec<(usize, i64)>>::new();
                for (index, (row, body)) in source_phi.iter().zip(&source_phi_bodies).enumerate() {
                    if body.is_zero() {
                        continue;
                    }
                    for sign in [1_i64, -1] {
                        let signed = (Atom::num(sign) * body).collect_compact_factors();
                        source_body_aliases
                            .entry((row.0, signed))
                            .or_default()
                            .push((index, sign));
                    }
                }
                let mut phi_alias_pairs = 0usize;
                let mut parent_alias_rows = 0usize;
                let mut parent_zero_rows = 0usize;
                for (parent_index, (parent_row, parent_body)) in
                    parent_phi.iter().zip(&parent_phi_bodies).enumerate()
                {
                    if parent_body.is_zero() {
                        parent_zero_rows += 1;
                        continue;
                    }
                    let Some(aliases) =
                        source_body_aliases.get(&(parent_row.0, parent_body.clone()))
                    else {
                        continue;
                    };
                    parent_alias_rows += 1;
                    for (source_index, sign) in aliases {
                        let source_row = &source_phi[*source_index];
                        let residual = normalize_phi_body(
                            &(parent_body - Atom::num(*sign) * &source_phi_bodies[*source_index]),
                        )?;
                        eyre::ensure!(residual.is_zero());
                        phi_alias_pairs += 1;
                        debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                            stage = "native_d09_finite_phi_body_identity",
                            given = %given.subgraph().string_label(), file.cut = ?parent_row.0,
                            source_index = *source_index, parent_index, parent_equals_source_times = *sign,
                            file.source_monomial = %source_row.1.to_canonical_string(),
                            file.parent_monomial = %parent_row.1.to_canonical_string(),
                            file.source_pointed_arguments = ?source_row.3,
                            file.parent_pointed_arguments = ?parent_row.3,
                            file.common_phi_body = %parent_body.to_canonical_string(),
                            file.source_coefficient = %source_row.2.to_canonical_string(),
                            file.parent_coefficient = %parent_row.2.to_canonical_string(),
                            exact_finite_phi_scalar_identity = true,
                            full_key_domain_energy_association_checked = false,
                            denominator_coefficient_relation_checked = false,
                            emitted_mass_cancellation_proved = true, dimensional_R2_condition_proved = false,
                            "Finite complete Phi scalars agree up to an exact sign; denominator weights remain separate"
                        );
                    }
                }
                debug_tags!(#uv, #local, #four_d, #cff, #summary;
                    stage = "native_d09_finite_phi_alias_diagnostic",
                    given = %given.subgraph().string_label(),
                    source_exact_zero_bodies = source_phi_bodies.iter().filter(|body| body.is_zero()).count(),
                    parent_exact_zero_bodies = parent_zero_rows,
                    parent_exact_signed_body_matches = parent_alias_rows, exact_body_identity_pairs = phi_alias_pairs,
                    parent_unmatched_bodies = parent_phi_bodies.len() - parent_zero_rows - parent_alias_rows,
                    coefficients_combined = false, energy_algebra_applied_to_phi = false,
                    emitted_mass_cancellation_proved = true, dimensional_R2_condition_proved = false,
                    "Body aliases supplement the separate once-only summed cancellation certificate"
                );
                debug_tags!(#uv, #local, #four_d, #cff, #summary;
                    stage = "native_d09_emitted_phi_transport", given = %given.subgraph().string_label(),
                    source_monomials = source_phi.len(), parent_monomials = parent_phi.len(),
                    unchanged_component_coordinates = component_coordinates.len(), frozen_kernels = frozen_kernels.len(),
                    actual_complete_monomial_phi = true, independent_native_mass_probe = true,
                    emitted_mass_cancellation_proved = true, dimensional_R2_condition_proved = false,
                    "Both actual CFF lanes obey the exact summed mass condition with factorized Phi observations"
                );
                eyre::ensure!(prefix == &retained_prefix);
                eyre::ensure!(local == &retained_local);
                eyre::ensure!(integrated.physical_finite_counterterm()? == retained_finite);
                eyre::ensure!(
                    projection_graph.global_prefactor.projector == retained_global_projector
                );
                debug_tags!(#uv, #local, #four_d, #cff, #summary;
                    stage = "native_d09_actual_source_cff", given = %given.subgraph().string_label(),
                    families = source_result.numerators().len(),
                    open_parent_families = open_parent.numerators().len(),
                    root_bytes = source_result.iter().map(|(_, root)| root.as_view().get_byte_size()).sum::<usize>(),
                    body_bytes = source_result.numerators().iter().map(|entry| entry.rhs.expression().as_view().get_byte_size()).sum::<usize>(),
                    actual_source_constructed = true, projection_only_open_external_matrix = true,
                    analytic_backend_repeated = false, emitted_mass_cancellation_proved = true,
                    "Actual D09 source and parent CFF pass the matched emitted mass condition"
                );
            }
        }
        debug_tags!(#uv, #integrated, #summary;
            stage = "native_d09_rational_finalizer_mass_projection",
            source = %diagram.name(), given = %given.subgraph().string_label(),
            domains = domains.len(), frozen_domains,
            per_full_domain_source_local_zero = true, actual_finalized_finite_mass_zero = true,
            mass_projector_unit_calibration = true, mass_perturbation_nonzero_control = true,
            signed_parent_carrier_transport = true, reconstructed_kernels_nonzero = true,
            child_epsilon_depth_checked_separately = true,
            emitted_cff_condition_proved = true, dimensional_R2_condition_proved = false,
            "Actual D09 operands pass the rational and emitted-CFF matched mass conditions"
        );
        Ok(())
    }

    pub(in crate::uv) fn certify_parent_mass_projection_input(
        &self,
        diagram: &feynkit_graph::FeynmanDiagram,
        ctx: &UVCtx<'_>,
        parent: &OwnedForestNode,
        given: &OwnedForestNode,
        prefix: &Full4dCts,
    ) -> Result<Local4dCts> {
        let reference = self;
        let graph = ctx.graph;
        let child = given;
        eyre::ensure!(
            parent.renormalization_scheme() == ApproximationType::OS && parent.dod() == 1
        );
        eyre::ensure!(graph.n_loops(parent.subgraph()) == 2);
        eyre::ensure!(
            graph.paired_edges(parent.subgraph()) == graph.paired_edges(&graph.full_filter())
        );
        eyre::ensure!(ctx.settings.generate_integrated && !ctx.settings.add_marker);
        eyre::ensure!(child.subgraph().is_empty() || graph.n_loops(child.subgraph()) == 1);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let exact_zero = Self::certify_parent_tensor_zero;
        let frame_zero =
            |label: &str, difference: &SymbolicTensor<IndexedStructure>| -> Result<()> {
                exact_zero(
                    label,
                    &Self::reduce_parent_external_frame(label, difference)?,
                )
            };
        assert_eq!(GS.emr_mom, feynkit_graph::symbols::momentum());
        assert_eq!(GS.hedgeaind, feynkit_graph::symbols::hedge_index());
        let native_projector = diagram.projector().clone().try_into_indexed()?;
        assert!(
            graph
                .global_prefactor
                .projector
                .same_value(&native_projector)
        );
        let mut native_ports = Vec::new();
        let projector_source = cooking.uncook(diagram.projector().expression().as_view());
        projector_source.visitor(&mut |view| {
            if let AtomView::Fun(wavefunction) = view
                && [GS.u, GS.ubar].contains(&wavefunction.get_symbol())
                && wavefunction.get_nargs() == 2
                && let Ok(edge) = usize::try_from(wavefunction.get(0))
            {
                native_ports.push((EdgeIndex(edge), wavefunction.get(1).to_owned()));
            }
            true
        });
        assert_eq!(native_ports.len(), 2);
        for (edge, native_port) in native_ports {
            let runtime_ports = graph.underlying[graph[&edge].1.any_hedge()]
                .edge_spin_slots()
                .map(|slot| slot.to_atom())
                .collect::<Vec<_>>();
            assert_eq!(runtime_ports, [native_port]);
        }
        let crown = graph.dummy_less_full_crown(parent.subgraph());
        let boundary = crown.included_iter().collect::<Vec<_>>();
        assert_eq!(boundary.len(), 2);
        let incoming = *boundary
            .iter()
            .find(|h| graph.flow(**h) == Flow::Sink)
            .unwrap();
        let owner = graph[&incoming];
        let mass = graph[owner].particle.mass_atom(&graph.model);
        let external = boundary.iter().map(|h| graph[h]).collect::<Vec<_>>();
        let (projector, _) = reference.mass_projector(diagram, &Atom::var(GS.dim))?;
        // The generic reference promotes its source color ports to Nc.
        // This runtime source instead retains its native SU(3) representations.
        let projector = projector
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
        let projector =
            projector.with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
        let projector = projector.try_into_indexed()?;
        let velocity = ctx.settings.on_shell_velocity()?;
        let lmb = compatible_lmb(ctx, parent, child, prefix)?;
        let source = grow(prefix, ctx, parent, child)?;
        let bare_model_source = graph
            .numerator(parent.subgraph(), &graph.empty_subgraph())?
            .to_d_dim(GS.dim)?
            .state
            .expr;
        eyre::ensure!(
            !bare_model_source.expression().contains_symbol(AGS.gamma0)
                && !bare_model_source.expression().contains_symbol(ETS.delta),
            "the canonical model itself has no fixed subtraction frame"
        );
        if child.renormalization_scheme() == ApproximationType::OS {
            eyre::ensure!(graph.n_loops(child.subgraph()) == 1 && child.dod() == 1);
            eyre::ensure!(prefix.0.has_on_shell_ancestry);
            eyre::ensure!(
                velocity == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero],
                "continuation of the certified primitive OS child is initially rest-only"
            );
        } else {
            eyre::ensure!(
                !source.expression().contains_symbol(AGS.gamma0)
                    && !source.expression().contains_symbol(ETS.delta),
                "a non-OS prefix must not introduce an external fixed frame"
            );
        }
        let at_point = |tensor: &SymbolicTensor<IndexedStructure>,
                        sign: i32|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let atom = tensor.expression();
            let point = GS.on_shell_momentum(owner, &(Atom::num(sign) * &mass), &velocity);
            let routed = GS.erase_uv_momentum_provenance(atom).replace_multiple(
                graph.uv_wrapped_replacement(parent.subgraph(), &lmb, &[W_.x___]),
            );
            let expression = GS
                .resolve_os_momenta(&routed.replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() > 0
                        && usize::try_from(momentum.get(0))
                            .is_ok_and(|edge| external.contains(&EdgeIndex(edge)))
                    {
                        let indices = momentum
                            .iter()
                            .skip(1)
                            .map(|index| index.to_owned())
                            .collect::<Vec<_>>();
                        **output = GS.indexed_momentum(&point, &indices);
                    }
                }))
                .replace(mass.to_pattern())
                .with(1);
            Ok(tensor.with_checked_expression(expression)?)
        };
        let projector = at_point(&projector, 1)?;
        let phi = |tensor: &SymbolicTensor<IndexedStructure>| -> Result<SymbolicTensor<IndexedStructure>> {
            Ok(projector.clone().multiply(at_point(tensor, 1)?)?)
        };
        let [positive, negative] =
            graph.on_shell_branches(&graph.model, &source, parent.subgraph(), &lmb, &velocity)?;
        for (label, value) in [
            ("source", phi(&source)?),
            ("positive", phi(&positive)?),
            ("negative", phi(&negative)?),
        ] {
            let admitted = &value;
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_parent_mass_channel", label,
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.velocity = ?velocity, rank = admitted.rank(),
                file.value = %cooking.uncook(value.expression().as_view()).to_canonical_string(),
                "Admitting each actual mass-projected branch before equality"
            );
            eyre::ensure!(
                admitted.structure().canonical().is_scalar(),
                "{label}: mass projector left unmatched native spin/color ports"
            );
        }
        let outgoing = *boundary
            .iter()
            .find(|h| graph.flow(**h) == Flow::Source)
            .unwrap();
        let row = graph.underlying[outgoing]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom();
        let spinor_ports = |atom: &Atom| {
            let mut matcher = SlotMatcher::default();
            let mut ports = BTreeSet::new();
            atom.visitor(&mut |view| {
                if matcher.parse::<Bispinor, Aind>(view).is_ok() {
                    ports.insert(view.to_owned());
                    false
                } else {
                    true
                }
            });
            ports
        };
        // Only this isolated external block is represented in its finite
        // physical basis. Every indexed Lorentz port must be fixed by a
        // delta coordinate 0..3; no internal source or loop vector enters.
        let physical_external_block = |block: &Atom| -> Result<Atom> {
            eyre::ensure!(
                !block.contains_symbol(GS.emr_mom)
                    && !block.contains_symbol(GS.den)
                    && !block.contains_symbol(SPENSO_TAG.trace),
                "the finite block cannot contain internal momentum, poles or traces"
            );
            let mut matcher = SlotMatcher::default();
            let mut lorentz = BTreeSet::new();
            let mut fixed = BTreeSet::new();
            let mut invalid = false;
            block.visitor(&mut |view| {
                if matcher.parse::<Minkowski, Aind>(view).is_ok() {
                    lorentz.insert(view.to_owned());
                    return false;
                }
                if let AtomView::Fun(function) = view {
                    let head = function.get_symbol();
                    invalid |= head.has_tag(&SPENSO_TAG.tensor)
                        && ![ETS.metric, AGS.gamma, AGS.gamma0, GS.delta_vec].contains(&head);
                    if head == GS.delta_vec && function.get_nargs() == 2 {
                        invalid |=
                            !matches!(matcher.concrete_component(function.get(0)), Some(Ok(0..=3)));
                        fixed.insert(function.get(1).to_owned());
                    }
                    if head == CS.fundamental_rep {
                        invalid |= function.get(0) != Atom::num(3).as_view();
                    }
                }
                true
            });
            eyre::ensure!(
                !invalid && lorentz == fixed,
                "external projector must use only its explicit physical 0..3 basis"
            );
            let without_lorentz_ports = block.replace_map(|view, _, output| {
                if let AtomView::Fun(port) = view
                    && port.get_symbol() == LibraryRep::from(Minkowski {}).symbol()
                {
                    **output = Atom::one();
                }
            });
            eyre::ensure!(
                !without_lorentz_ports.contains_symbol(GS.dim),
                "the physical block cannot carry a contracted formal-d Lorentz factor"
            );
            Ok(block.map_mink_dim(4))
        };
        let finite_projector =
            physical_external_block(&cooking.uncook(projector.expression().as_view()))?;
        let block_index = idenso::bis!(4, Atom::from(Aind::new_dummy()));
        eyre::ensure!(
            !cooking
                .uncook(source.expression().as_view())
                .contains(block_index.as_view())
        );
        let identity = function!(
            spenso::network::library::symbolic::ETS.metric,
            &row,
            &block_index
        );
        let gamma = on_shell_tests::gamma_velocity_at_port(graph, &identity, outgoing.0, &velocity);
        let gamma = physical_external_block(&gamma)?;
        for sign in [1, -1] {
            let p = (&identity + Atom::num(sign) * &gamma) / Atom::num(2);
            let expected = if sign == 1 {
                finite_projector
                    .replace(row.to_pattern())
                    .with(block_index.clone())
            } else {
                Atom::Zero
            };
            // This is only the already witnessed finite external block.
            // Native exact library components recognize transposed chain
            // spellings without assigning four dimensions to the source X.
            use spenso::{
                iterators::IteratableTensor,
                network::{ExecutionResult, Sequential, SmallestDegree},
            };
            let difference = physical_external_block(&(&finite_projector * p - expected))?;
            let mut network = difference.as_view().parse_into_net()?;
            let library = crate::utils::TENSORLIB.read().unwrap();
            network.execute::<Sequential, SmallestDegree, _, _, _>(
                &*library,
                &*crate::utils::FUN_LIB,
            )?;
            let stored_entries = match network.result_tensor(&*library)? {
                ExecutionResult::Zero => 0,
                ExecutionResult::One => {
                    return Err(eyre!("finite external projector difference is one"));
                }
                ExecutionResult::Val(tensor) => {
                    assert_eq!(tensor.size()?, 4 * 4 * 3 * 3);
                    let mut entries = 0;
                    for (_, value) in tensor.iter_flat() {
                        eyre::ensure!(
                            value.is_zero(),
                            "finite ordered external Phi P_{sign} has a nonzero exact component: {value}"
                        );
                        entries += 1;
                    }
                    entries
                }
            };
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_parent_finite_external_components", sign,
                graph = %graph.name, child = %child.subgraph().string_label(),
                component_space = 144, stored_entries,
                file.difference = %difference.to_canonical_string(),
                "All native exact spin and SU(3) components of the external projector difference vanish"
            );
        }

        // Extract the identity coefficient of the ACTUAL projector action.
        // Its new spinor port did not exist in X. Setting just the matrices
        // from the graph output to that port to zero recovers X_sigma/2.
        // Reconstruction below certifies this decomposition before the
        // independent source/routing comparison. No child tensor is replaced.
        let certify_branches = |label: &str,
                                input_native: &SymbolicTensor<IndexedStructure>,
                                branches: [&SymbolicTensor<IndexedStructure>; 2]|
         -> Result<()> {
            let untouched_native = input_native.clone();
            let input_source = cooking.uncook(input_native.expression().as_view());
            let input = &input_source;
            let untouched = input.clone();
            let original_ports = spinor_ports(input);
            for (sign, branch_native) in [1, -1].into_iter().zip(branches) {
                let branch_source = cooking.uncook(branch_native.expression().as_view());
                let branch = &branch_source;
                let new_ports = spinor_ports(branch)
                    .difference(&original_ports)
                    .cloned()
                    .collect::<Vec<_>>();
                let [fresh] = new_ports.as_slice() else {
                    return Err(eyre!(
                        "{label}: actual projector must introduce one fresh spinor interface"
                    ));
                };
                let mut selected = 0;
                let identity_part = branch.replace_map(|view, _, output| {
                    if let AtomView::Fun(matrix) = view {
                        let ordered = matrix.get_nargs() >= 2
                            && matrix.get(0) == row.as_view()
                            && matrix.get(1) == fresh.as_view();
                        // Native gamma0 is symmetric; ordinary gamma endpoint
                        // order remains strict, including the spatial block.
                        let temporal = matrix.get_symbol() == AGS.gamma0
                            && matrix.get_nargs() == 2
                            && (ordered
                                || (matrix.get(1) == row.as_view()
                                    && matrix.get(0) == fresh.as_view()));
                        if temporal
                            || (matrix.get_symbol() == AGS.gamma
                                && matrix.get_nargs() == 3
                                && ordered)
                        {
                            selected += 1;
                            **output = Atom::Zero;
                        }
                    }
                });
                eyre::ensure!(
                    selected > 0,
                    "{label}: no actual ordered external projector block found"
                );
                let evaluated = Atom::num(2) * identity_part;
                eyre::ensure!(
                    !evaluated.contains(fresh.as_view()),
                    "{label}: projector coefficient extraction retained the new spinor port"
                );
                let evaluated_native =
                    input_native.with_checked_expression(cooking.try_cook(evaluated.as_view())?)?;
                let row_native = cooking.try_cook(row.as_view())?;
                let row_native = Slot::<LibraryRep, AbstractIndex>::try_from(row_native.as_view())?;
                let fresh_native = cooking.try_cook(fresh.as_view())?;
                let fresh_native =
                    Slot::<LibraryRep, AbstractIndex>::try_from(fresh_native.as_view())?;
                assert_eq!(row_native.rep(), fresh_native.rep());
                let row_position = evaluated_native
                    .structure()
                    .logical_slots()
                    .iter()
                    .position(|slot| slot.to_atom() == row_native.to_atom())
                    .ok_or_else(|| {
                        eyre!("{label}: certified free row is not an actual logical port")
                    })?;
                let renamed_native = evaluated_native.reindex_interface_ports(
                    &std::collections::HashMap::from([(row_position, fresh_native.aind())]),
                )?;
                let renamed = cooking.uncook(renamed_native.expression().as_view());
                assert_eq!(
                    renamed.replace(fresh.to_pattern()).with(row.clone()),
                    evaluated,
                    "only the certified free endpoint may be renamed"
                );
                let identity =
                    function!(spenso::network::library::symbolic::ETS.metric, &row, fresh);
                let block4 =
                    on_shell_tests::gamma_velocity_at_port(graph, &identity, outgoing.0, &velocity);
                let block4 = physical_external_block(&block4)?;
                // Embed ONLY the witnessed physical external block into the
                // unchanged formal-d interface of the actual source.
                let block_d = block4.map_mink_dim(GS.dim);
                let spinor = [LibraryRep::from(Bispinor {})];
                let attachment = ContractSettings {
                    representations: Some(&spinor),
                    rank_one: false,
                    collect_chains: false,
                    collect_traces: false,
                    expand: false,
                    ..Default::default()
                };
                // Remove the helper's auxiliary identity metric before sewing
                // this small block to the unchanged factorized source. Whole-
                // expression contraction can choose a different source layout.
                let cooked_block = cooking.try_cook(block_d.as_view())?;
                assert_eq!(cooking.uncook(cooked_block.as_view()), block_d);
                let block = SymbolicTensor::infer_indexed(cooked_block)?;
                let normalized_block = block.contract(attachment)?;
                assert_eq!(
                    block.structure().canonical().external_structure(),
                    normalized_block
                        .structure()
                        .canonical()
                        .external_structure()
                );
                let normalized_block_source =
                    cooking.uncook(normalized_block.expression().as_view());
                let projected_native = normalized_block.multiply(renamed_native)?;
                let projected_native = projected_native
                    .with_checked_expression(Atom::num(sign) * projected_native.expression())?;
                let template_native = evaluated_native.try_add(&projected_native)?;
                let template_native = template_native
                    .with_checked_expression(template_native.expression() / Atom::num(2))?;
                assert_eq!(
                    branch_native.structure().logical_slots(),
                    template_native.structure().logical_slots()
                );
                let template = cooking.uncook(template_native.expression().as_view());
                let literal_attachment = branch_native.same_value(&template_native);
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_vertex_actual_projector_template", label, sign, literal_attachment,
                    graph = %graph.name, child = %child.subgraph().string_label(),
                    file.actual_branch = %branch.to_canonical_string(),
                    file.actual_source = %evaluated.to_canonical_string(),
                    file.external_block_source = %block_d.to_canonical_string(),
                    file.external_block = %normalized_block_source.to_canonical_string(),
                    file.template = %template.to_canonical_string(),
                    "Certifying actual ordered projector factor and evaluated source separately"
                );
                // Identical operands already give an exact attachment proof.
                // Still admit the tensor and its interface, without repeating
                // expensive normalization of the complete source. Otherwise
                // contract only spinor metrics and canonicalize internal indices.
                let mut attached = Vec::new();
                for value in [branch_native, &template_native]
                    .into_iter()
                    .take(if literal_attachment { 1 } else { 2 })
                {
                    let source = cooking.uncook(value.expression().as_view());
                    let cooked = cooking.try_cook(source.as_view())?;
                    assert_eq!(cooking.uncook(cooked.as_view()), source);
                    let tensor = value.with_checked_expression(value.expression().clone())?;
                    let logical = tensor.structure().logical_slots();
                    let slots = tensor.structure().canonical().external_structure();
                    let normalized = if literal_attachment {
                        tensor
                    } else {
                        tensor.contract(attachment)?.canonize()?
                    };
                    assert_eq!(
                        slots,
                        normalized.structure().canonical().external_structure()
                    );
                    assert_eq!(logical, normalized.structure().logical_slots());
                    attached.push(normalized);
                }
                if !literal_attachment {
                    assert_eq!(
                        attached[0].structure().canonical().external_structure(),
                        attached[1].structure().canonical().external_structure()
                    );
                    eyre::ensure!(
                        attached[0].expression() == attached[1].expression(),
                        "actual branch differs from the native ordered external-block template"
                    );
                }
                let actual_source = at_point(&evaluated_native, sign)?;
                let expected_source = at_point(input_native, sign)?;
                exact_zero(
                    "actual evaluated X_sigma equals independently routed source",
                    &actual_source.try_sub(&expected_source)?,
                )?;
            }
            assert_eq!(
                *input, untouched,
                "formal-d graph source stays byte-exact during external-block proof"
            );
            assert!(input_native.same_value(&untouched_native));
            Ok(())
        };
        certify_branches("complete parent input", &source, [&positive, &negative])?;

        // Compare physical rational terms in one certified chart. Owner
        // metadata is retained below as a separate witness, not mistaken
        // for a scalar denominator. Only the existing denominator owner
        // expands its complete small polynomial; numerators retain d and
        // their factorized source circuits.
        let common_chart_zero = |label: &str,
                                 left: &SymbolicTensor<IndexedStructure>,
                                 right: &SymbolicTensor<IndexedStructure>|
         -> Result<()> {
            let mut coefficients = BTreeMap::<Vec<Atom>, SymbolicTensor<IndexedStructure>>::new();
            let mut witnesses = Vec::new();
            for (side, sign, expression) in [("left", 1, left), ("right", -1, right)] {
                let raw_terms = FourDTerm::from_tensor(expression)?;
                for (term, value) in raw_terms.iter().enumerate() {
                    for denominator in &value.denominators {
                        eyre::ensure!(
                            denominator.source_edge.0 < graph.underlying.n_edges()
                                && parent
                                    .subgraph()
                                    .includes(&graph[&denominator.source_edge].1),
                            "{label}: denominator lost its native parent source owner"
                        );
                        witnesses.push(serde_json::json!({
                            "side":side,"term":term,"source_edge":denominator.source_edge.0,
                            "momentum":denominator.momentum.to_canonical_string(),
                            "mass_squared":denominator.mass_squared.to_canonical_string(),
                            "full_polynomial":denominator.full_expr.to_canonical_string(),
                        }));
                    }
                }
                // Admit compact source contractions before routing their
                // vector arguments. Then delta_i(Q) is an explicit dot,
                // so Q -> Q4-Q2 stays a linear vector contraction rather
                // than putting that sum into rank-one scalar metadata.
                let admitted = expression;
                let dotted = admitted.to_dots()?;
                assert_eq!(admitted.structure(), dotted.structure());
                let physical_expression = GS.resolve_os_momenta(
                    &GS.erase_uv_momentum_provenance(dotted.expression())
                        .replace_multiple(graph.uv_wrapped_replacement(
                            parent.subgraph(),
                            &lmb,
                            &[W_.x___],
                        )),
                );
                let physical = dotted.with_checked_expression(physical_expression)?;
                for term in FourDTerm::from_tensor(&physical)? {
                    let mut polynomials = term
                        .denominators
                        .iter()
                        .map(|denominator| denominator.polynomial_in_lmb(&lmb))
                        .collect::<Result<Vec<_>>>()?;
                    polynomials.sort(); // repeated entries retain propagator powers
                    let numerator = term
                        .numerator
                        .with_checked_expression(Atom::num(sign) * term.numerator.expression())?;
                    if let Some(sum) = coefficients.get_mut(&polynomials) {
                        *sum = sum.try_add(&numerator)?;
                    } else {
                        coefficients.insert(polynomials, numerator);
                    }
                }
            }
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_parent_common_chart", label,
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.left = %cooking.uncook(left.expression().as_view()).to_canonical_string(), file.right = %cooking.uncook(right.expression().as_view()).to_canonical_string(),
                file.owners = %serde_json::to_string(&witnesses)?,
                file.lmb = ?lmb, denominator_classes = coefficients.len(),
                "Preserved native owners before formal-d denominator-class comparison"
            );
            for (class, (denominators, numerator)) in coefficients.into_iter().enumerate() {
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_vertex_parent_denominator_class", label, class,
                    file.denominators = ?denominators,
                    file.numerator = %cooking.uncook(numerator.expression().as_view()).to_canonical_string(),
                    "Only identical full denominator polynomials share a numerator certificate"
                );
                frame_zero(&format!("{label}, denominator class {class}"), &numerator)?;
            }
            Ok(())
        };
        let ordinary = t_raw(&source, ctx, parent, child, &lmb)?;
        if child.renormalization_scheme() != ApproximationType::OS {
            eyre::ensure!(
                !ordinary.expression().contains_symbol(AGS.gamma0)
                    && !ordinary.expression().contains_symbol(ETS.delta),
                "ordinary vacuum Taylor grading cannot introduce a fixed frame"
            );
        }
        let on_shell = positive.try_add(&negative)?;
        let on_shell_resolved =
            on_shell.with_checked_expression(GS.resolve_os_momenta(on_shell.expression()))?;
        let overlap = t_raw(&on_shell_resolved, ctx, parent, child, &lmb)?;
        let [ordinary_positive, ordinary_negative] =
            graph.on_shell_branches(&graph.model, &ordinary, parent.subgraph(), &lmb, &velocity)?;
        certify_branches(
            "ordinary parent coefficient",
            &ordinary,
            [&ordinary_positive, &ordinary_negative],
        )?;
        common_chart_zero(
            "actual open-tensor U1 O equals O U1",
            &overlap,
            &ordinary_positive.try_add(&ordinary_negative)?,
        )?;
        let actual = uv_limit(prefix, ctx, parent, child, parent, child)?;
        common_chart_zero(
            "actual parent U1/O/U1O composition",
            actual.atom(),
            &overlap.try_sub(&ordinary)?.try_sub(&on_shell)?,
        )?;
        let mut vacuum_sectors = actual.active_sectors().iter().map(|s| s.integration_atom());
        let mut vacuum = vacuum_sectors
            .next()
            .ok_or_else(|| eyre!("actual parent has no analytic sector"))?
            .clone();
        for sector in vacuum_sectors {
            vacuum = vacuum.try_add(sector)?;
        }
        common_chart_zero(
            "actual analytic input is only the signed vacuum remainder",
            &vacuum,
            &overlap.try_sub(&ordinary)?,
        )?;
        // The certified finite Phi P block, actual branch templates,
        // source gates, and U1/O commutation imply Phi[-(U1+O-U1O)X]
        // = -Phi[X], without evaluating mixed-temporal formal-d traces.
        assert!(actual.recursive_completion().is_empty());
        assert_eq!(actual.active_sectors().len(), prefix.sectors().count());
        debug_tags!(#uv, #integrated, #inspect;
            stage = "native_mass_parent_input_certificate", source = %diagram.name(),
            given = %child.subgraph().string_label(), given_scheme = %child.renormalization_scheme(),
            file.velocity = ?velocity, file.parent_lmb = ?lmb,
            actual_parent_compared = true, actual_analytic_input_compared = true,
            integrated_parent_finite_computed = false, public_dispatch_admitted = false,
            "Actual parent/projector/source identities; actual finite integration is a separate next gate"
        );
        Ok(actual)
    }

    fn certify_parent_mass_projection(&self, subtraction_velocity: [&str; 4]) -> Result<()> {
        let reference = self;
        let mut certified = reference.certify_integrated_qg_children(&[6, 6, 21])?;
        certified.extend(reference.certify_integrated_qg_children(&[6, 21, 21])?);
        assert_eq!(certified.len(), 4);
        let exact_zero = Self::certify_parent_tensor_zero;
        for record in certified {
            let graph = &record.graph;
            let child = &record.child;
            let mut settings = record.settings.clone();
            let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
            let mut parents = graph
                .classified_spinneys(&uncut, &graph.model, &settings, &graph.loop_momentum_basis)?
                .into_iter()
                .filter(|spinney| {
                    spinney.n_components() == 1
                        && spinney.dod == 1
                        && spinney.renormalization_scheme == ApproximationType::OS
                        && graph.n_loops(spinney.filter()) == 2
                        && spinney.filter().includes(child.subgraph())
                })
                .collect::<Vec<_>>();
            assert_eq!(parents.len(), 1);
            let parent = OwnedForestNode {
                spinney: parents.remove(0),
                topo_order: 2,
            };
            let guard_error = validate_local_counterterm_schemes(
                graph,
                &super::super::CutStructure::empty(graph),
                &settings,
            )
            .expect_err("this private proof must not admit the public OS/MUV forest");
            let guard_error = format!("{guard_error:#}");
            assert!(
                guard_error.contains("on-shell/OS component")
                    && guard_error.contains("unsupported nested or disconnected scheme policy")
            );
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_prefix_parent_proof_start",
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.child_edges = ?graph.paired_edges(child.subgraph()),
                file.child_lmb = ?child.lmb(), public_guard_error = %guard_error,
                "Inspecting actual certified child before private parent composition"
            );

            let prefix = Full4dCts::recursion_input(
                &record.local,
                &record.integrated,
                ApproximationType::MUV,
                false,
                child.lmb(),
            )?;
            let sectors = prefix.sectors().collect::<Vec<_>>();
            assert_eq!(
                sectors.len(),
                2,
                "one local child plus one actual finite completion"
            );
            let active = sectors
                .iter()
                .filter(|s| s.frozen_lmbs().is_empty())
                .collect::<Vec<_>>();
            let frozen = sectors
                .iter()
                .filter(|s| !s.frozen_lmbs().is_empty())
                .collect::<Vec<_>>();
            assert_eq!(active.len(), 1);
            assert_eq!(frozen.len(), 1);
            assert_eq!(active[0].active_components.len(), 1);
            assert_eq!(&active[0].active_components[0].0, child.subgraph());
            assert_eq!(&active[0].active_components[0].1, child.subgraph());
            assert_eq!(active[0].active_components[0].2.loop_edges.len(), 1);
            assert!(frozen[0].active_components.is_empty());
            assert_eq!(frozen[0].frozen_lmbs(), std::slice::from_ref(child.lmb()));
            let physical_finite = record.integrated.physical_finite_counterterm()?;
            let bookkeeping_finite = physical_finite.with_checked_expression(
                physical_finite.expression() * Atom::var(GS.integrated_loop_scale).pow(4),
            )?;
            exact_zero(
                "consumed one-loop measure marker",
                &frozen[0].atom.try_sub(&bookkeeping_finite)?,
            )?;
            let child_integral = record
                .integrated
                .physical_pole()?
                .try_sub(&physical_finite)?;
            let frozen_physical = frozen[0].atom.with_checked_expression(
                frozen[0]
                    .atom
                    .expression()
                    .replace(GS.integrated_loop_scale)
                    .with(1),
            )?;
            let completed_integral = child_integral.try_add(&frozen_physical)?;
            let completed_integral = completed_integral.with_checked_expression(
                graph
                    .model
                    .expand_couplings(completed_integral.expression())
                    .replace(GS.m_uv_vacuum)
                    .with(1)
                    .replace(GS.m_uv_expansion)
                    .with(1)
                    .replace(GS.mu_r_sq)
                    .with(1),
            )?;
            let pole_difference = completed_integral.try_sub(&record.pole_reference)?;
            let pole_difference = pole_difference.with_checked_expression(
                pole_difference
                    .expression()
                    .series(GS.dim_epsilon, 0, 1)?
                    .to_atom(),
            )?;
            exact_zero(
                "completed native vertex equals independent pole tensor",
                &pole_difference,
            )?;

            // Use actual component/cograph ownership; no synthetic replacement of
            // the completed local child enters the parent operation.
            let reduced = parent.reduced_subgraph(child);
            let cograph_edges = graph.paired_edges(&reduced);
            assert_eq!(cograph_edges.len(), 2);
            let mut cograph_pdgs = cograph_edges
                .iter()
                .map(|edge| {
                    graph
                        .model
                        .particle_by_id(graph[*edge].particle.particle().unwrap())
                        .map(|particle| particle.pdg_code.abs())
                })
                .collect::<Result<Vec<_>, _>>()?;
            cograph_pdgs.sort_unstable();
            assert_eq!(cograph_pdgs, [6, 21]);
            let child_frame = &active[0].active_components[0].2;
            for edge in &cograph_edges {
                assert!(
                    child_frame.edge_signatures[*edge]
                        .internal
                        .iter()
                        .all(|sign| sign.is_zero()),
                    "cograph edge must be independent of the actual child hard integration coordinate"
                );
            }
            let cograph_numerator = graph
                .numerator(&reduced, child.subgraph())?
                .to_d_dim(GS.dim)?
                .state
                .expr;
            let cograph_denominator = graph.denominator(&reduced, &graph.model, |_| 1);
            let denominator_terms = FourDTerm::from_view(cograph_denominator.pow(-1).as_view())?;
            assert_eq!(denominator_terms.len(), 1);
            assert_eq!(denominator_terms[0].denominators.len(), 2);
            assert_eq!(
                denominator_terms[0]
                    .denominators
                    .iter()
                    .map(|d| d.source_edge)
                    .collect::<BTreeSet<_>>(),
                cograph_edges.iter().copied().collect()
            );

            // The child boundary gluon becomes an INTERNAL summed parent slot.
            // Check its native formal representation before the parent finalizer.
            let finite_tensor = &frozen[0].atom;
            let mut slots = SlotMatcher::default();
            let mut child_lorentz_ports = 0;
            for slot in finite_tensor.structure().logical_slots() {
                let atom = slot.to_atom();
                if slots
                    .parse::<Minkowski, AbstractIndex>(atom.as_view())
                    .is_ok()
                {
                    let AtomView::Fun(port) = atom.as_view() else {
                        unreachable!()
                    };
                    assert_eq!(port.get(0), Atom::var(GS.dim).as_view());
                    child_lorentz_ports += 1;
                }
            }
            assert_eq!(child_lorentz_ports, 1);
            let attached = prefix.atom().clone().multiply(cograph_numerator)?;
            let attached =
                attached.with_checked_expression(attached.expression() / &cograph_denominator)?;
            for slot in attached.structure().logical_slots() {
                let atom = slot.to_atom();
                assert!(
                    slots
                        .parse::<Minkowski, AbstractIndex>(atom.as_view())
                        .is_err(),
                    "the actual child gluon port must close against its cograph partner"
                );
            }

            settings.on_shell_subtraction_velocity = subtraction_velocity.map(str::to_owned);
            let ctx = UVCtx::new(graph, &graph.model, &settings);
            let lmb = compatible_lmb(&ctx, &parent, child, &prefix)?;
            // Audit the actual hard grading before extracting the parent's U1.
            // No frozen localization kernel is part of this coefficient.
            let t = Atom::var(GS.rescale);
            let scaled = frozen[0].atom.with_checked_expression(graph.uv_rescaled(
                parent.subgraph(),
                2,
                &lmb,
                &lmb,
                frozen[0].atom.expression(),
            ))?;
            let expected_scaled = frozen[0].atom.with_checked_expression(
                frozen[0]
                    .atom
                    .expression()
                    .replace(GS.m_uv_vacuum)
                    .with(Atom::var(GS.m_uv_vacuum) / &t)
                    .replace(GS.integrated_loop_scale)
                    .with(Atom::var(GS.integrated_loop_scale) * &t)
                    * t.pow(-8),
            )?;
            exact_zero(
                "completed-child mass and consumed-measure hard grading",
                &scaled.try_sub(&expected_scaled)?,
            )?;
            let actual = self.certify_parent_mass_projection_input(
                &record.diagram,
                &ctx,
                &parent,
                child,
                &prefix,
            )?;
            let velocity = settings.on_shell_velocity()?;
            assert_eq!(actual.active_sectors().len(), 2);
            assert!(actual.recursive_completion().is_empty());
            assert_eq!(
                actual
                    .active_sectors()
                    .iter()
                    .filter(|s| !s.frozen_lmbs().is_empty())
                    .count(),
                1
            );
            let output_frozen = actual
                .active_sectors()
                .iter()
                .find(|sector| !sector.frozen_lmbs().is_empty())
                .expect("the actual parent retains its one frozen completion");
            assert_eq!(
                output_frozen.frozen_lmbs(),
                std::slice::from_ref(child.lmb()),
                "the actual parent output retains the complete child localization basis"
            );
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_prefix_parent_mass_projection",
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.velocity = ?velocity, file.parent_lmb = ?lmb,
                active_sectors = actual.active_sectors().len(), frozen_sectors = 1,
                public_dispatch_admitted = false,
                "Actual native child prefix satisfies matched parent mass-projection identities"
            );
        }
        Ok(())
    }
}

impl QcdMassReference {
    /// Actual proper-local W only. Two calibrated scalar traces, not open R2.
    pub(in crate::uv) fn certify_d09_proper_local_vacuum_angular_traces(
        &self,
        diagram: &feynkit_graph::FeynmanDiagram,
        ctx: &UVCtx<'_>,
        given: &OwnedForestNode,
        prefix: &Full4dCts,
        actual: &(OwnedForestNode, Local4dCts, IntegratedCts),
        backend: &vakint::VakintSettings,
    ) -> Result<[SymbolicTensor<IndexedStructure>; 4]> {
        use crate::uv::Spinney;
        let (parent, local, _) = actual;
        let graph = ctx.graph;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        eyre::ensure!(diagram.name() == "D09" && graph.name == "D09");
        eyre::ensure!(parent.renormalization_scheme() == ApproximationType::OS);
        eyre::ensure!(given.renormalization_scheme() == ApproximationType::OS);
        eyre::ensure!(parent.dod() == 1 && given.dod() == 1);
        eyre::ensure!(graph.n_loops(parent.subgraph()) == 2);
        eyre::ensure!(graph.n_loops(given.subgraph()) == 1);
        eyre::ensure!(graph.paired_edges(given.subgraph()) == [EdgeIndex(2), EdgeIndex(3)]);
        eyre::ensure!(ctx.settings.generate_integrated && !ctx.settings.add_marker);
        eyre::ensure!(ctx.settings.vakint.normalization == "MSbar");
        eyre::ensure!(ctx.settings.vakint.additional_normalization == "1");
        eyre::ensure!(matches!(
            &backend.integral_normalization_factor,
            vakint::LoopNormalizationFactor::MSbar
        ));
        let velocity = ctx.settings.on_shell_velocity()?;
        eyre::ensure!(velocity == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]);
        eyre::ensure!(local.recursive_completion().is_empty());
        let selected = local
            .active_sectors()
            .iter()
            .filter(|sector| sector.frozen_lmbs.is_empty() && sector.active_components.len() == 2)
            .collect::<Vec<_>>();
        let [sector] = selected.as_slice() else {
            return Err(eyre!(
                "D09 needs exactly one original two-active proper-local W sector"
            ));
        };
        let sector = *sector;
        let retained_local = local.clone();
        let retained_w = sector.integration_atom().clone();
        eyre::ensure!(
            !retained_w.expression().contains_symbol(GS.dim_epsilon)
                && !retained_w
                    .expression()
                    .contains_symbol(GS.integrated_loop_scale),
            "the proper-local UV gate excludes integrated child coefficients"
        );
        eyre::ensure!(
            sector.active_components[0]
                .2
                .loop_edges
                .iter()
                .copied()
                .eq([EdgeIndex(2)])
        );
        eyre::ensure!(
            sector.active_components[1]
                .2
                .loop_edges
                .iter()
                .copied()
                .eq([EdgeIndex(4)])
        );
        eyre::ensure!(sector.inherited_loop_edges == [EdgeIndex(2), EdgeIndex(4)]);
        eyre::ensure!(sector.scheme_integrand.is_some() && !retained_w.same_value(&sector.atom));
        let actual_open = &sector.atom;
        let (closing_duals, open_basis): (
            [SymbolicTensor<IndexedStructure>; 2],
            [SymbolicTensor<IndexedStructure>; 2],
        ) = {
            use crate::numerator::AppliedFeynmanRule;
            use color_eyre::eyre::{ensure, eyre};
            use idenso::color::CS;
            use idenso::dirac::AGS;
            use idenso::representations::Bispinor;
            use idenso::{CookMode, CookSettings};
            use spenso::network::library::symbolic::ETS;
            use spenso::structure::representation::LibraryRep;
            use spenso::structure::slot::IsAbstractSlot;
            use std::collections::HashMap;
            use symbolica::atom::{Atom, AtomCore, AtomView};

            ensure!(diagram.name() == "D09");
            ensure!(actual_open.rank() == 4);
            let cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            let (projector, reference_p) = self.mass_projector(diagram, &Atom::var(GS.dim))?;
            let projector = projector
                .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
            let projector = projector
                .with_checked_expression(projector.expression().replace(CS.nc).with(3))?
                .try_into_indexed()?;
            ensure!(projector.rank() == 4);
            let AtomView::Fun(reference_p) = reference_p.as_view() else {
                return Err(eyre!("D09 dual has no native incoming momentum owner"));
            };
            ensure!(reference_p.get_symbol() == GS.emr_mom && reference_p.get_nargs() == 1);
            let incoming_edge =
                usize::try_from(reference_p.get(0)).map_err(|error| eyre!(error))?;

            // Only the actual projector's one incoming momentum is removed/evaluated.
            // Its original bispinor/color ports and all representations stay attached.
            let mut momentum_count = 0usize;
            let mut bad_momentum = false;
            projector.expression().visitor(&mut |view| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    momentum_count += 1;
                    bad_momentum |= momentum.get_nargs() != 2
                        || usize::try_from(momentum.get(0)).ok() != Some(incoming_edge);
                }
                true
            });
            ensure!(momentum_count > 0 && !bad_momentum);
            let zero_projector = projector.with_checked_expression(
                projector.expression().replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                    {
                        **output = Atom::Zero;
                    }
                }),
            )?;
            let closing_identity = zero_projector
                .with_checked_expression(Atom::num(12) * zero_projector.expression())?;
            ensure!(!closing_identity.expression().contains_symbol(GS.emr_mom));

            // The independent temporal word uses the SAME literal source spin metric
            // endpoints. Replace only that one bis4 metric; its native color closure
            // is untouched. No anonymous trace or new external port is introduced.
            let bis = LibraryRep::from(Bispinor {}).symbol();
            let identity_source = cooking.uncook(closing_identity.expression().as_view());
            let mut replaced_spin_metric = 0usize;
            let temporal_source = identity_source.replace_map(|view, _, output| {
                let AtomView::Fun(metric) = view else { return };
                if metric.get_symbol() != ETS.metric || metric.get_nargs() != 2 {
                    return;
                }
                let spin_metric = (0..2).all(|index| {
                    matches!(metric.get(index), AtomView::Fun(rep)
                if rep.get_symbol() == bis && rep.get_nargs() == 2
                    && rep.get(0) == Atom::num(4).as_view())
                });
                if spin_metric {
                    replaced_spin_metric += 1;
                    **output = function!(AGS.gamma0, metric.get(0), metric.get(1));
                }
            });
            ensure!(replaced_spin_metric == 1);
            let admitted_temporal = AppliedFeynmanRule::new(temporal_source)?.expr;
            ensure!(admitted_temporal.rank() == 4);
            ensure!(
                admitted_temporal
                    .structure()
                    .canonical()
                    .external_structure()
                    == closing_identity
                        .structure()
                        .canonical()
                        .external_structure()
            );
            let closing_temporal =
                closing_identity.with_checked_expression(admitted_temporal.expression().clone())?;

            // Independently bind gamma0 to the p-linear source projector at the unit
            // rest point. energy_delta embeds that fixed component into the ORIGINAL
            // formal-D Lorentz representation; no scalar D->4 replacement occurs.
            let closing_p = projector.try_sub(&zero_projector)?;
            let closing_p =
                closing_p.with_checked_expression(Atom::num(12) * closing_p.expression())?;
            let rest_p = closing_p.with_checked_expression(closing_p.expression().replace_map(
                |view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                    {
                        **output = GS.energy_delta(momentum.get(1));
                    }
                },
            ))?;
            let rest_difference = rest_p.try_sub(&closing_temporal)?;
            let rest_difference = Self::reduce_parent_external_frame(
                "D09 independent temporal dual equals its native rest momentum source",
                &rest_difference,
            )?;
            Self::certify_parent_tensor_zero(
                "D09 independent temporal dual equals its native rest momentum source",
                &rest_difference,
            )?;

            // Reverse each actual endpoint pair, retaining the literal inverse and
            // zero interface, exactly as in the native proper-prefix dual owner.
            let transpose = |word: &SymbolicTensor<IndexedStructure>|
     -> color_eyre::Result<SymbolicTensor<IndexedStructure>> {
        let slots = word.structure().logical_slots();
        ensure!(slots.len() == 4);
        let mut map = HashMap::new();
        let mut inverse = HashMap::new();
        for (position, slot) in slots.iter().enumerate() {
            let opposite = slots.iter().filter(|other|
                slot.rep().matches(&other.rep()) && other.aind() != slot.aind()
            ).collect::<Vec<_>>();
            ensure!(opposite.len() == 1,
                "D09 dual requires one actual opposite spin/color port");
            map.insert(position, opposite[0].aind());
            inverse.insert(position, slot.aind());
        }
        let open = word.reindex_interface_ports(&map)?;
        ensure!(open.reindex_interface_ports(&inverse)?.same_value(word));
        ensure!(open.structure().canonical().external_structure()
            == actual_open.structure().canonical().external_structure());
        let zero = word.with_checked_expression(Atom::Zero)?;
        let open_zero = zero.reindex_interface_ports(&map)?;
        ensure!(open_zero.rank() == 4);
        ensure!(open_zero.reindex_interface_ports(&inverse)?.same_value(&zero));
        Ok(open)
    };
            let open_basis = [transpose(&closing_identity)?, transpose(&closing_temporal)?];
            let closing_duals = [
                closing_identity
                    .with_checked_expression(closing_identity.expression() / Atom::num(12))?,
                closing_temporal
                    .with_checked_expression(closing_temporal.expression() / Atom::num(12))?,
            ];

            // Full finite 2x2 Gram calibration: Tr_spin,color(I)=4*3,
            // gamma0^2=I, and both cross traces vanish. Use the existing native
            // component oracle on this closed bis4/cof3 product only. Symmetric
            // gamma0 endpoints can compress its self-pair to a squared chain;
            // this calibration does not require symbolic power contraction or
            // replace any active-D momentum, denominator or vacuum numerator.
            for (row, dual) in closing_duals.iter().enumerate() {
                for (column, basis) in open_basis.iter().enumerate() {
                    let trace = dual.clone().multiply(basis.clone())?;
                    ensure!(
                        trace.is_scalar(),
                        "D09 calibrated dual left native ports open"
                    );
                    let expected = Atom::num(if row == column { 1 } else { 0 });
                    let trace = cooking.uncook(trace.expression().as_view());
                    let component = on_shell_tests::scalar_component(&trace)?;
                    ensure!(
                        (component - expected).collect_compact_factors().is_zero(),
                        "D09 unit-rest identity/temporal dual pairing ({row}, {column})"
                    );
                }
            }
            (closing_duals, open_basis)
        };

        // The point belongs to the actual incoming source and its signed
        // opposite boundary, not to a new momentum or an external projector.
        let (_, reference_p) = self.mass_projector(diagram, &Atom::var(GS.dim))?;
        let basis = diagram.loop_momentum_basis();
        let routed_p = basis.route_expression(&reference_p);
        let boundary = graph.dummy_less_full_crown(parent.subgraph());
        let mut signs = BTreeMap::new();
        let mut incoming = None;
        for hedge in boundary.included_iter() {
            let edge = graph[&hedge];
            let routed = basis.route_expression(&function!(GS.emr_mom, edge.0));
            let sign = if (&routed - &routed_p).expand().is_zero() {
                1
            } else if (&routed + &routed_p).expand().is_zero() {
                -1
            } else {
                return Err(eyre!("D09 UV dual has an unmatched external routing"));
            };
            eyre::ensure!(signs.insert(edge, sign).is_none());
            if graph.flow(hedge) == Flow::Sink {
                eyre::ensure!(sign == 1 && incoming.replace(edge).is_none());
            }
        }
        eyre::ensure!(signs.len() == 2);
        let incoming =
            incoming.ok_or_else(|| eyre!("D09 UV input has no physical incoming owner"))?;
        let mass = graph[incoming].particle.mass_atom(ctx.model);
        for (edge, pdg) in [(2, 6), (3, 21), (4, 6), (5, 6), (6, 21)] {
            let data = &graph[EdgeIndex(edge)];
            eyre::ensure!(
                ctx.model
                    .particle_by_id(data.particle.particle().unwrap())?
                    .pdg_code
                    .abs()
                    == pdg
            );
            if pdg == 6 {
                eyre::ensure!(data.particle.mass_atom(ctx.model) == mass);
            } else {
                eyre::ensure!(data.particle.mass_atom(ctx.model).is_zero());
            }
        }
        let local_prefixes = prefix
            .sectors()
            .filter(|value| value.frozen_lmbs.is_empty() && value.active_components.len() == 1)
            .collect::<Vec<_>>();
        let [local_prefix] = local_prefixes.as_slice() else {
            return Err(eyre!(
                "D09 requires the actual one-active local child prefix"
            ));
        };
        let local_prefix = *local_prefix;
        eyre::ensure!(
            local_prefix.active_components[0]
                .2
                .loop_edges
                .iter()
                .copied()
                .eq([EdgeIndex(2)])
        );
        let isolated = Full4dCts(FourDSectors::new(vec![local_prefix.clone()], Vec::new())?);
        let lmb = compatible_lmb(ctx, parent, given, &isolated)?;
        eyre::ensure!(
            lmb.loop_edges
                .iter()
                .copied()
                .eq([EdgeIndex(2), EdgeIndex(4)])
        );

        // All point/representation changes happen on a numerator AFTER the
        // actual operators. Negative denominator polynomials stay dimensional.
        let prepare_numerator = |numerator: &SymbolicTensor<IndexedStructure>,
                                 physical_four: bool|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let mut invalid_denominator = false;
            numerator.expression().visitor(&mut |view| {
                if let AtomView::Pow(power) = view
                    && power.get_base().contains_symbol(GS.den)
                    && i64::try_from(power.get_exp()).is_ok_and(|power| power < 0)
                {
                    invalid_denominator = true;
                }
                true
            });
            let expression = numerator.expression().replace_map(|view, _, output| {
                if let AtomView::Fun(denominator) = view
                    && denominator.get_symbol() == GS.den
                {
                    if denominator.get_nargs() != 4 {
                        invalid_denominator = true;
                    } else {
                        // Positive denominator factors belong to the NUMERATOR.
                        // Their fourth argument is the complete scalar value.
                        **output = denominator.get(3).to_owned();
                    }
                }
            });
            eyre::ensure!(
                !invalid_denominator,
                "D09 numerator retained a negative or malformed denominator shell"
            );
            let dotted = numerator.with_checked_expression(expression)?.to_dots()?;
            let expression = GS
                .erase_uv_momentum_provenance(dotted.expression())
                .replace_multiple(graph.uv_wrapped_replacement(
                    parent.subgraph(),
                    &lmb,
                    &[W_.x___],
                ));
            let mut unexpected = false;
            let expression = expression.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    let Ok(edge) = usize::try_from(momentum.get(0)) else {
                        unexpected = true;
                        return;
                    };
                    let edge = EdgeIndex(edge);
                    if let Some(sign) = signs.get(&edge) {
                        let point =
                            GS.on_shell_momentum(incoming, &(Atom::num(*sign) * &mass), &velocity);
                        let indices = momentum
                            .iter()
                            .skip(1)
                            .map(|index| index.to_owned())
                            .collect::<Vec<_>>();
                        **output = GS.indexed_momentum(&point, &indices);
                    } else if ![EdgeIndex(2), EdgeIndex(4)].contains(&edge) {
                        unexpected = true;
                    }
                }
            });
            eyre::ensure!(!unexpected, "D09 UV numerator has an unowned momentum");
            let expression = GS
                .resolve_os_momenta(&expression)
                .replace(mass.to_pattern())
                .with(1);
            let routed = dotted.with_checked_expression(expression)?;
            if !physical_four {
                return Ok(routed);
            }
            let source = cooking.uncook(routed.expression().as_view());
            let outside_representations = source.replace_map(|view, _, output| {
                if view.get_symbol() == Some(LibraryRep::from(Minkowski {}).symbol()) {
                    **output = Atom::one();
                }
            });
            eyre::ensure!(
                !outside_representations.contains_symbol(GS.dim),
                "D09 physical-four lift refuses an already contracted active scalar D"
            );
            let mut invalid_dimension = false;
            let mapped = routed.map_representation_dimensions(|rep, dimension| {
                if rep == LibraryRep::from(Minkowski {}) {
                    if dimension != Atom::var(GS.dim).as_view()
                        && dimension != Atom::num(4).as_view()
                    {
                        invalid_dimension = true;
                    }
                    Atom::num(4)
                } else {
                    dimension.to_owned()
                }
            })?;
            eyre::ensure!(!invalid_dimension);
            eyre::ensure!(mapped.structure().logical_slots() == routed.structure().logical_slots());
            Ok(mapped)
        };
        let rational_zero = |label: &str,
                             left: &SymbolicTensor<IndexedStructure>,
                             right: &SymbolicTensor<IndexedStructure>,
                             physical_four: bool|
         -> Result<()> {
            let mut classes = BTreeMap::<Vec<Atom>, SymbolicTensor<IndexedStructure>>::new();
            for (sign, expression) in [(1, left), (-1, right)] {
                for term in FourDTerm::from_tensor(expression)? {
                    let mut key = term
                        .denominators
                        .iter()
                        .map(|denominator| {
                            denominator.polynomial_in_lmb(&lmb).map(|polynomial| {
                                polynomial.replace(GS.m_uv_expansion).with(GS.m_uv_vacuum)
                            })
                        })
                        .collect::<Result<Vec<_>>>()?;
                    key.sort();
                    let numerator = prepare_numerator(&term.numerator, physical_four)?;
                    let numerator = numerator
                        .with_checked_expression(Atom::num(sign) * numerator.expression())?;
                    if let Some(sum) = classes.get_mut(&key) {
                        *sum = sum.try_add(&numerator)?;
                    } else {
                        classes.insert(key, numerator);
                    }
                }
            }
            for numerator in classes.into_values() {
                let difference = Self::reduce_parent_external_frame(label, &numerator)?;
                Self::certify_parent_tensor_zero(label, &difference)?;
            }
            Ok(())
        };

        // Exact source absorption. Keep the parent's hard rescaling of mUV
        // and its terminal denominator policy inside the existing t_raw owner.
        let root = OwnedForestNode {
            spinney: Spinney::empty(graph),
            topo_order: 0,
        };
        let bare = Full4dCts::from_factorized_local(&Local4dCts::root());
        let child_lmb = compatible_lmb(ctx, given, &root, &bare)?;
        eyre::ensure!(child_lmb.loop_edges.iter().copied().eq([EdgeIndex(2)]));
        let child_source = grow(&bare, ctx, given, &root)?;
        let signed_child_u = -t_raw(&child_source, ctx, given, &root, &child_lmb)?;
        eyre::ensure!(
            signed_child_u.structure().canonical().external_structure()
                == local_prefix
                    .atom
                    .structure()
                    .canonical()
                    .external_structure()
        );
        let co_graph = grow_factor(ctx, parent, given)?;
        let absorbed_source = finalize(signed_child_u.multiply(co_graph)?)?;
        let actual_source = grow(&isolated, ctx, parent, given)?;
        let ordinary = t_raw(&actual_source, ctx, parent, given, &lmb)?;
        let absorbed = t_raw(&absorbed_source, ctx, parent, given, &lmb)?;
        for physical_four in [false, true] {
            rational_zero(
                "actual D09 parent-U absorbs child O-UO at unit rest",
                &ordinary,
                &absorbed,
                physical_four,
            )?;
        }
        let [positive, negative] = graph.on_shell_branches(
            ctx.model,
            &actual_source,
            parent.subgraph(),
            &lmb,
            &velocity,
        )?;
        let on_shell = positive.try_add(&negative)?;
        let overlap = t_raw(
            &on_shell.with_checked_expression(GS.resolve_os_momenta(on_shell.expression()))?,
            ctx,
            parent,
            given,
            &lmb,
        )?;
        let [positive, negative] =
            graph.on_shell_branches(ctx.model, &ordinary, parent.subgraph(), &lmb, &velocity)?;
        let [absorbed_positive, absorbed_negative] = graph.on_shell_branches(
            ctx.model,
            &absorbed_source,
            parent.subgraph(),
            &lmb,
            &velocity,
        )?;
        let absorbed_on_shell = absorbed_positive.try_add(&absorbed_negative)?;
        let absorbed_overlap = t_raw(
            &absorbed_on_shell
                .with_checked_expression(GS.resolve_os_momenta(absorbed_on_shell.expression()))?,
            ctx,
            parent,
            given,
            &lmb,
        )?;
        for physical_four in [false, true] {
            rational_zero(
                "actual D09 proper-local parent UO=OU",
                &overlap,
                &positive.try_add(&negative)?,
                physical_four,
            )?;
            rational_zero(
                "actual D09 original proper-local W is UO-U",
                &retained_w,
                &overlap.try_sub(&ordinary)?,
                physical_four,
            )?;
            rational_zero(
                "actual D09 original proper-local W equals absorbed-child W",
                &retained_w,
                &absorbed_overlap.try_sub(&absorbed)?,
                physical_four,
            )?;
        }

        // Independently calibrate the source phase/color/coupling before any
        // new integral is read. D09's three quark and two gluon propagators
        // and four actual GC_11=iG vertices give i*G^4*CF^2 at zero momenta.
        let coupling = Atom::var(symbol!("UFO::G"));
        eyre::ensure!(
            graph
                .model
                .expand_couplings(&Atom::var(symbol!("UFO::GC_11")))
                == Atom::i() * &coupling
        );
        let empty = graph.empty_subgraph();
        let source = graph
            .numerator(parent.subgraph(), &empty)?
            .to_d_dim(GS.dim)?
            .color_simplify()?
            .state
            .expr;
        let zero_momenta = source
            .expression()
            .replace_map(|view, _, output| {
                if matches!(view, AtomView::Fun(momentum) if momentum.get_symbol() == GS.emr_mom) {
                    **output = Atom::Zero;
                }
            })
            .replace(mass.to_pattern())
            .with(1);
        let source = AppliedFeynmanRule::new(cooking.uncook(zero_momenta.as_view()))?.expr;
        eyre::ensure!(
            source.structure().canonical().external_structure()
                == actual_open.structure().canonical().external_structure()
        );
        let source_trace = integrated::simplify(&closing_duals[0].clone().multiply(source)?)?;
        let source_trace = source_trace
            .with_checked_expression(graph.model.expand_couplings(source_trace.expression()))?;
        let expected_phase =
            Atom::i() * coupling.pow(4) * Atom::num((4, 3)).pow(2) * Atom::var(GS.dim).pow(2);
        Self::certify_parent_tensor_zero(
            "D09 independent zero-momentum native phase and CF-squared source",
            &source_trace.with_checked_expression(source_trace.expression() - &expected_phase)?,
        )?;

        let terms = FourDTerm::from_tensor(&retained_w)?;
        eyre::ensure!(terms.len() == 4);
        let expected_powers = [
            [1usize, 2, 1, 1, 1],
            [1, 2, 1, 1, 2],
            [1, 1, 1, 1, 1],
            [1, 1, 1, 1, 2],
        ];
        let k = function!(GS.emr_mom, 2);
        let l = function!(GS.emr_mom, 4);
        let rep = (Minkowski {}).new_rep(GS.dim).to_symbolic([]);
        let square = |q: &Atom| {
            function!(
                SPENSO_TAG.dot,
                GS.indexed_momentum(q, std::slice::from_ref(&rep)),
                GS.indexed_momentum(q, std::slice::from_ref(&rep))
            )
        };
        let mut descriptors = Vec::new();
        for (ordinal, term) in terms.iter().enumerate() {
            let mut powers = [0usize; 5];
            for denominator in &term.denominators {
                eyre::ensure!((2..=6).contains(&denominator.source_edge.0));
                powers[denominator.source_edge.0 - 2] += 1;
                let expected_momentum = match denominator.source_edge.0 {
                    2 => k.clone(),
                    3 => -&k,
                    _ => l.clone(),
                };
                eyre::ensure!(denominator.momentum == expected_momentum);
                eyre::ensure!(denominator.mass_squared == Atom::var(GS.m_uv_expansion).pow(2));
                let expected_polynomial = square(if denominator.source_edge.0 <= 3 {
                    &k
                } else {
                    &l
                }) - Atom::var(GS.m_uv_vacuum).pow(2);
                eyre::ensure!(
                    (denominator.polynomial_in_lmb(&lmb)? - expected_polynomial)
                        .expand()
                        .is_zero(),
                    "D09 W changed the complete vacuum polynomial or prescription"
                );
            }
            eyre::ensure!(powers == expected_powers[ordinal]);
            descriptors.push(serde_json::json!({"term":ordinal,"source_edge_powers":powers}));
        }

        let integrator = integrated::Integrated::new(crate::utils::vakint()?, backend);
        let mut traces = Vec::new();
        for physical_four in [false, true] {
            for (dual_index, dual) in closing_duals.iter().enumerate() {
                let mut sum = SymbolicTensor::infer_indexed(Atom::Zero)?;
                for term in &terms {
                    let numerator = prepare_numerator(&term.numerator, physical_four)?;
                    let closed = integrated::simplify(&dual.clone().multiply(numerator)?)?;
                    eyre::ensure!(closed.is_scalar());
                    let closed = if physical_four {
                        closed.expand_dots()?
                    } else {
                        closed
                    };
                    if physical_four {
                        let mut unsupported = false;
                        closed.expression().visitor(&mut |view| {
                            if let AtomView::Fun(momentum) = view
                                && momentum.get_symbol() == GS.emr_mom
                            {
                                unsupported |= momentum.get_nargs() != 2
                                    || !matches!(
                                        SlotMatcher::default().concrete_component(momentum.get(1)),
                                        Some(Ok(0..=3))
                                    );
                                return false;
                            }
                            true
                        });
                        eyre::ensure!(
                            !unsupported,
                            "physical-four numerator did not expose finite momentum components"
                        );
                    }
                    let denominator = Atom::mul_many(term.denominators.iter().map(|pole| {
                        function!(
                            GS.den,
                            pole.source_edge.0,
                            &pole.momentum,
                            &pole.mass_squared,
                            &pole.full_expr
                        )
                        .pow(-1)
                    }));
                    sum = sum.try_add(
                        &closed.with_checked_expression(closed.expression() * denominator)?,
                    )?;
                }
                let scalar_sector = FourDSector {
                    atom: sum.clone(),
                    ..sector.clone()
                }
                .with_scheme_integrand(sum);
                eyre::ensure!(scalar_sector.active_components == sector.active_components);
                eyre::ensure!(scalar_sector.frozen_lmbs == sector.frozen_lmbs);
                eyre::ensure!(scalar_sector.inherited_loop_edges == sector.inherited_loop_edges);
                eyre::ensure!(scalar_sector.provenance_components == sector.provenance_components);
                eyre::ensure!(scalar_sector.provenance_branches == sector.provenance_branches);
                let scalar_local = Local4dCts(FourDSectors::new(vec![scalar_sector], Vec::new())?);
                let result = integrator.run(&scalar_local, ctx, parent, given, parent, given)?;
                let laurent = result
                    .physical_pole()?
                    .try_sub(&result.physical_finite_counterterm()?)?;
                eyre::ensure!(laurent.is_scalar());
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_d09_proper_local_vacuum_scalar_trace", physical_four, dual_index,
                    file.trace = %cooking.uncook(laurent.expression().as_view()).to_canonical_string(),
                    source_absorption_proved = true, original_terms = terms.len(),
                    file.denominator_inventory = %serde_json::to_string(&descriptors)?,
                    negative_denominators_remain_dimensional = true,
                    finite_child_included = false, bare_included = false,
                    open_matrix_proved = false, emitted_R2_proved = false,
                    "Actual proper-local W scalar input uses original frames and the existing Integrated owner"
                );
                traces.push(laurent);
            }
        }
        eyre::ensure!(
            local == &retained_local && sector.integration_atom().same_value(&retained_w)
        );
        eyre::ensure!(open_basis.iter().all(|value| value.rank() == 4));
        traces
            .try_into()
            .map_err(|_| eyre!("D09 needs exactly four complete scalar traces"))
    }
}
