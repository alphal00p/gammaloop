use super::*;

#[test]
fn contracted_sixteen_trace_matches_clifford_components_and_dimension_specialization() {
    use super::super::trace_kernel::tests::clifford_trace;
    use spenso::network::library::symbolic::ETS;

    let r = test_initialize();
    let settings = GammaSimplifySettings::default();
    let mut outputs = Vec::new();
    for mink in [&r.mink4, &r.mink_d] {
        let p = momenta(&mink.to_symbolic([]));
        let slots: [_; 4] = std::array::from_fn(|i| mink.pattern(Atom::num(100 + i)));
        let [a, b, c, d] = &slots;
        let input = trace!(r.bis4.to_symbolic([]); [
            a, &p[0], b, &p[1], c, &p[2], d, &p[3],
            a, &p[4], d, &p[5], c, &p[6], b, &p[7],
        ].map(|argument| gamma!(argument)));
        let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        assert_ne!(result, input);
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((result).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            result
        );
        result.visitor(&mut |node| {
            assert!(!slots.iter().any(|slot| slot.as_view() == node));
            true
        });
        outputs.push(result);
    }
    let dimension = minkowski_dimension(r.mink_d.to_symbolic([]).as_view())
        .unwrap()
        .to_owned();
    let specialized = outputs[1]
        .replace(dimension.to_pattern())
        .with(Atom::num(4).to_pattern())
        .expand();
    assert_eq!(specialized, outputs[0].expand());

    // Reuse the independent blade-product oracle, specializing every metric
    // to a Euclidean orthonormal 4D basis. No trace recurrence is used here.
    let p = momenta(&r.mink4.to_symbolic([]));
    let mut seed = 71u64;
    let mut nonzero = false;
    for _ in 0..3 {
        let vectors: [[i64; 4]; 8] = std::array::from_fn(|_| {
            std::array::from_fn(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((seed >> 32) % 3) as i64 - 1
            })
        });
        let mut expected = 0;
        for assignment in 0usize..256 {
            let basis: [[i64; 4]; 4] = std::array::from_fn(|slot| {
                std::array::from_fn(|axis| i64::from(axis == ((assignment >> (2 * slot)) & 3)))
            });
            expected += clifford_trace(
                &[
                    basis[0], vectors[0], basis[1], vectors[1], basis[2], vectors[2], basis[3],
                    vectors[3], basis[0], vectors[4], basis[3], vectors[5], basis[2], vectors[6],
                    basis[1], vectors[7],
                ],
                false,
            );
        }
        nonzero |= expected != 0;
        let actual = outputs[0].replace_map(|node, _, out| {
            if let AtomView::Fun(metric) = node
                && metric.get_symbol() == ETS.metric
            {
                let indices: Vec<_> = metric
                    .iter()
                    .map(|argument| p.iter().position(|p| p.as_view() == argument).unwrap())
                    .collect();
                assert_eq!(indices.len(), 2);
                **out = Atom::num(
                    vectors[indices[0]]
                        .iter()
                        .zip(vectors[indices[1]])
                        .map(|(a, b)| a * b)
                        .sum::<i64>(),
                );
            }
        });
        assert_eq!(actual, Atom::num(expected));
    }
    assert!(
        nonzero,
        "the component oracle must not only test zero traces"
    );
}

#[test]
fn contracted_sixteen_mixed_components_keep_normalizer_semantics() {
    fn component(value: AtomView<'_>, out: &mut Settable<Atom>) {
        if let AtomView::Fun(vector) = value
            && let Some(AtomView::Fun(slot)) = vector.iter().last()
            && slot.get_symbol() == *MINKOWSKI_SYMBOL
            && slot.get_nargs() == 2
        {
            **out = Atom::num(1);
        }
    }
    let r = test_initialize();
    let rep = r.mink_d.to_symbolic([]);
    let p = symbolica::function!(
        spenso::vector_symbol!("idenso::contracted_sixteen_callback_p", norm = component),
        &rep
    );
    let q = symbolica::function!(
        spenso::vector_symbol!("idenso::contracted_sixteen_callback_q", norm = component),
        &rep
    );
    let a = slot!(r.mink_d, contracted_sixteen_a).into_atom();
    let b = slot!(r.mink_d, contracted_sixteen_b).into_atom();
    let c = slot!(r.mink_d, contracted_sixteen_c).into_atom();
    let input = trace!(r.bis4.to_symbolic([]);
        [&a, &a, &b, &p, &c, &q].into_iter()
            .chain(std::iter::repeat_n(&p, 10))
            .map(|argument| gamma!(argument)));
    let dimension = mink_slot_dimension(a.as_view()).unwrap();
    // The surviving components p(b), p(c), q(b), q(c) all normalize to one.
    // Treating their pre-normalization metrics as polynomial variables fails.
    let expected = (dimension
        * g!(&p, &p).pow(5).as_view()
        * (Atom::num(8) - Atom::num(4) * g!(&b, &c) * g!(&p, &q)))
    .expand();
    let settings = GammaSimplifySettings::default();
    let source = SymbolicTensor::<PartialStructure>::infer(input.clone()).unwrap();
    let result = DiracSimplifier::new(&settings)
        .evaluate_terminal_trace::<false>(input.as_view())
        .unwrap();
    assert_eq!(result.expand(), expected.expand());
    // The callback deliberately mixes a scalar with surviving b,c ports.
    // Keep its raw kernel oracle while rejecting the changed typed interface.
    assert!(source.simplify_gamma(settings).is_err());
    assert!(source.with_rewritten_expression(result).is_err());
}

#[test]
fn contracted_sixteen_unsupported_coefficients_keep_atom_fallback() {
    test_initialize();
    let settings = GammaSimplifySettings::default();
    for (dimension, unit) in [
        (
            Atom::num(4),
            Atom::var(symbolica::symbol!("contracted_sixteen_spin")),
        ),
        (Atom::num(4), Atom::num(1) / Atom::num(2)),
        (Atom::num(6), Atom::num(4)),
    ] {
        let rep = symbolica::function!(*MINKOWSKI_SYMBOL, &dimension);
        let a = symbolica::function!(
            *MINKOWSKI_SYMBOL,
            &dimension,
            Atom::var(symbolica::symbol!("contracted_sixteen_fallback_a"))
        );
        let p = &momenta(&rep)[0];
        let spin = symbolica::function!(*BISPINOR_SYMBOL, &unit);
        let input = trace!(&spin;
            [&a, &a].into_iter().chain(std::iter::repeat_n(p, 14))
                .map(|argument| gamma!(argument)));
        let expected = &unit * &dimension * g!(p, p).pow(7);
        // These raw trace units intentionally exceed typed gamma admission.
        // Keep the explicit sparse oracle and its Atom-emission fallback covered
        // through the same local trace evaluator, without a public raw route.
        let result = DiracSimplifier {
            settings: &settings,
            output: trace_kernel::TraceOutput::Expanded,
        }
        .evaluate_terminal_trace::<false>(input.as_view())
        .expect("the local trace owner accepts these coefficient domains");
        assert_eq!(result.expand(), expected.expand());
        assert_eq!(
            settings.rewrite_expression(result.clone(), trace_kernel::TraceOutput::Expanded),
            result
        );
    }
}

#[test]
fn sixteen_trace_preserves_open_words_and_metadata() {
    let r = test_initialize();
    let rep = r.mink_d.to_symbolic([]);
    let metadata = symbolica::parse_lit!((contracted_sixteen_x + contracted_sixteen_y) ^ 3);
    let vectors: Vec<_> = (0..16)
        .map(|i| {
            FunctionBuilder::new(T.rank_one_tensor_symbol("idenso::contracted_sixteen_metadata"))
                .add_arg(&metadata)
                .add_arg(Atom::num(i))
                .add_arg(&rep)
                .finish()
        })
        .collect();
    let start = slot!(r.bis4, contracted_sixteen_start).into_atom();
    let end = slot!(r.bis4, contracted_sixteen_end).into_atom();
    let open = chain!(&start, &end; vectors.iter().map(|p| gamma!(p)));
    let closed = trace!(r.bis4.to_symbolic([]); vectors.iter().map(|p| gamma!(p)));
    let settings = GammaSimplifySettings::default();
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((open).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        open
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((closed).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(settings.without_trace_evaluation())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        closed
    );
}
