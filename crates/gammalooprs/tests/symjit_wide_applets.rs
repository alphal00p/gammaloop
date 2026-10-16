use symjit::{Application, Complex, Composer, Config, Defuns, Slot, Storage, Translator};

#[test]
fn retained_applets_use_call_arity_when_nested_and_reloaded()
-> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    for arity in [127, 128, 255, 256, 257, 513, 1023, 1024, 1025, 2047, 2048] {
        for direct in [false, true] {
            for complex in [false, true] {
                for optimization in [0, 3] {
                    let mut config = Config::default();
                    config.set_complex(complex);
                    config.set_direct(direct);
                    config.set_opt_level(optimization);
                    // SymJIT 2.27 defaults to 64 call slots; reserve the tested arity.
                    config.set_num_args(arity as u32);
                    let mut small = Translator::new(config.clone());
                    small.set_num_params(2);
                    small.append_mul(
                        &Slot::Out(0),
                        &[Slot::Param(0), Slot::Param(1)],
                        if complex { 0 } else { 2 },
                    )?;
                    let mut small_definitions = Defuns::new();
                    small_definitions.add_applet("small", small.compile()?.seal()?);
                    let mut inner_config = config.clone();
                    inner_config.set_defuns(small_definitions.clone());
                    let mut inner = Translator::new(inner_config);
                    inner.set_num_params(arity);
                    let parameters = (0..arity).map(Slot::Param).collect::<Vec<_>>();
                    inner.append_add(
                        &Slot::Temp(0),
                        &parameters,
                        if complex { 0 } else { arity },
                    )?;
                    inner.append_fun(
                        &Slot::Temp(1),
                        "small",
                        &[Slot::Param(0), Slot::Param(arity - 1)],
                        false,
                    )?;
                    inner.append_add(
                        &Slot::Out(0),
                        &[Slot::Temp(0), Slot::Temp(1)],
                        if complex { 0 } else { 2 },
                    )?;
                    let mut definitions = small_definitions;
                    definitions.add_applet("retained_wide", inner.compile()?.seal()?);
                    config.set_defuns(definitions);
                    let mut outer = Translator::new(config.clone());
                    outer.set_num_params(1);
                    let mut arguments = vec![Slot::Param(0)];
                    for i in 1..arity {
                        // Symbolica likewise indexes its appended constants directly.
                        outer.append_constant(Complex::new(i as f64, 0.0))?;
                        arguments.push(Slot::Const(i - 1));
                    }
                    outer.append_mul(
                        &Slot::Temp(0),
                        &[Slot::Param(0), Slot::Param(0)],
                        if complex { 0 } else { 2 },
                    )?;
                    outer.append_fun(&Slot::Temp(1), "retained_wide", &arguments, false)?;
                    outer.append_fun(
                        &Slot::Temp(2),
                        "small",
                        &[Slot::Temp(1), Slot::Param(0)],
                        false,
                    )?;
                    outer.append_add(
                        &Slot::Out(0),
                        &[Slot::Temp(0), Slot::Temp(2)],
                        if complex { 0 } else { 2 },
                    )?;
                    let evaluator = outer.compile()?;
                    let mut bytes = Vec::new();
                    evaluator.save(&mut bytes)?;
                    let restored = Application::load(&mut bytes.as_slice(), &config)?;
                    let sum = (arity * (arity - 1) / 2) as f64;
                    for evaluator in [evaluator, restored] {
                        if complex {
                            let result = evaluator.evaluate_single(&[Complex::new(3.0, 2.0)]);
                            assert_eq!(
                                result,
                                Complex::new(
                                    3.0 * sum + 5.0 * (arity + 1) as f64,
                                    2.0 * sum + 12.0 * (arity + 1) as f64
                                ),
                                "arity={arity} direct={direct} optimization={optimization}"
                            );
                        } else {
                            let result = evaluator.evaluate_single(&[3.0]);
                            assert_eq!(
                                result,
                                3.0 * sum + 9.0 * (arity + 1) as f64,
                                "arity={arity} direct={direct} optimization={optimization}"
                            );
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn retained_applet_argument_reservation_is_clamped() {
    for direct in [false, true] {
        for (requested, expected) in [(0, 64), (64, 64), (1025, 1025), (2048, 2048), (2049, 2048)] {
            let mut config = Config::default();
            config.set_direct(direct);
            assert_eq!(config.num_args(), 64);
            config.set_num_args(requested);
            assert_eq!(
                config.num_args(),
                expected,
                "direct={direct} requested={requested}"
            );
        }
    }
}

#[test]
fn indirect_retained_applets_reject_calls_beyond_reserved_capacity()
-> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    // Stock 2.27's direct translator silently accepts out-of-reservation calls;
    // never execute such a program. GammaLoop checks this boundary before JIT.
    for (arity, reservation) in [(65, 64), (1025, 1024), (2049, 2049)] {
        let mut config = Config::default();
        config.set_direct(false);
        config.set_num_args(reservation);
        let capacity = config.num_args();
        let mut inner = Translator::new(config.clone());
        inner.set_num_params(arity);
        inner.append_assign(&Slot::Out(0), &Slot::Param(arity - 1))?;
        let mut definitions = Defuns::new();
        definitions.add_applet("too_wide", inner.compile()?.seal()?);
        config.set_defuns(definitions);
        let mut outer = Translator::new(config);
        outer.set_num_params(arity);
        let arguments = (0..arity).map(Slot::Param).collect::<Vec<_>>();
        let failure = outer
            .append_fun(&Slot::Out(0), "too_wide", &arguments, false)
            .and_then(|()| outer.compile().map(|_| ()))
            .expect_err("indirect translation must reject a call above its reserved capacity");
        assert!(
            failure.to_string().contains(&format!(
                "function call with {arity} args, but only {capacity} are available"
            )),
            "unexpected rejection: {failure}"
        );
    }
    Ok(())
}
