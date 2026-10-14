use symjit::{Application, Complex, Composer, Config, Defuns, Slot, Storage, Translator};

#[test]
fn retained_applets_use_call_arity_when_nested_and_reloaded()
-> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    for arity in [127, 128, 255, 256, 257, 513, 1023, 1024] {
        for direct in [false, true] {
            for complex in [false, true] {
                for optimization in [0, 3] {
                    let mut config = Config::default();
                    config.set_complex(complex);
                    config.set_direct(direct);
                    config.set_opt_level(optimization);
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
fn retained_applets_reject_calls_beyond_the_upstream_argument_limit() {
    // Stock SymJIT supports 1024 retained-call arguments, including hoisted
    // constants. Check both frontends reject a larger signature explicitly.
    for direct in [false, true] {
        let failure = std::panic::catch_unwind(|| {
            let mut config = Config::default();
            config.set_direct(direct);
            let mut inner = Translator::new(config.clone());
            inner.set_num_params(1025);
            let arguments = (0..1025).map(Slot::Param).collect::<Vec<_>>();
            inner.append_add(&Slot::Out(0), &arguments, 1025).unwrap();
            let mut definitions = Defuns::new();
            definitions.add_applet("too_wide", inner.compile().unwrap().seal().unwrap());
            config.set_defuns(definitions);
            let mut outer = Translator::new(config);
            outer.set_num_params(1025);
            outer
                .append_fun(&Slot::Out(0), "too_wide", &arguments, false)
                .unwrap();
            outer.compile().unwrap();
        })
        .expect_err("upstream must reject a retained call above its 1024-argument limit");
        let message = failure
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| failure.downcast_ref::<&str>().copied())
            .unwrap_or("");
        assert!(
            message.contains("n <= SLICE_CAP"),
            "unexpected rejection: {message}"
        );
    }
}
