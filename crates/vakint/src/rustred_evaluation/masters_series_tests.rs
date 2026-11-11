//! Frozen V21 scalar-substitution oracle and Laurent-series parity gates.
use super::*;

// Keep this oracle independent of the production series assembler.
impl MasterCatalog {
    fn scalar_reference_finalize(
        &self,
        raw: Atom,
        settings: &VakintSettings,
        loop_count: usize,
        mode: RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        let _stage = StageSpan::new("vakint.masters.finalize");
        // User normalization factors can introduce finite literals after the
        // frontend staged the original input. Bind those remaining literals
        // before rational algebra; an exact zero prefactor must not launch an
        // unbounded search for its first nonzero Laurent coefficient.
        let raw = if mode == RustRedMasterMode::Numerical {
            let _stage = StageSpan::new("vakint.masters.stage_finite_literals");
            crate::precision::legacy::stage_table(
                S.rustred_input.call(loop_count as i64).as_view(),
                raw.as_view(),
                vk_symbol!(&settings.epsilon_symbol),
                "Input",
            )?
        } else {
            raw
        };
        let direct =
            self.scalar_reference_finalize_independent(raw.clone(), settings, loop_count, mode);
        if matches!(&direct, Err(error) if !Self::lacks_master_data(error))
            || mode == RustRedMasterMode::Retain
            || self.basis_aliases.is_empty()
        {
            return direct;
        }
        // Prefer authenticated exact relations even when independently rounded
        // coefficient rows are available: a large prefactor must not amplify
        // their rounding differences in a mathematically zero combination.
        // Keep direct data as a fallback when an alias has additional orders.
        // Corruption, parsing and algebra errors are never hidden by this path.
        let alias_substitution = StageSpan::new("vakint.masters.alias_substitution");
        let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        let replacements: Vec<_> = self
            .basis_aliases
            .iter()
            .map(|(master, basis)| {
                let basis = basis
                    .replace(Atom::var(S.form_epsilon).to_pattern())
                    .with(epsilon.clone());
                Replacement::new(master.clone(), basis)
            })
            .collect();
        let normalized = raw.replace_multiple(&replacements);
        drop(alias_substitution);
        let normalized = {
            let _stage = StageSpan::new("vakint.masters.alias_combine_exact");
            combine_exact(normalized)
        };
        if normalized == raw {
            return direct;
        }
        match self.scalar_reference_finalize_independent(normalized, settings, loop_count, mode) {
            Err(error) if Self::lacks_master_data(&error) => direct,
            result => result,
        }
    }

    fn scalar_reference_finalize_independent(
        &self,
        raw: Atom,
        settings: &VakintSettings,
        loop_count: usize,
        mode: RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        let _stage = StageSpan::new("vakint.masters.finalize_independent");
        if mode == RustRedMasterMode::Retain {
            return Ok(self.retain(
                raw,
                &BTreeMap::new(),
                None,
                crate::rustred_evaluation::postprocessing_threads(settings),
            ));
        }
        // The scalar/tensor admission boundary has already checked the source
        // representations. Their dimensions label tensor interfaces rather
        // than scalar epsilon dependence, so hold that metadata fixed during
        // both valuation and Laurent expansion, then restore it on return.
        let cook = StageSpan::new("vakint.masters.cook");
        let cooking = idenso::CookSettings::reversible().with_representation_payloads(false, true);
        let raw = cooking.try_cook(raw.as_view()).map_err(failure)?;
        drop(cook);
        let raw = {
            let _stage = StageSpan::new("vakint.masters.initial_combine_exact");
            combine_exact(raw)
        };
        self.require_coefficients(&raw, &BTreeMap::new())?;
        let order_setup = StageSpan::new("vakint.masters.order_setup");
        let epsilon = vk_symbol!(&settings.epsilon_symbol);
        let target = settings
            .number_of_terms_in_epsilon_expansion
            .checked_sub(i64::try_from(loop_count).map_err(failure)?)
            .and_then(|value| value.checked_sub(1))
            .ok_or_else(|| failure("epsilon expansion depth overflow"))?;
        if !self.basis_aliases.contains_key(&raw)
            && let Some(master) = self.masters.get(&raw)
        {
            // Aliases retain their original validation order, so fallback
            // cannot hide malformed lower-order coefficient formulas.
            // A bare master needs its highest requested coefficient directly;
            // no other term or finite coefficient can cancel this dependency.
            // General expressions retain the exact combination below.
            let first = master.leading_power()?;
            if target >= first && !master.coefficients.contains_key(&target) {
                let known = master.known_through(first);
                return Err(crate::truncation_guard::MissingMasterOrders {
                    backend: self.guard_backend(),
                    requested_order: Some(target),
                    masters: vec![crate::truncation_guard::MissingMasterOrder {
                        master: self.describe_master(&raw),
                        first_missing_order: known + 1,
                        last_needed_order: Some(target),
                        first_affected_result_order: known + 1,
                    }],
                }
                .into());
            }
        }
        drop(order_setup);
        let symbols = {
            let _stage = StageSpan::new("vakint.masters.master_discovery");
            Self::referenced_masters(&raw)
                .into_iter()
                .collect::<Vec<_>>()
        };
        let terms = {
            let _stage = StageSpan::new("vakint.masters.coefficient_list");
            raw.coefficient_list::<i32>(&symbols)
        };
        let mut required = BTreeMap::<Atom, i64>::new();
        let mut known = Atom::Zero;
        for (monomial, coefficient) in terms {
            let coefficient = {
                let _stage = StageSpan::new("vakint.masters.coefficient_combine_exact");
                combine_exact(coefficient)
            };
            if coefficient.is_zero() {
                continue;
            }
            let valuation_timing = StageSpan::new("vakint.masters.coefficient_valuation");
            let mut multiplicities = BTreeMap::new();
            Self::master_multiplicities(monomial.as_view(), &mut multiplicities)?;
            let series = coefficient
                .series(epsilon, Atom::Zero.as_view(), SeriesDepth::relative(1))
                .map_err(failure)?;
            let exponent = series.get_trailing_exponent();
            if !exponent.is_integer() {
                return Err(failure("fractional epsilon power in reduction coefficient").into());
            }
            let valuation = exponent
                .numerator()
                .to_i64()
                .ok_or_else(|| failure("epsilon valuation exceeds integer range"))?;
            drop(valuation_timing);
            let order_planning = StageSpan::new("vakint.masters.coefficient_order_planning");
            let minimum = multiplicities
                .iter()
                .try_fold(0_i64, |sum, (symbol, power)| {
                    let leading = self.masters[symbol].leading_power()?;
                    sum.checked_add(
                        leading
                            .checked_mul(*power)
                            .ok_or_else(|| failure("master product order overflow"))?,
                    )
                    .ok_or_else(|| failure("master product order overflow"))
                })?;
            for symbol in multiplicities.keys() {
                let own_minimum = self.masters[symbol].leading_power()?;
                let last = target
                    .checked_sub(valuation)
                    .and_then(|value| value.checked_sub(minimum))
                    .and_then(|value| value.checked_add(own_minimum))
                    .ok_or_else(|| failure("required master epsilon order overflow"))?;
                required
                    .entry(symbol.clone())
                    .and_modify(|current| *current = (*current).max(last))
                    .or_insert(last);
            }
            drop(order_planning);
            let _stage = StageSpan::new("vakint.masters.known_coefficient_assembly");
            known += monomial * coefficient;
        }
        let mut exact_values = BTreeMap::new();
        for (symbol, last) in required {
            let _stage = StageSpan::new("vakint.masters.master_expansion");
            let master = &self.masters[&symbol];
            let first = master.leading_power()?;
            // Substitute the available coefficients the planned depth needs
            // and an explicit unknown tail ε^k·Oep(k, master) beyond them.
            // Missing data and any misplanned depth thus survive truncation
            // as `Oep` instead of becoming zero (see `truncation_guard`).
            let through = last.min(master.known_through(first));
            let mut expansion = crate::truncation_guard::unknown_tail(
                epsilon,
                through
                    .checked_add(1)
                    .ok_or_else(|| failure("required master epsilon order overflow"))?,
                &symbol,
            );
            for power in first..=through {
                let placeholder = self.coefficient_symbol(&symbol, power);
                if let Some(coefficient) = master.coefficients.get(&power)
                    && let Some(formula) = &coefficient.exact
                {
                    let exact = {
                        let _stage = StageSpan::new("vakint.masters.formula_parse");
                        vk_parse!(formula).map_err(failure)?
                    };
                    // Keep finite embedded constants behind a named coefficient
                    // until their precision has been checked in numerical mode.
                    let finite_precision = {
                        let _stage = StageSpan::new("vakint.masters.formula_precision");
                        Self::finite_precision(&exact)
                    };
                    if finite_precision.is_none() {
                        let exact = {
                            let _stage = StageSpan::new("vakint.masters.formula_normalization");
                            crate::hpl::normalize_constants(exact.as_view())?
                        };
                        // Numerical evaluation needs the original coefficient
                        // identity for precision diagnostics. Exact zeros must
                        // still disappear before epsilon availability checks.
                        exact_values.insert(placeholder.clone(), exact);
                    }
                }
                expansion += placeholder * Atom::var(epsilon).pow(Atom::num(power));
            }
            let _stage = StageSpan::new("vakint.masters.master_substitution");
            known = known.replace(symbol.to_pattern()).with(expansion);
        }
        // Exact zero coefficients and exact relations must cancel before an
        // absent/high-order or low-precision coefficient can cause an error.
        let coefficient_substitution = StageSpan::new("vakint.masters.coefficient_substitution");
        for (source, value) in &exact_values {
            if mode == RustRedMasterMode::Symbolic || value.is_zero() {
                known = known.replace(source.to_pattern()).with(value.clone());
            }
        }
        drop(coefficient_substitution);
        let result = {
            let _stage = StageSpan::new("vakint.masters.epsilon_series");
            known
                .series(epsilon, Atom::Zero.as_view(), Rational::from(target))
                .map_err(failure)?
                .to_atom()
        };
        let mut result = {
            let _stage = StageSpan::new("vakint.masters.final_combine_exact");
            combine_exact(result)
        };
        // Keeping exact coefficient identities improves diagnostics, but an
        // unavailable order can multiply a provably zero combination of known
        // lower orders. Resolve that algebraic question before rejecting it.
        let numerical_availability = StageSpan::new("vakint.masters.numerical_availability");
        if mode == RustRedMasterMode::Numerical
            && (crate::truncation_guard::has_unknown_tail(result.as_view())
                || matches!(
                    self.check_coefficients(&result, false),
                    Err(VakintError::RustRedEvaluation(
                        RustRedEvaluationError::MasterCapability { .. }
                    ))
                ))
        {
            let substitution = StageSpan::new("vakint.masters.numerical_exact_substitution");
            let replacements = exact_values
                .into_iter()
                .map(|(source, value)| Replacement::new(source, value))
                .collect::<Vec<_>>();
            result = result.replace_multiple(&replacements);
            drop(substitution);
            let _stage = StageSpan::new("vakint.masters.numerical_exact_combine_exact");
            result = combine_exact(result);
        }
        drop(numerical_availability);
        let validation = StageSpan::new("vakint.masters.remaining_validation");
        crate::truncation_guard::ensure_known_orders(
            &result,
            epsilon,
            Some(target),
            &self.guard_backend(),
            |master| self.describe_master(master),
        )?;
        self.check_coefficients(&result, mode == RustRedMasterMode::Numerical)?;
        drop(validation);
        if mode == RustRedMasterMode::Numerical {
            result = self.resolve_numerical(result, settings)?;
        }
        if mode == RustRedMasterMode::Symbolic {
            let _stage = StageSpan::new("vakint.masters.bind_period_references");
            // An exact formula can still depend on snapshot-local finite
            // periods. Bind those leaves before returning a symbolic result,
            // so loading another release cannot silently change its value.
            result = self.bind_period_references(result)?;
        }
        let _stage = StageSpan::new("vakint.masters.uncook");
        Ok(cooking.uncook(result.as_view()))
    }
}

fn settings(target: i64, loops: usize) -> VakintSettings {
    VakintSettings {
        epsilon_symbol: "ep".into(),
        number_of_terms_in_epsilon_expansion: target + loops as i64 + 1,
        use_dot_product_notation: true,
        ..Default::default()
    }
}

fn add_master(catalog: &mut MasterCatalog, name: &str, loops: usize, rows: &[(i64, &str)]) -> Atom {
    let symbol = catalog.symbol(name, &[1]);
    catalog.masters.insert(
        symbol.clone(),
        Master {
            identity: name.into(),
            loop_count: loops,
            graph: "series-parity.dot".into(),
            metadata: serde_json::Value::Null,
            data: None,
            parsed_basis: OnceLock::new(),
            parsed_names: OnceLock::new(),
            coefficients: rows
                .iter()
                .map(|&(power, value)| {
                    (
                        power,
                        Coefficient {
                            exact: Some(value.into()),
                            numerical: None,
                            significant_digits: None,
                            absolute_error: None,
                            error_is_bounded: false,
                        },
                    )
                })
                .collect(),
        },
    );
    symbol
}

fn check(
    catalog: &MasterCatalog,
    raw: Atom,
    settings: &VakintSettings,
    loops: usize,
    mode: RustRedMasterMode,
) -> Result<Atom, String> {
    let expected = catalog
        .scalar_reference_finalize(raw.clone(), settings, loops, mode)
        .map_err(|error| error.to_string());
    let actual = catalog
        .finalize(raw.clone(), settings, loops, mode)
        .map_err(|error| error.to_string());
    assert_eq!(actual, expected, "mode={mode:?}, input={raw}");
    actual
}

#[test]
fn master_series_matches_scalar_products_poles_and_leading_zeros() {
    let (mut catalog, _) = tests::catalog(Some("2"), 30);
    let a = add_master(
        &mut catalog,
        "series-a",
        1,
        &[
            (-1, "1"),
            (0, "2"),
            (1, "3"),
            (2, "5"),
            (3, "7"),
            (4, "11"),
            (5, "13"),
            (6, "17"),
        ],
    );
    let b = add_master(
        &mut catalog,
        "series-b",
        2,
        &[
            (-2, "0"),
            (-1, "0"),
            (0, "3"),
            (1, "4"),
            (2, "5"),
            (3, "6"),
            (4, "7"),
            (5, "8"),
            (6, "9"),
        ],
    );
    let eps = vk_parse!("ep").unwrap();
    let cases = [
        Atom::Zero,
        vk_parse!("ep^-3/(1-ep)").unwrap(),
        a.pow(2) / eps.pow(2),
        &a * &b / eps.pow(2),
        a.pow(2) * &b * vk_parse!("(1+ep)/(1-ep)").unwrap(),
        &a * vk_parse!("(1+𝑖)*2^(1/2)/(1-ep)").unwrap() + &b * vk_parse!("(2-𝑖)/(1+ep)").unwrap(),
        &a * vk_parse!("mUV^(-2*ep)*mu^(2*ep)/(1-ep)").unwrap(),
    ];
    for target in [-3, -1, 0, 2] {
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            for raw in &cases {
                let _ = check(&catalog, raw.clone(), &settings(target, 3), 3, mode).unwrap();
            }
        }
    }
    let result = check(
        &catalog,
        a.pow(2) / eps.pow(2),
        &settings(0, 3),
        3,
        RustRedMasterMode::Symbolic,
    )
    .unwrap();
    assert_eq!(
        result.expand(),
        vk_parse!("ep^-4+4*ep^-3+10*ep^-2+22*ep^-1+43").unwrap()
    );
}

#[test]
fn master_series_preserves_missing_orders_exact_cancellations_and_alias_errors() {
    let (mut catalog, a) = tests::catalog(Some("2"), 30);
    let b = add_master(&mut catalog, "series-alias", 1, &[(-1, "1"), (0, "2")]);
    catalog.basis_aliases.insert(b.clone(), a.clone());
    let unknown = vk_parse!("RustRedMaster(series_unknown,1)").unwrap();
    let zero = vk_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap();
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        for (raw, target, success) in [
            (&a + &zero * &unknown, 0, true),
            (&zero * &unknown, 2, true),
            (&a - &b, 3, true),
            (&a / vk_parse!("ep").unwrap(), 0, false),
            (a.clone(), 1, false),
            (&a + &unknown, 0, false),
            (&a * vk_parse!("ep^(1/2)").unwrap(), 0, false),
        ] {
            let result = check(&catalog, raw, &settings(target, 1), 1, mode);
            assert_eq!(result.is_ok(), success, "{result:?}");
        }
    }
    // A direct formula error must not be hidden by an otherwise valid alias.
    catalog
        .masters
        .get_mut(&b)
        .unwrap()
        .coefficients
        .get_mut(&0)
        .unwrap()
        .exact = Some("broken(".into());
    assert!(check(&catalog, b, &settings(0, 1), 1, RustRedMasterMode::Symbolic).is_err());
}

#[test]
fn master_series_keeps_cross_master_singular_cancellation_and_float_order() {
    let (mut catalog, _) = tests::catalog(Some("2"), 30);
    let a = add_master(
        &mut catalog,
        "singular-a",
        1,
        &[(-1, "1"), (0, "exp(1/ep)")],
    );
    let b = add_master(
        &mut catalog,
        "singular-b",
        1,
        &[(-1, "-1"), (0, "-exp(1/ep)")],
    );
    // The singular functions cancel after ordered substitution, before any
    // series is requested. Expanding either master separately would fail.
    assert_eq!(
        check(
            &catalog,
            &a + &b,
            &settings(0, 1),
            1,
            RustRedMasterMode::Symbolic
        )
        .unwrap(),
        Atom::Zero
    );
    let a = add_master(&mut catalog, "float-a", 1, &[(-1, "1"), (0, "2")]);
    let b = add_master(&mut catalog, "float-b", 1, &[(-1, "1"), (0, "2")]);
    let c = add_master(&mut catalog, "float-c", 1, &[(-1, "1"), (0, "2")]);
    let large = vk_parse!("1000000000000000000000000000000.0").unwrap();
    let raw = &large * a + vk_parse!("0.1").unwrap() * b - large * c;
    let _ = check(
        &catalog,
        raw,
        &settings(0, 1),
        1,
        RustRedMasterMode::Symbolic,
    )
    .unwrap();
}

#[test]
fn master_series_preserves_dimension_metadata_and_embedded_master_references() {
    use spenso::{
        network::library::symbolic::ETS,
        structure::representation::{BaseRepName, Minkowski, initialize},
    };
    let (catalog, a) = tests::catalog(Some("2"), 30);
    initialize();
    let dimension = vk_parse!("4-2*ep").unwrap();
    let slot = |index| Minkowski::selfless_symbol().call((&dimension, index));
    let metric = ETS.metric(slot(vk_symbol!("series_mu")), slot(vk_symbol!("series_nu")));
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let _ = check(
            &catalog,
            &a * &metric / dimension.clone(),
            &settings(0, 1),
            1,
            mode,
        )
        .unwrap();
    }
    // A coefficient can contain a master below an opaque function. Its poles
    // must still trigger the original series error, never stay unexpanded.
    let raw = vk_symbol!("series_opaque").call((&a,)) + a;
    assert!(
        check(
            &catalog,
            raw,
            &settings(0, 1),
            1,
            RustRedMasterMode::Symbolic
        )
        .is_err()
    );
}

#[test]
#[ignore = "captured four-loop proper-UV expression; compares complete scalar and series finalizers"]
fn master_series_matches_frozen_scalar_qcd_fixture() {
    use std::io::Read;
    let engine = crate::Vakint::new().unwrap();
    let data = engine.rustred_data().unwrap();
    let catalog = &data.masters;
    let provenance: serde_json::Value = serde_json::from_slice(include_bytes!(
        "../../tests/resources/rustred/qcd-l4-358-master-finalization.expr.provenance.json"
    ))
    .unwrap();
    let compressed =
        include_bytes!("../../tests/resources/rustred/qcd-l4-358-master-finalization.expr.gz");
    assert_eq!(
        format!("{:x}", Sha256::digest(compressed)),
        provenance["compressed"]["sha256"].as_str().unwrap()
    );
    let mut input = String::new();
    flate2::read::GzDecoder::new(compressed.as_slice())
        .read_to_string(&mut input)
        .unwrap();
    assert_eq!(
        format!("{:x}", Sha256::digest(input.as_bytes())),
        provenance["decoded"]["sha256"].as_str().unwrap()
    );
    let raw = vk_parse!(&input).unwrap();
    assert_eq!(raw.as_view().nterms(), 1050);
    let settings = VakintSettings {
        epsilon_symbol: "gammalooprs::ε".into(),
        ..settings(0, 3)
    };
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let _ = check(catalog, raw.clone(), &settings, 3, mode).unwrap();
    }
}

#[test]
#[ignore = "replay the 81 captured exact-algebra inputs against frozen scalar finalization"]
fn master_series_replays_captured_inputs_against_scalar_finalizer() {
    let engine = crate::Vakint::new().unwrap();
    let data = engine.rustred_data().unwrap();
    let catalog = &data.masters;
    let directory = std::env::var("VAKINT_COMBINE_CAPTURE_DIR")
        .expect("set VAKINT_COMBINE_CAPTURE_DIR to the frozen 81 canonical .expr payloads");
    let mut paths = std::fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "expr")
        })
        .collect::<Vec<_>>();
    paths.sort();
    assert_eq!(paths.len(), 81);
    let mut succeeded = 0;
    for path in paths {
        let input = std::fs::read_to_string(&path).unwrap();
        let raw = vk_parse!(&input).unwrap();
        let epsilon = if raw.contains_symbol(vk_symbol!("gammalooprs::ε")) {
            "gammalooprs::ε"
        } else {
            "ep"
        };
        let settings = VakintSettings {
            epsilon_symbol: epsilon.into(),
            ..settings(0, 3)
        };
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            succeeded += usize::from(check(catalog, raw.clone(), &settings, 3, mode).is_ok());
        }
    }
    assert!(
        succeeded > 0,
        "the replay must exercise successful complete finalization"
    );
}
