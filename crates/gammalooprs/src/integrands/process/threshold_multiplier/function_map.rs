use std::collections::{BTreeMap, BTreeSet};

use color_eyre::eyre::{Context, Result, eyre};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Indeterminate, Symbol},
    evaluate::FunctionMap,
    symbol,
};

use crate::integrands::process::param_builder::FnMapEntry;

/// Scalar definitions shared by multiplier expressions before binding their cut-local inputs.
pub(super) struct ThresholdMultiplierFunctions {
    entries: Vec<FnMapEntry>,
}

impl ThresholdMultiplierFunctions {
    pub(super) fn parse(
        definitions: &BTreeMap<String, String>,
        parse_atom: impl Fn(&str) -> Result<Atom>,
    ) -> Result<Self> {
        let mut entries = Vec::<FnMapEntry>::new();
        let default = (
            "dE(m_a,m_b,Ea,Eb,ax,ay,az,bx,by,bz)",
            "((m_a-m_b)*(m_a+m_b)+(ax-bx)*(ax+bx)+(ay-by)*(ay+by)+(az-bz)*(az+bz))/(Ea+Eb)",
        );
        for (definition_index, (signature, body)) in definitions
            .iter()
            .map(|(signature, body)| (signature.as_str(), body.as_str()))
            .chain(std::iter::once(default))
            .enumerate()
        {
            let header = parse_atom(signature)
                .with_context(|| format!("invalid threshold function signature `{signature}`"))?;
            let (name, arguments) = match header.as_view() {
                AtomView::Var(variable) => (variable.get_symbol(), Vec::new()),
                AtomView::Fun(function) => (
                    function.get_symbol(),
                    function
                        .iter()
                        .map(|argument| match argument {
                            AtomView::Var(variable) => Ok(variable.get_symbol()),
                            _ => Err(eyre!(
                                "threshold function `{signature}` requires distinct scalar argument names; tags and static edge/frame formals are not supported"
                            )),
                        })
                        .collect::<Result<Vec<_>>>()?,
                ),
                _ => {
                    return Err(eyre!(
                        "threshold function signature `{signature}` must be a name or scalar function call"
                    ));
                }
            };
            if Self::reserved(name) {
                return Err(eyre!("reserved threshold function name `{name}`"));
            }
            if entries
                .iter()
                .any(|entry| entry.lhs.as_fun_view().unwrap().get_symbol() == name)
            {
                if definition_index == definitions.len() {
                    // An explicit user definition replaces the default, without depending on
                    // Symbolica's first-registration-wins behavior.
                    continue;
                }
                return Err(eyre!(
                    "duplicate threshold function head `{name}` in `{signature}`"
                ));
            }
            for (index, argument) in arguments.iter().enumerate() {
                if arguments[..index].contains(argument) || Self::reserved(*argument) {
                    return Err(eyre!(
                        "duplicate or reserved argument `{argument}` in threshold function `{signature}`"
                    ));
                }
            }

            // Private formal symbols keep a caller's argument from shadowing a captured
            // graph input in another function, in both Symbolica translation modes.
            let formals = arguments
                .iter()
                .enumerate()
                .map(|(index, _)| {
                    symbol!(&format!(
                        "gammalooprs::threshold_multiplier_arg::{}::arg_{index}",
                        name.get_name()
                    ))
                })
                .collect::<Vec<_>>();
            let body = parse_atom(body)
                .with_context(|| format!("invalid body of threshold function `{signature}`"))?;
            let body = body.replace_map(|view, _, output| {
                if let AtomView::Var(variable) = view
                    && let Some(index) = arguments
                        .iter()
                        .position(|argument| *argument == variable.get_symbol())
                {
                    **output = Atom::var(formals[index]);
                }
            });
            entries.push(FnMapEntry {
                lhs: FunctionBuilder::new(name)
                    .add_args(formals.iter().map(|formal| Atom::var(*formal)))
                    .finish(),
                rhs: body,
                args: formals.into_iter().map(Indeterminate::from).collect(),
                tags: Vec::new(),
            });
        }
        Ok(Self { entries })
    }

    pub(super) fn bind(
        &self,
        source: Atom,
        normalize: impl Fn(&Atom) -> Result<Atom>,
    ) -> Result<(Atom, FunctionMap, Vec<FnMapEntry>)> {
        let mut entries = Vec::new();
        let source = self.bind_expression(
            &source,
            &normalize,
            &mut vec![0; self.entries.len()],
            &mut entries,
        )?;
        let mut functions = FunctionMap::new();
        for entry in &entries {
            functions
                .add_function(
                    entry.lhs.as_fun_view().unwrap().get_symbol(),
                    entry.args.clone(),
                    entry.rhs.clone(),
                )
                .map_err(|error| eyre!("failed to bind threshold function: {error}"))?;
        }
        Ok((source, functions, entries))
    }

    fn bind_expression(
        &self,
        source: &Atom,
        normalize: &impl Fn(&Atom) -> Result<Atom>,
        states: &mut [u8],
        entries: &mut Vec<FnMapEntry>,
    ) -> Result<Atom> {
        let mut dependencies = BTreeSet::new();
        let mut error = None;
        let source = source.replace_map(|view, _, output| {
            let (name, arity) = match view {
                AtomView::Var(variable) => (variable.get_symbol(), 0),
                AtomView::Fun(function) => (function.get_symbol(), function.get_nargs()),
                _ => return,
            };
            let Some((index, entry)) = self
                .entries
                .iter()
                .enumerate()
                .find(|(_, entry)| entry.lhs.as_fun_view().unwrap().get_symbol() == name)
            else {
                return;
            };
            if arity != entry.args.len() {
                error = Some(eyre!(
                    "threshold function `{name}` expects {} arguments, got {arity}",
                    entry.args.len()
                ));
                return;
            }
            dependencies.insert(index);
            if matches!(view, AtomView::Var(_)) {
                **output = FunctionBuilder::new(name).finish();
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        for index in dependencies {
            let entry = &self.entries[index];
            match states[index] {
                1 => {
                    return Err(eyre!(
                        "cyclic threshold function dependency at `{}`",
                        entry.lhs
                    ));
                }
                2 => continue,
                _ => {}
            }
            states[index] = 1;
            let rhs = self.bind_expression(&entry.rhs, normalize, states, entries)?;
            entries.push(FnMapEntry {
                rhs,
                ..entry.clone()
            });
            states[index] = 2;
        }
        normalize(&source)
    }

    fn reserved(symbol: Symbol) -> bool {
        let name = symbol.get_name().rsplit("::").next().unwrap();
        symbol.is_builtin()
            || name.starts_with("threshold_multiplier_")
            || matches!(
                name,
                "E" | "Q3" | "Q" | "P" | "K" | "eta" | "eset" | "cind" | "effective" | "star"
            )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{parse, try_parse};

    #[test]
    fn nested_aliases_and_scalar_functions_capture_inputs_without_shadowing() {
        crate::initialisation::test_initialise().unwrap();
        let functions = ThresholdMultiplierFunctions::parse(
            &BTreeMap::from([
                ("captured".into(), "x+1".into()),
                ("nested(x)".into(), "x*captured".into()),
                ("outer(y)".into(), "nested(y+1)+captured()".into()),
            ]),
            |source| try_parse!(source).map_err(|error| eyre!("{error}")),
        )
        .unwrap();
        let (source, map, entries) = functions
            .bind(parse!("outer(2)+captured"), |atom| Ok(atom.clone()))
            .unwrap();
        assert_eq!(entries.len(), 3);
        assert!(matches!(entries[0].lhs.as_view(), AtomView::Fun(_)));
        for direct in [false, true] {
            let mut evaluator = source
                .evaluator(&[parse!("x")])
                .function_map(map.clone())
                .direct_translation(direct)
                .build()
                .unwrap()
                .map_coeff(&|coefficient| coefficient.re.to_f64());
            assert_eq!(evaluator.evaluate_single(&[5.0]), 30.0);
            assert_eq!(evaluator.evaluate_single(&[7.0]), 40.0);
        }
    }

    #[test]
    fn only_reachable_functions_bind_cut_geometry() {
        crate::initialisation::test_initialise().unwrap();
        let functions = ThresholdMultiplierFunctions::parse(
            &BTreeMap::from([
                ("local".into(), "2".into()),
                ("foreign".into(), "eta(star,eset(7,8))".into()),
                ("unreachable_cycle".into(), "unreachable_cycle()".into()),
            ]),
            |source| try_parse!(source).map_err(|error| eyre!("{error}")),
        )
        .unwrap();
        let normalize = |atom: &Atom| {
            if atom.contains_symbol(symbol!("eta")) {
                Err(eyre!("foreign eta is absent from this cut"))
            } else {
                Ok(atom.clone())
            }
        };
        let (_, _, entries) = functions.bind(parse!("local"), normalize).unwrap();
        assert_eq!(entries.len(), 1);
        assert!(functions.bind(parse!("foreign"), normalize).is_err());
        assert!(
            functions
                .bind(parse!("unreachable_cycle"), normalize)
                .is_err()
        );
    }

    #[test]
    fn invalid_definitions_and_calls_fail_before_evaluation() {
        crate::initialisation::test_initialise().unwrap();
        for definitions in [
            BTreeMap::from([("f(a)".into(), "a".into()), ("f(b)".into(), "b".into())]),
            BTreeMap::from([("f".into(), "1".into()), ("f()".into(), "2".into())]),
            BTreeMap::from([("f(a,a)".into(), "a".into())]),
            BTreeMap::from([("f(1)".into(), "1".into())]),
            BTreeMap::from([("E(a)".into(), "a".into())]),
            BTreeMap::from([("sin(a)".into(), "a".into())]),
            BTreeMap::from([("f(star)".into(), "star".into())]),
        ] {
            assert!(
                ThresholdMultiplierFunctions::parse(&definitions, |source| {
                    try_parse!(source).map_err(|error| eyre!("{error}"))
                })
                .is_err()
            );
        }
        let functions = ThresholdMultiplierFunctions::parse(
            &BTreeMap::from([
                ("f(x)".into(), "x".into()),
                ("a".into(), "b()".into()),
                ("b()".into(), "a".into()),
                ("unknown".into(), "missing_input+missing_function(1)".into()),
            ]),
            |source| try_parse!(source).map_err(|error| eyre!("{error}")),
        )
        .unwrap();
        for source in ["f", "f()", "f(1,2)", "a"] {
            assert!(
                functions
                    .bind(parse!(source), |atom| Ok(atom.clone()))
                    .is_err()
            );
        }
        let (source, map, _) = functions
            .bind(parse!("unknown"), |atom| Ok(atom.clone()))
            .unwrap();
        assert!(
            source
                .evaluator(&[] as &[Atom])
                .function_map(map)
                .build()
                .is_err()
        );
    }

    #[test]
    fn default_energy_difference_and_explicit_override_are_composable() {
        crate::initialisation::test_initialise().unwrap();
        for definitions in [
            BTreeMap::new(),
            BTreeMap::from([("dE(a,b)".into(), "a-b".into())]),
        ] {
            let functions = ThresholdMultiplierFunctions::parse(&definitions, |source| {
                try_parse!(source).map_err(|error| eyre!("{error}"))
            })
            .unwrap();
            let source = if definitions.is_empty() {
                // The first energy has mass 3 and momentum (4,0,0); the second is
                // at rest with mass/energy 2. This also checks unequal masses.
                parse!("dE(3,2,5,2,4,0,0,0,0,0)")
            } else {
                parse!("dE(5,2)")
            };
            let (source, map, entries) = functions.bind(source, |atom| Ok(atom.clone())).unwrap();
            assert_eq!(entries.len(), 1);
            let mut evaluator = source
                .evaluator(&[] as &[Atom])
                .function_map(map)
                .build()
                .unwrap()
                .map_coeff(&|coefficient| coefficient.re.to_f64());
            assert_eq!(evaluator.evaluate_single(&[]), 3.0);
        }
    }
}
