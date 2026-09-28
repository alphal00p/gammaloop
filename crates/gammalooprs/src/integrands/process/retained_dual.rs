// This module deliberately depends only on std and Symbolica: standalone exports
// embed the same lowering code as the runtime evaluator builder.
use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    sync::atomic::{AtomicUsize, Ordering},
};

use symbolica::{
    atom::{EvaluationInfo, FunctionBuilder, Indeterminate, SymbolBuilder},
    domains::{dual::HyperDual, float::Complex},
    evaluate::{
        Dualizer, ExpressionEvaluator, FunctionMap, FunctionRegistrationOptions, InliningPolicy,
        InstructionList, OptimizationSettings, VectorInstruction, Vectorize,
    },
    prelude::*,
};

type Number = Complex<Rational>;

#[derive(Clone)]
pub(crate) struct RetainedFunctionDefinition {
    pub lhs: Atom,
    pub rhs: Atom,
    pub tags: Vec<Atom>,
    pub args: Vec<Indeterminate>,
    pub retained: bool,
}

// Resolve ordinary wrappers at their lexical arguments. A numerator body is
// resolved only while constructing its family, never at its residue call sites.
fn resolve_definitions(
    source: &Atom,
    definitions: &[RetainedFunctionDefinition],
    bound: &[Atom],
    stop_at_retained: bool,
) -> Result<Atom, String> {
    let mut result = source.clone();
    for _ in 0..=definitions.len() {
        let next = result.replace_map(|view, _, out| {
            if bound.iter().any(|parameter| parameter.as_view() == view) {
                out.set_from_view(&view);
                return;
            }
            let Some(definition) = definitions.iter().find(|definition| {
                if stop_at_retained && definition.retained {
                    return false;
                }
                match (definition.lhs.as_view(), view) {
                    (AtomView::Fun(lhs), AtomView::Fun(call)) => {
                        lhs.get_symbol() == call.get_symbol()
                            && call.get_nargs() == definition.tags.len() + definition.args.len()
                            && call
                                .iter()
                                .zip(&definition.tags)
                                .all(|(a, b)| a == b.as_view())
                    }
                    (AtomView::Var(_), _) => definition.lhs.as_view() == view,
                    _ => false,
                }
            }) else {
                return;
            };
            let rhs = if let AtomView::Fun(call) = view {
                let replacements = definition
                    .args
                    .iter()
                    .zip(call.iter().skip(definition.tags.len()))
                    .map(|(formal, actual)| {
                        Replacement::new(formal.as_view().to_pattern(), actual.to_owned())
                    })
                    .collect::<Vec<_>>();
                definition.rhs.replace_multiple(&replacements)
            } else {
                definition.rhs.clone()
            };
            out.set_from_view(&rhs.as_view());
        });
        if next == result {
            return Ok(result);
        }
        result = next;
    }
    Err("Cyclic function-map definitions while lowering retained dual functions".into())
}

struct Family {
    definition: RetainedFunctionDefinition,
    body: Atom,
    captures: Vec<Atom>,
}

#[derive(Clone, PartialEq, Eq, Hash)]
struct Request {
    family: usize,
    zeros: Vec<(usize, usize)>,
}

struct SharedDualizer<'a> {
    scalar: Dualizer<HyperDual<Number>>,
    shape: &'a [Vec<usize>],
    families: &'a [Family],
    requests: &'a HashMap<Symbol, Request>,
    function_map: &'a FunctionMap,
    optimization: &'a OptimizationSettings,
    scope: usize,
    cache: RefCell<HashMap<Request, ExpressionEvaluator<Number>>>,
}

impl Vectorize<Number> for SharedDualizer<'_> {
    fn duplicate_constants(&self) -> bool {
        false
    }

    fn get_dimension(&self) -> usize {
        self.shape.len()
    }

    fn map_instruction(
        &self,
        instruction: &VectorInstruction,
        instructions: &mut InstructionList<Number>,
    ) -> Vec<VectorInstruction> {
        self.scalar.map_instruction(instruction, instructions)
    }

    fn vectorize_function(
        &self,
        symbol: Symbol,
        tags: &[Atom],
        argument_count: usize,
    ) -> Result<Option<ExpressionEvaluator<Number>>, String> {
        let Some(request) = self.requests.get(&symbol) else {
            return self.scalar.vectorize_function(symbol, tags, argument_count);
        };
        if let Some(evaluator) = self.cache.borrow().get(request) {
            return Ok(Some(evaluator.clone()));
        }
        let family = &self.families[request.family];
        let dim = self.shape.len();
        let parameter_symbol = symbol!("gammalooprs::retained_dual_argument");
        let parameters = (0..argument_count * dim)
            .map(|index| function!(parameter_symbol, self.scope, symbol.get_id(), index))
            .collect::<Vec<_>>();
        let scalar_arguments = family
            .definition
            .args
            .iter()
            .map(|argument| argument.as_view().to_owned())
            .chain(family.captures.iter().cloned())
            .collect::<Vec<_>>();
        // The final argument anchors constant calls to an ordinary caller input.
        // Its components are deliberately absent from the function body.
        if argument_count != scalar_arguments.len() + 1 {
            return Err("Retained dual function has an inconsistent closure signature".into());
        }
        let replacements = scalar_arguments
            .iter()
            .zip(parameters.chunks(dim))
            .map(|(argument, components)| {
                Replacement::new(argument.to_pattern(), components[0].clone())
            })
            .collect::<Vec<_>>();
        let mut components = vec![Atom::Zero; dim];
        components[0] = family.body.replace_multiple(&replacements);
        let indices = self
            .shape
            .iter()
            .enumerate()
            .map(|(index, powers)| (powers.clone(), index))
            .collect::<HashMap<_, _>>();
        let zeros = request.zeros.iter().copied().collect::<HashSet<_>>();
        let variables = parameters
            .iter()
            .cloned()
            .map(Indeterminate::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| error.to_string())?;
        let mut order = (1..dim).collect::<Vec<_>>();
        order.sort_by_key(|&index| self.shape[index].iter().sum::<usize>());
        for component in order {
            let alpha = &self.shape[component];
            let axis = alpha
                .iter()
                .position(|&power| power != 0)
                .ok_or("Dual shape has a non-leading scalar component")?;
            let mut beta = alpha.clone();
            beta[axis] -= 1;
            let parent = &components[*indices
                .get(&beta)
                .ok_or("Dual shape is not ancestor closed")?];
            let mut derivative = Atom::Zero;
            for (index, gamma) in self.shape.iter().enumerate() {
                let mut successor = gamma.clone();
                successor[axis] += 1;
                let Some(&next) = indices.get(&successor) else {
                    continue;
                };
                for argument in 0..scalar_arguments.len() {
                    if zeros.contains(&(argument, next)) {
                        continue;
                    }
                    derivative += parent.derivative(&variables[argument * dim + index])
                        * &parameters[argument * dim + next]
                        * Atom::num(successor[axis]);
                }
            }
            components[component] = derivative / Atom::num(alpha[axis]);
            if components[component].contains_symbol(Symbol::DERIVATIVE) {
                return Err(format!(
                    "Cannot lower derivative {:?} of retained numerator {}",
                    alpha, family.definition.lhs
                ));
            }
        }
        let mut function_map = self.function_map.clone();
        let mut calls = Vec::with_capacity(dim);
        for (component, body) in components.into_iter().enumerate() {
            if body == Atom::Zero {
                calls.push(Atom::Zero);
                continue;
            }
            let component_symbol = symbol!(format!("{}_component_{component}", symbol.get_name()));
            function_map
                .add_tagged_function_with_options(
                    component_symbol,
                    vec![Atom::num(component)],
                    variables.clone(),
                    body,
                    FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
                )
                .map_err(|error| error.to_string())?;
            calls.push(
                FunctionBuilder::new(component_symbol)
                    .add_arg(component)
                    .add_args(&parameters)
                    .finish(),
            );
        }
        let evaluator = Atom::evaluator_multiple(&calls, &parameters)
            .function_map(function_map)
            .optimization_settings(self.optimization.clone())
            .build()
            .map_err(|error| error.to_string())?;
        self.cache
            .borrow_mut()
            .insert(request.clone(), evaluator.clone());
        Ok(Some(evaluator))
    }
}

pub(crate) fn build_dual_evaluator(
    expressions: &[Atom],
    parameters: &[Atom],
    function_map: &FunctionMap,
    definitions: &[RetainedFunctionDefinition],
    shape: Vec<Vec<usize>>,
    zero_components: Vec<(usize, usize)>,
    optimization: OptimizationSettings,
) -> Result<ExpressionEvaluator<Number>, String> {
    let scalar = Dualizer::new(
        HyperDual::<Number>::new(shape.clone()),
        zero_components.clone(),
    );
    if !definitions.iter().any(|definition| definition.retained) {
        return Atom::evaluator_multiple(expressions, parameters)
            .function_map(function_map.clone())
            .optimization_settings(optimization)
            .build()
            .map_err(|error| error.to_string())?
            .vectorize(&scalar);
    }
    if parameters.is_empty() {
        let expressions = expressions
            .iter()
            .flat_map(|expression| {
                std::iter::once(expression.clone()).chain((1..shape.len()).map(|_| Atom::Zero))
            })
            .collect::<Vec<_>>();
        return Atom::evaluator_multiple(&expressions, parameters)
            .function_map(function_map.clone())
            .optimization_settings(optimization)
            .build()
            .map_err(|error| error.to_string());
    }
    static NEXT_SCOPE: AtomicUsize = AtomicUsize::new(0);
    let scope = loop {
        let scope = NEXT_SCOPE.fetch_add(1, Ordering::Relaxed);
        let name = format!("gammalooprs::retained_dual_{scope}_0");
        if !symbolica::state::State::symbol_iter().any(|(symbol, existing)| {
            existing == name || symbol.get_aliases().iter().any(|alias| alias == &name)
        }) {
            break scope;
        }
    };
    let families = definitions
        .iter()
        .filter(|definition| definition.retained)
        .map(|definition| {
            let mut bound = parameters.to_vec();
            bound.extend(
                definition
                    .args
                    .iter()
                    .map(|argument| argument.as_view().to_owned()),
            );
            let body = resolve_definitions(&definition.rhs, definitions, &bound, false)?;
            let captures = parameters
                .iter()
                .filter(|parameter| {
                    !definition
                        .args
                        .iter()
                        .any(|arg| arg.as_view() == parameter.as_view())
                        && body.contains(*parameter)
                })
                .cloned()
                .collect();
            Ok(Family {
                definition: definition.clone(),
                body,
                captures,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let zeros = zero_components.iter().copied().collect::<HashSet<_>>();
    let mut symbols = HashMap::<Request, Symbol>::new();
    let mut requests = HashMap::<Symbol, Request>::new();
    let mut expressions_out = Vec::with_capacity(expressions.len());
    for expression in expressions {
        let expression = resolve_definitions(expression, definitions, parameters, true)?;
        let transformed = expression.replace_map_bottom_up(|view, _, out| {
            let AtomView::Fun(call) = view else {
                return;
            };
            let Some((family_index, family)) = families.iter().enumerate().find(|(_, family)| {
                family.definition.lhs.as_fun_view().is_some_and(|lhs| {
                    call.get_symbol() == lhs.get_symbol()
                        && call.get_nargs()
                            == family.definition.tags.len() + family.definition.args.len()
                        && call
                            .iter()
                            .zip(&family.definition.tags)
                            .all(|(a, b)| a == b.as_view())
                })
            }) else {
                return;
            };
            let mut arguments = call
                .iter()
                .skip(family.definition.tags.len())
                .map(|a| a.to_owned())
                .collect::<Vec<_>>();
            arguments.extend(family.captures.iter().cloned());
            arguments.push(parameters[0].clone());
            let mut local_zeros = Vec::new();
            for (argument_index, argument) in arguments.iter().enumerate() {
                // Close the supported monomials under multiplication. This is
                // conservative for every analytic expression of these inputs.
                let mut supported = vec![false; shape.len()];
                supported[0] = true;
                let mut depends_on = vec![false; parameters.len()];
                let mut pending = vec![argument.clone()];
                while let Some(value) = pending.pop() {
                    value.visitor(&mut |view| {
                        if let Some(index) =
                            parameters.iter().position(|atom| atom.as_view() == view)
                        {
                            depends_on[index] = true;
                            return false;
                        }
                        if let AtomView::Fun(call) = view
                            && let Some(request) = requests.get(&call.get_symbol())
                        {
                            let nested = &families[request.family];
                            let arguments = nested
                                .definition
                                .args
                                .iter()
                                .map(|arg| arg.as_view().to_owned())
                                .chain(nested.captures.iter().cloned());
                            // An inner closure's anchor and unused formals do not
                            // introduce derivatives into the enclosing argument.
                            for (formal, actual) in arguments.zip(call.iter()) {
                                if nested.body.contains(&formal) {
                                    pending.push(actual.to_owned());
                                }
                            }
                            return false;
                        }
                        true
                    });
                }
                for (parameter, depends_on) in depends_on.into_iter().enumerate() {
                    if depends_on {
                        for (component, support) in supported.iter_mut().enumerate().skip(1) {
                            *support |= !zeros.contains(&(parameter, component));
                        }
                    }
                }
                loop {
                    let mut next = supported.clone();
                    for left in shape
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| supported[*i])
                        .map(|(_, powers)| powers)
                    {
                        for right in shape
                            .iter()
                            .enumerate()
                            .filter(|(i, _)| supported[*i])
                            .map(|(_, powers)| powers)
                        {
                            let sum = left
                                .iter()
                                .zip(right)
                                .map(|(a, b)| a + b)
                                .collect::<Vec<_>>();
                            if let Some(component) = shape.iter().position(|powers| powers == &sum)
                            {
                                next[component] = true;
                            }
                        }
                    }
                    if next == supported {
                        break;
                    }
                    supported = next;
                }
                local_zeros.extend(supported.iter().enumerate().skip(1).filter_map(
                    |(component, &supported)| (!supported).then_some((argument_index, component)),
                ));
            }
            let request = Request {
                family: family_index,
                zeros: local_zeros,
            };
            let placeholder = *symbols.entry(request.clone()).or_insert_with(|| {
                let name = format!("gammalooprs::retained_dual_{scope}_{}", requests.len());
                let symbol = SymbolBuilder::new(symbolica::wrap_symbol!(name.as_str()))
                    .with_evaluation_info(EvaluationInfo::new())
                    .build()
                    .unwrap();
                requests.insert(symbol, request);
                symbol
            });
            let value = FunctionBuilder::new(placeholder)
                .add_args(&arguments)
                .finish();
            out.set_from_view(&value.as_view());
        });
        expressions_out.push(transformed);
    }
    let dualizer = SharedDualizer {
        scalar,
        shape: &shape,
        families: &families,
        requests: &requests,
        function_map,
        optimization: &optimization,
        scope,
        cache: RefCell::new(HashMap::new()),
    };
    Atom::evaluator_multiple(&expressions_out, parameters)
        .function_map(function_map.clone())
        .optimization_settings(optimization.clone())
        .build()
        .map_err(|error| error.to_string())?
        .vectorize(&dualizer)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bound_function_parameters_remain_opaque_when_resolving_aliases() {
        let alias = parse!("retained_dual_bound::alias");
        let parameter = function!(symbol!("retained_dual_bound::p"), &alias);
        let definitions = [RetainedFunctionDefinition {
            lhs: alias,
            rhs: parse!("retained_dual_bound::x"),
            tags: Vec::new(),
            args: Vec::new(),
            retained: false,
        }];
        assert_eq!(
            resolve_definitions(
                &parameter,
                &definitions,
                std::slice::from_ref(&parameter),
                true
            )
            .unwrap(),
            parameter,
        );
    }

    #[test]
    fn retained_mixed_higher_derivatives_match_scalar_dualizer() {
        let x = parse!("retained_dual_higher::x");
        let y = parse!("retained_dual_higher::y");
        let z = parse!("retained_dual_higher::z");
        let symbol = symbol!("retained_dual_higher::f");
        let body = z.clone().pow(3) * y.clone().sqrt();
        let definitions = [RetainedFunctionDefinition {
            lhs: function!(symbol, &z),
            rhs: body.clone(),
            tags: Vec::new(),
            args: vec![z.clone().try_into().unwrap()],
            retained: true,
        }];
        let parameters = [x.clone(), y];
        let shape = vec![
            vec![0, 0],
            vec![1, 0],
            vec![0, 1],
            vec![2, 0],
            vec![1, 1],
            vec![0, 2],
        ];
        let optimization = OptimizationSettings::new().cores(1);
        let mut retained = build_dual_evaluator(
            &[function!(symbol, &x)],
            &parameters,
            &FunctionMap::new(),
            &definitions,
            shape.clone(),
            Vec::new(),
            optimization.clone(),
        )
        .unwrap()
        .map_coeff(&|number| number.re.to_f64());
        let mut scalar = body
            .replace(z.to_pattern())
            .with(x)
            .evaluator(&parameters)
            .optimization_settings(optimization)
            .build()
            .unwrap()
            .vectorize(&Dualizer::new(HyperDual::<Number>::new(shape), Vec::new()))
            .unwrap()
            .map_coeff(&|number| number.re.to_f64());
        let input = [
            2.0, 1.0, 0.25, 0.125, -0.5, 0.75, 4.0, 0.5, -1.0, 0.25, 0.125, -0.25,
        ];
        let mut actual = [0.0; 6];
        let mut expected = [0.0; 6];
        retained.evaluate(&input, &mut actual);
        scalar.evaluate(&input, &mut expected);
        for (actual, expected) in actual.into_iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn retained_dual_closes_nested_bodies_and_prunes_zero_mass_seeds() {
        let x = parse!("retained_dual_test::x");
        let y = parse!("retained_dual_test::y");
        let mass = parse!("retained_dual_test::mass");
        let z = parse!("retained_dual_test::z");
        let inner = symbol!("retained_dual_test::inner");
        let outer = symbol!("retained_dual_test::outer");
        let definitions = vec![
            RetainedFunctionDefinition {
                lhs: function!(inner, &z),
                rhs: &z * &z + mass.clone().sqrt(),
                tags: vec![],
                args: vec![z.clone().try_into().unwrap()],
                retained: true,
            },
            RetainedFunctionDefinition {
                lhs: function!(outer, &z),
                rhs: function!(inner, &z) + &y,
                tags: vec![],
                args: vec![z.clone().try_into().unwrap()],
                retained: true,
            },
        ];
        let mut function_map = FunctionMap::new();
        for definition in &definitions {
            function_map
                .add_function_with_options(
                    definition.lhs.as_fun_view().unwrap().get_symbol(),
                    definition.args.clone(),
                    definition.rhs.clone(),
                    FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
                )
                .unwrap();
        }
        let expression = function!(outer, &x + Atom::num(1)) + function!(outer, &x * Atom::num(2));
        let nested = function!(outer, function!(inner, &x));
        let guarded = function!(Symbol::IF, &x, x.clone().pow(-1), &nested);
        let evaluator = build_dual_evaluator(
            &[expression, nested, guarded],
            &[x, y, mass],
            &function_map,
            &definitions,
            vec![vec![0], vec![1]],
            vec![(1, 1), (2, 1)],
            OptimizationSettings::new().cores(1),
        )
        .unwrap();
        let exported = evaluator.export_instructions();
        assert_eq!(exported.sub_evaluators.len(), 4);
        let encoded = bincode::encode_to_vec(&evaluator, bincode::config::standard()).unwrap();
        let (evaluator, _): (ExpressionEvaluator<Number>, _) =
            bincode::decode_from_slice(&encoded, bincode::config::standard()).unwrap();
        let mut eager = evaluator.map_coeff(&|number| number.re.to_f64());
        let cases = [
            (
                [3.0, 1.0, 5.0, 0.0, 0.0, 0.0],
                [62.0, 32.0, 86.0, 108.0, 1.0 / 3.0, -1.0 / 9.0],
            ),
            (
                [0.0, 1.0, 5.0, 0.0, 0.0, 0.0],
                [11.0, 2.0, 5.0, 0.0, 5.0, 0.0],
            ),
        ];
        for (input, expected) in cases {
            let mut output = [0.0; 6];
            eager.evaluate(&input, &mut output);
            assert_eq!(output, expected);
        }
        let mut jit = eager
            .jit_compile(symbolica::evaluate::JITCompilationSettings::new().optimization_level(0))
            .unwrap();
        for (input, expected) in cases {
            let mut output = [0.0; 6];
            jit.evaluate(&input, &mut output);
            assert_eq!(output, expected);
        }
        let directory =
            std::env::temp_dir().join(format!("ltd_retained_dual_{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        for assembly in [false, true] {
            let stem = if assembly {
                "dual_assembly"
            } else {
                "dual_cpp"
            };
            let cpp = directory.join(format!("{stem}.cpp"));
            let library = directory.join(format!("{stem}.so"));
            let mut compiled = eager
                .export_cpp::<f64>(
                    &cpp,
                    stem,
                    symbolica::evaluate::ExportSettings::new().inline_asm(if assembly {
                        symbolica::evaluate::InlineASM::X64
                    } else {
                        symbolica::evaluate::InlineASM::None
                    }),
                )
                .unwrap()
                .compile(&library, symbolica::evaluate::CompileOptions::default())
                .unwrap()
                .load()
                .unwrap();
            for (input, expected) in cases {
                let mut output = [0.0; 6];
                compiled.evaluate(&input, &mut output);
                assert_eq!(output, expected);
            }
        }
    }
}
