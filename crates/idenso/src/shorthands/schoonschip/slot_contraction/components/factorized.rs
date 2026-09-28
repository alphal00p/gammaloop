//! State merging on the collector's existing port graph.
//!
//! Remaining factors are opaque terminals of the same component reducer used
//! for a completed monomial. Only the selected factor contributes alternatives;
//! substitutions remain port bindings until a distinct output variable is emitted.

use super::*;

// A state contains occurrence-local port bindings, not rebuilt expressions.
// Variables and the incidence scratch still belong to ComponentSum.
#[derive(Clone, PartialEq, Eq, Hash)]
struct Remaining<'a> {
    factors: Vec<(usize, Vec<Argument<'a>>)>,
    residual: Vec<(usize, u16)>,
}

type FactorTerm = (Vec<(usize, u16)>, Rational);

/// Exact output may retain pending contractions when a frontier budget is met.
/// Algebraic validity and completion are separate facts.
pub(crate) struct FactorizedContraction {
    pub(crate) root: Atom,
    pub(crate) aliases: Vec<(Atom, Atom)>,
    pub(crate) complete: bool,
}

impl SlotContraction {
    pub(crate) fn contract_factorized(
        &self,
        source: AtomView<'_>,
        order: Option<&[usize]>,
    ) -> Option<FactorizedContraction> {
        if source.needs_normalization()
            || self.metric.get_evaluation_info().is_some()
            || !InterfaceInference::normalization_is_intrinsic(source)
        {
            return None;
        }
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(self, Intake::FactoredContraction, &mut slots);
        sum.preserve_scalar_factors = true;
        let factors = match source {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            _ => vec![source],
        };
        let mut selected = Vec::new();
        let mut remaining = Vec::new();
        for (position, &factor) in factors.iter().enumerate() {
            let structure = InterfaceInference::replacement_interface(factor).ok()?;
            if !structure.open_positions().is_empty() {
                return None;
            }
            let mut ports = Vec::new();
            for slot in structure.logical_slots() {
                let wanted = crate::tensor::composition::port_atom(slot);
                let mut found = None;
                factors[position].visitor(&mut |value| {
                    if value == wanted.as_view() {
                        found = Some(value);
                    }
                    found.is_none()
                });
                ports.push(found?);
            }
            if sum.owns_factor(factor) {
                selected.push(position);
            }
            remaining.push((
                position,
                ports.iter().copied().map(Argument::Original).collect(),
            ));
            sum.opaque_factors.push((factor, ports, structure));
        }
        if let Some(order) = order {
            let mut seen = vec![false; sum.opaque_factors.len()];
            for &position in order {
                if *seen.get(position)? {
                    return None;
                }
                seen[position] = true;
            }
            if seen.iter().any(|&entry| !entry) {
                return None;
            }
            selected
                .sort_by_key(|position| order.iter().position(|entry| entry == position).unwrap());
        }
        // Work on connected factors independently. Scalar spectators and
        // disconnected foreign sums must not enter another component's states.
        let mut incidence: AHashMap<_, Vec<usize>> = AHashMap::new();
        for (position, (_, ports, _)) in sum.opaque_factors.iter().enumerate() {
            for &port in ports {
                incidence.entry(port).or_default().push(position);
            }
        }
        let mut seen = vec![false; remaining.len()];
        let mut roots = Vec::new();
        let mut definitions = Vec::new();
        let mut complete = true;
        for start in 0..remaining.len() {
            if seen[start] {
                continue;
            }
            seen[start] = true;
            let mut component = vec![start];
            let mut position = 0;
            while position < component.len() {
                for port in &sum.opaque_factors[component[position]].1 {
                    for &neighbor in &incidence[port] {
                        if !seen[neighbor] {
                            seen[neighbor] = true;
                            component.push(neighbor);
                        }
                    }
                }
                position += 1;
            }
            let selected = selected
                .iter()
                .copied()
                .filter(|position| component.contains(position))
                .collect::<Vec<_>>();
            if selected.is_empty() {
                roots.push(Atom::mul_many(
                    component
                        .iter()
                        .map(|&position| sum.opaque_factors[position].0),
                ));
                continue;
            }
            let result = sum.contract_states(
                Remaining {
                    factors: component
                        .into_iter()
                        .map(|position| remaining[position].clone())
                        .collect(),
                    residual: Vec::new(),
                },
                &selected,
            )?;
            roots.push(result.root);
            definitions.extend(result.aliases);
            complete &= result.complete;
        }
        Some(FactorizedContraction {
            root: Atom::mul_many(roots),
            aliases: definitions,
            complete,
        })
    }
}

impl<'a> ComponentSum<'a, '_> {
    fn owns_factor(&self, source: AtomView<'_>) -> bool {
        match source {
            AtomView::Add(sum) => sum.iter().any(|term| self.owns_factor(term)),
            AtomView::Mul(product) => product.iter().any(|factor| self.owns_factor(factor)),
            AtomView::Pow(power) => self.owns_factor(power.get_base_exp().0),
            AtomView::Fun(function) => {
                let head = function.get_symbol();
                head == self.contractor.metric || head.has_tag(&self.contractor.tags.rank1)
            }
            _ => false,
        }
    }

    fn opaque_factor(&mut self, position: usize, arguments: &[Argument<'a>]) -> Option<()> {
        let tensor = self.tensors.len();
        self.tensors
            .push((TensorSource::Factor(position), arguments.to_vec()));
        for (port, &argument) in arguments.iter().enumerate() {
            if let Argument::Original(slot) = argument {
                let (space, index) = self.resolve_endpoint(slot)?;
                let node = self.endpoint(slot, space, index)?;
                self.nodes[node]
                    .terminals
                    .push(Terminal::Tensor(tensor, port));
            }
        }
        Some(())
    }

    fn residual_factor(&mut self, variable: usize, exponent: u16) -> Option<()> {
        match self.variables[variable].clone() {
            Variable::Vector(head, slot) => {
                let (space, index) = self.resolve_endpoint(slot)?;
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            Variable::Metric([first, second]) => {
                self.metric(
                    Argument::Original(first),
                    Argument::Original(second),
                    exponent,
                )?;
            }
            Variable::Tensor(source, arguments) => {
                if exponent != 1 {
                    return None;
                }
                let tensor = self.tensors.len();
                for (position, &argument) in arguments.iter().enumerate() {
                    if let Argument::Original(slot) = argument
                        && matches!(self.slots.classify(slot), SlotMatch::Explicit(_))
                    {
                        let (space, index) = self.resolve_endpoint(slot)?;
                        let node = self.endpoint(slot, space, index)?;
                        self.nodes[node]
                            .terminals
                            .push(Terminal::Tensor(tensor, position));
                    }
                }
                self.tensors.push((source, arguments));
            }
            value => self.variable(value, exponent),
        }
        Some(())
    }

    fn reset_components(&mut self, coefficient: &Rational) {
        self.nodes.clear();
        self.occurrences.clear();
        self.factors.clear();
        self.tensors.clear();
        self.metrics = 0;
        self.coefficient = coefficient.clone();
    }

    fn factor_terms(&mut self, position: usize) -> Option<Vec<FactorTerm>> {
        let source = self.opaque_factors[position].0;
        let (node, _) = self.compile_input(source, 0)?;
        self.input_terms.clear();
        self.distribute(&mut vec![(node, 1)], &Rational::one(), 0)?;
        let mut terms = std::mem::take(&mut self.input_terms)
            .into_iter()
            .filter(|(_, (_, coefficient))| !coefficient.is_zero())
            .collect::<Vec<_>>();
        terms.sort_unstable_by_key(|(_, (position, _))| *position);
        Some(
            terms
                .into_iter()
                .map(|(factors, (_, coefficient))| (factors, coefficient))
                .collect(),
        )
    }

    fn contract_states(
        &mut self,
        remaining: Remaining<'a>,
        order: &[usize],
    ) -> Option<FactorizedContraction> {
        let mut states = vec![(remaining, Atom::num(1))];
        let mut definitions = Vec::new();
        let mut emitted: AHashMap<usize, Atom> = AHashMap::new();
        let mut weights = AHashMap::<Atom, Atom>::new();
        let mut generated = 0usize;
        let mut definition_bytes = 0usize;
        let mut complete = true;
        for &selected in order {
            let terms = self.factor_terms(selected)?;
            // Predict growth of the whole frontier, not just the local
            // template's term iterator. This includes graph-state storage and
            // existing definitions, not a hard heap bound on Atom buffers.
            // Flush exact remaining factors before an exponential next level;
            // explicit materialization may expand them.
            let state_bytes = states
                .iter()
                .map(|(remaining, _)| {
                    std::mem::size_of::<Remaining<'a>>()
                        + remaining
                            .factors
                            .iter()
                            .map(|(_, arguments)| {
                                std::mem::size_of::<(usize, Vec<Argument<'a>>)>()
                                    + arguments.len() * std::mem::size_of::<Argument<'a>>()
                            })
                            .sum::<usize>()
                        + remaining.residual.len() * std::mem::size_of::<(usize, u16)>()
                })
                .max()
                .unwrap_or(0);
            let count = states.len().saturating_mul(terms.len());
            let predicted = count.saturating_mul(state_bytes.saturating_mul(4).saturating_add(256));
            const MAX_FRONTIER_BYTES: usize = if cfg!(test) {
                64 * 1024
            } else {
                64 * 1024 * 1024
            };
            if generated.saturating_add(count) > Self::MAX_GENERATED_TERMS
                || predicted.saturating_add(definition_bytes) > MAX_FRONTIER_BYTES
            {
                complete = false;
                break;
            }
            generated += count;
            let mut positions: AHashMap<Remaining<'a>, usize> = AHashMap::new();
            let mut next: Vec<(Remaining<'a>, Vec<Atom>)> = Vec::new();
            for (remaining, weight) in states {
                let source = remaining
                    .factors
                    .iter()
                    .find(|(position, _)| *position == selected)?;
                self.overrides = self.opaque_factors[selected]
                    .1
                    .iter()
                    .copied()
                    .zip(source.1.iter().copied())
                    .collect();
                for (factors, coefficient) in &terms {
                    self.reset_components(coefficient);
                    for &(factor, exponent) in factors {
                        let value = self.input_atoms[factor];
                        if self.scalar_spectator(value) {
                            self.variable(Variable::Scalar(value), exponent);
                        } else {
                            self.factor(value, exponent)?;
                        }
                    }
                    for &(variable, exponent) in &remaining.residual {
                        self.residual_factor(variable, exponent)?;
                    }
                    for (position, arguments) in &remaining.factors {
                        if *position != selected {
                            self.opaque_factor(*position, arguments)?;
                        }
                    }
                    self.reduce_components()?;
                    let mut key = Remaining {
                        factors: Vec::new(),
                        residual: Vec::new(),
                    };
                    let mut coefficient = vec![Atom::num(self.coefficient.clone()), weight.clone()];
                    for &(variable, exponent) in &self.factors {
                        match &self.variables[variable] {
                            Variable::Tensor(TensorSource::Factor(position), arguments) => {
                                key.factors.push((*position, arguments.clone()));
                            }
                            value @ (Variable::Scalar(_) | Variable::Dot(_, _)) => {
                                let atom = if let Some(value) = emitted.get(&variable) {
                                    value.clone()
                                } else {
                                    let atom = self.emit_variable(value)?;
                                    emitted.insert(variable, atom.clone());
                                    atom
                                };
                                coefficient.push(atom.pow(exponent));
                            }
                            _ => key.residual.push((variable, exponent)),
                        }
                    }
                    key.factors.sort_unstable_by_key(|(position, _)| *position);
                    key.residual.sort_unstable();
                    let coefficient = Atom::mul_many(coefficient);
                    if coefficient.is_zero() {
                        continue;
                    }
                    if let Some(&position) = positions.get(&key) {
                        next[position].1.push(coefficient);
                    } else {
                        positions.insert(key.clone(), next.len());
                        next.push((key, vec![coefficient]));
                    }
                }
            }
            self.overrides.clear();
            states = Vec::with_capacity(next.len());
            for (key, coefficients) in next {
                let body = Atom::add_many(coefficients);
                if body.is_zero() {
                    continue;
                }
                // Different remaining states often carry the same coefficient.
                // Keep one literal definition for that value, and avoid aliases
                // whose body is just a number or an existing weight handle.
                if matches!(body.as_view(), AtomView::Num(_)) {
                    states.push((key, body));
                    continue;
                }
                if let Some(handle) = weights.get(&body) {
                    states.push((key, handle.clone()));
                    continue;
                }
                let body =
                    SymbolicTensor::checked_parts(body, PartialStructure::from_logical_slots([]))
                        .ok()?;
                let handle = body.alias_handle().ok()?;
                definition_bytes = definition_bytes
                    .saturating_add(body.expression.as_view().get_byte_size())
                    .saturating_add(handle.expression.as_view().get_byte_size());
                weights.insert(body.expression.clone(), handle.expression.clone());
                weights.insert(handle.expression.clone(), handle.expression.clone());
                definitions.push((handle.expression.clone(), body.expression));
                states.push((key, handle.expression));
            }
        }
        let mut terms = Vec::with_capacity(states.len());
        for (remaining, weight) in states {
            let mut factors = vec![weight];
            for (position, arguments) in remaining.factors {
                factors.push(self.emit_factor(position, &arguments)?);
            }
            for (variable, exponent) in remaining.residual {
                let atom = if let Some(value) = emitted.get(&variable) {
                    value.clone()
                } else {
                    let atom = self.emit_variable(&self.variables[variable])?;
                    emitted.insert(variable, atom.clone());
                    atom
                };
                factors.push(atom.pow(exponent));
            }
            terms.push(Atom::mul_many(factors));
        }
        Some(FactorizedContraction {
            root: Atom::add_many(terms),
            aliases: definitions,
            complete,
        })
    }

    pub(super) fn argument(&self, value: AtomView<'a>) -> Argument<'a> {
        self.overrides
            .iter()
            .find_map(|&(source, target)| (source == value).then_some(target))
            .unwrap_or(Argument::Original(value))
    }

    fn emit_argument(&self, argument: Argument<'a>) -> Atom {
        match argument {
            Argument::Original(value) => value.to_owned(),
            Argument::Vector(space, head) => {
                let (representation, dimension) = self.spaces[space];
                FunctionBuilder::new(head)
                    .add_arg(representation.to_symbolic([dimension]))
                    .finish()
            }
        }
    }

    pub(super) fn emit_factor(&self, position: usize, arguments: &[Argument<'a>]) -> Option<Atom> {
        let (source, ports, interface) = &self.opaque_factors[position];
        let replacements = ports
            .iter()
            .zip(arguments)
            .enumerate()
            .filter(|(_, (port, argument))| **argument != Argument::Original(**port))
            .map(|(position, (_, &argument))| (position, self.emit_argument(argument)))
            .collect::<std::collections::HashMap<_, _>>();
        if replacements.is_empty() {
            return Some(source.to_owned());
        }
        let source = SymbolicTensor::from_normalized_parts(source.to_owned(), interface.clone());
        let expression = crate::tensor::composition::rewrite_interface_ports(
            &source,
            &replacements,
            &mut |_, _| {},
        )
        .ok()?;
        Some(crate::shorthands::schoonschip::DotNormalizer::run(
            expression.as_view(),
        ))
    }

    pub(super) fn metric(
        &mut self,
        first: Argument<'a>,
        second: Argument<'a>,
        exponent: u16,
    ) -> Option<()> {
        match (first, second) {
            (Argument::Vector(space, a), Argument::Vector(other, b)) => {
                if space != other {
                    return None;
                }
                self.dot(space, a, b, exponent);
            }
            (Argument::Vector(space, head), Argument::Original(slot))
            | (Argument::Original(slot), Argument::Vector(space, head)) => {
                let (other, index) = self.resolve_endpoint(slot)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            (Argument::Original(first), Argument::Original(second)) => {
                let (space, first_index) = self.resolve_endpoint(first)?;
                let (other, second_index) = self.resolve_endpoint(second)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    if first_index == second_index {
                        return None;
                    }
                    self.factor(self.spaces[space].1, exponent / 2)?;
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let a = self.endpoint(first, space, first_index)?;
                    let b = self.endpoint(second, space, second_index)?;
                    let a = self.root(a);
                    let b = self.root(b);
                    self.metrics += 1;
                    self.nodes[a].metrics += 1;
                    self.nodes[b].parent = a;
                }
            }
        }
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::{input, setup};
    use super::*;
    use crate::shorthands::schoonschip::{Schoonschip, SchoonschipSettings};
    use ahash::AHashSet;
    use symbolica::atom::AliasedAtom;

    fn resolve(contracted: FactorizedContraction) -> Atom {
        let mut result = AliasedAtom::from(contracted.root);
        for (handle, body) in contracted.aliases {
            result.register_alias(handle, body);
        }
        result.into_inner()
    }

    #[test]
    fn factorized_states_use_the_component_reducer_in_both_orders() {
        let contractor = setup();
        for source in [
            "g(mink(4,a),mink(4,b))*(p(mink(4,a))+q(mink(4,a)))*r(mink(4,b))",
            "(g(mink(4,a),mink(4,b))*p(mink(4,c))+g(mink(4,a),mink(4,c))*q(mink(4,b)))*(p(mink(4,a))*q(mink(4,b))*r(mink(4,c))+q(mink(4,a))*r(mink(4,b))*s(mink(4,c)))",
            "(g(mink(4,a),mink(4,b))+p(mink(4,a))*q(mink(4,b)))*(g(mink(4,a),mink(4,b))+r(mink(4,a))*s(mink(4,b)))",
            "(g(mink(4,a),mink(4,b))+t(mink(4,a),mink(4,b)))*p(mink(4,a))*q(mink(4,b))",
        ] {
            let source = input(source);
            let expected = source
                .schoonschip_with_settings(
                    &SchoonschipSettings::default().with_expanded_contracted_sums(),
                )
                .expand();
            let count = match source.as_view() {
                AtomView::Mul(product) => product.iter().len(),
                _ => 1,
            };
            let reversed = (0..count).rev().collect::<Vec<_>>();
            for order in [None, Some(reversed.as_slice())] {
                let contracted = contractor
                    .contract_factorized(source.as_view(), order)
                    .unwrap();
                let mut bodies = AHashSet::new();
                let handles = contracted
                    .aliases
                    .iter()
                    .map(|(handle, _)| handle)
                    .collect::<AHashSet<_>>();
                for (_, body) in &contracted.aliases {
                    assert!(bodies.insert(body));
                    assert!(!matches!(body.as_view(), AtomView::Num(_)));
                    assert!(!handles.contains(body));
                }
                assert_eq!(resolve(contracted).expand(), expected, "{source}");
            }
        }
    }

    #[test]
    fn factorized_states_keep_foreign_and_disconnected_sums() {
        let contractor = setup();
        let spectator = input("(1+x)*(t(mink(4,c))+routing(z)*u(mink(4,c)))");
        let first = input("(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))");
        let second = input("(p(mink(4,b))+s(mink(4,b)))*q(mink(4,b))");
        let source = Atom::mul_many([&spectator, &first, &second]);
        let contracted = contractor
            .contract_factorized(source.as_view(), None)
            .unwrap();
        let result = resolve(contracted);
        let expected = Atom::mul_many([
            spectator,
            first.schoonschip_with_settings(
                &SchoonschipSettings::default().with_expanded_contracted_sums(),
            ),
            second.schoonschip_with_settings(
                &SchoonschipSettings::default().with_expanded_contracted_sums(),
            ),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn factorized_port_bindings_remain_occurrence_local() {
        let contractor = setup();
        let source =
            input("(p(mink(4,a))+q(mink(4,a)))*(p(mink(4,b))+q(mink(4,b)))*t(mink(4,a),mink(4,b))");
        let contracted = contractor
            .contract_factorized(source.as_view(), None)
            .unwrap();
        assert_eq!(
            resolve(contracted).expand(),
            source
                .schoonschip_with_settings(
                    &SchoonschipSettings::default().with_expanded_contracted_sums(),
                )
                .expand()
        );
        assert!(
            contractor
                .contract_factorized(source.as_view(), Some(&[0, 0, 1]))
                .is_none()
        );
    }

    #[test]
    fn factorized_frontier_budget_retains_an_exact_unexpanded_remainder() {
        let contractor = setup();
        let slots = (0..9).map(|i| format!("mink(4,a{i})")).collect::<Vec<_>>();
        let source = input(&format!(
            "t({})*{}",
            slots.join(","),
            slots
                .iter()
                .map(|slot| format!("(x*p({slot})+y*q({slot}))"))
                .collect::<Vec<_>>()
                .join("*")
        ));
        let contracted = contractor
            .contract_factorized(source.as_view(), None)
            .unwrap();
        assert!(!contracted.aliases.is_empty());
        assert!(!contracted.complete);
        use spenso::network::parsing::AtomStructureExt;
        assert!(
            contracted.root.has_repeated_explicit_indices(),
            "the bounded frontier must keep its unprocessed contractions"
        );
        let result = resolve(contracted);
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        assert_eq!(
            result
                .expand()
                .schoonschip_with_settings(&settings)
                .expand(),
            source
                .expand()
                .schoonschip_with_settings(&settings)
                .expand()
        );
    }
}
