//! Distributional integration by parts followed by Fermi-shell localization.
//!
//! Normal flows lower delta derivatives while preserving the other active
//! shell energies. Soft momentum charts prevent differentiating an integrable
//! pole into an unaccounted distribution. The complete coefficient, including
//! grouped thermal weights, is differentiated; its thermal derivatives generate
//! intersection contacts. Final ordinary deltas are localized together.

use std::collections::{BTreeMap, BTreeSet, btree_map::Entry};

use eyre::{Result, eyre};
use linnet::half_edge::involution::EdgeIndex;
use spenso::structure::concrete_index::ExpandedIndex;
use symbolica::{
    atom::{AliasedAtom, Atom, AtomCore, AtomView, Indeterminate, Symbol},
    function, symbol,
};

use crate::{
    graph::{Graph, LoopMomentumBasis, lmb::LMBwithEdges},
    integrands::process::param_builder::{ParamBuilder, ParamBuilderGraph},
    momentum::{SignOrZero, sample::LoopIndex},
    utils::GS,
    uv::uv_graph::UVE,
};

#[path = "fermi_surface_domain.rs"]
mod domain;
#[path = "fermi_surface_flow.rs"]
mod flow;

pub(crate) use domain::FermiSurfaceDomain;

/// `order` is the thermal derivative order: one denotes an ordinary delta.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct FermiDistribution {
    edge: EdgeIndex,
    orientation: Atom,
    order: usize,
}

type Sectors = BTreeMap<Vec<FermiDistribution>, Atom>;

struct FermiNormalFlow {
    velocity: Vec<[Atom; 3]>,
    /// V(E²)/E², with the soft factor removed within each momentum chart.
    energy_square_rates: BTreeMap<EdgeIndex, Atom>,
    domains: Vec<FermiSurfaceDomain>,
}

pub(crate) struct FermiSurfaceLocalizer<'a> {
    graph: &'a Graph,
}

impl<'a> FermiSurfaceLocalizer<'a> {
    pub(crate) fn new(graph: &'a Graph) -> Self {
        Self { graph }
    }

    /// This runs after tensor contraction, while thermal functions still carry
    /// their edge, derivative and orientation metadata. No numerator polynomial
    /// or difference of thermal products is expanded into integration channels.
    pub(crate) fn localize(
        &self,
        atoms: Vec<AliasedAtom>,
        builder: &ParamBuilder,
    ) -> Result<(Vec<AliasedAtom>, Vec<EdgeIndex>, Vec<FermiSurfaceDomain>)> {
        let hidden_distributions = builder
            .reps
            .iter()
            .any(|entry| Self::contains_delta(entry.rhs.as_view()));
        let mut localized_edges = BTreeSet::new();
        let mut domains = BTreeSet::new();
        let atoms = atoms
            .into_iter()
            .map(|atom| {
                if !hidden_distributions
                    && !Self::contains_delta(atom.get_root().as_view())
                    && !atom
                        .get_aliases()
                        .values()
                        .any(|body| Self::contains_delta(body.as_view()))
                {
                    return Ok(atom);
                }
                // Scalar aliases and function bodies can hide a thermal weight.
                // Resolve dependencies without distributing any products.
                let root = self.expand_functions(atom.clone().into_inner(), builder, false)?;
                if !Self::contains_delta(root.as_view()) {
                    return Ok(atom);
                }
                let sectors = self.extract(root.as_view())?;
                let (reduced, constraints) = self.reduce_derivatives(sectors, builder)?;
                domains.extend(constraints);
                let localized = reduced
                    .into_iter()
                    .map(|(support, coefficient)| {
                        if support.is_empty() || coefficient.is_zero() {
                            Ok(coefficient)
                        } else {
                            localized_edges.extend(support.iter().map(|factor| factor.edge));
                            self.localize_support(&support, coefficient, builder)
                        }
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(AliasedAtom::from(Atom::add_many(localized)))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok((
            atoms,
            localized_edges.into_iter().collect(),
            domains.into_iter().collect(),
        ))
    }

    fn distribution(atom: AtomView<'_>) -> Result<Option<FermiDistribution>> {
        let AtomView::Fun(call) = atom else {
            return Ok(None);
        };
        if call.get_symbol() != GS.thermal_distribution {
            return Ok(None);
        }
        let args = call.iter().collect::<Vec<_>>();
        if args.len() != 5 {
            return Err(eyre!("thermal distribution must have five arguments"));
        }
        if args[2] != Atom::num(0).as_view() {
            return Ok(None);
        }
        Ok(Some(FermiDistribution {
            edge: EdgeIndex(
                usize::try_from(args[0])
                    .map_err(|_| eyre!("Fermi distribution has a noninteger edge: {}", args[0]))?,
            ),
            order: usize::try_from(args[1]).map_err(|_| {
                eyre!(
                    "Fermi distribution has a noninteger derivative order: {}",
                    args[1]
                )
            })?,
            orientation: args[4].to_owned(),
        }))
    }

    fn contains_delta(atom: AtomView<'_>) -> bool {
        let mut found = false;
        atom.visitor(&mut |part| {
            found |= Self::distribution(part)
                .ok()
                .flatten()
                .is_some_and(|factor| factor.order > 0);
            !found
        });
        found
    }

    fn add_sector(sectors: &mut Sectors, mut support: Vec<FermiDistribution>, coefficient: Atom) {
        if coefficient.is_zero() {
            return;
        }
        support.sort();
        match sectors.entry(support) {
            Entry::Vacant(entry) => {
                entry.insert(coefficient);
            }
            Entry::Occupied(mut entry) => {
                *entry.get_mut() += coefficient;
                if entry.get().is_zero() {
                    let _ = entry.remove();
                }
            }
        }
    }

    /// Expand only the sum/product structure needed to identify explicit delta
    /// supports. Ordinary thermal numerators remain part of the coefficient.
    fn extract(&self, atom: AtomView<'_>) -> Result<Sectors> {
        if let Some(factor) = Self::distribution(atom)?
            && factor.order > 0
        {
            if self.graph[factor.edge].chemical_potential_atom().is_none() {
                return Ok(BTreeMap::new());
            }
            if !self.graph[factor.edge].is_fermion() {
                return Err(eyre!(
                    "Fermi-surface localization requires a fermion, edge {}",
                    factor.edge
                ));
            }
            return Ok(BTreeMap::from([(vec![factor], Atom::one())]));
        }
        if !Self::contains_delta(atom) {
            return Ok(BTreeMap::from([(Vec::new(), atom.to_owned())]));
        }
        match atom {
            AtomView::Add(sum) => {
                let mut sectors = BTreeMap::new();
                for term in sum.iter() {
                    for (support, coefficient) in self.extract(term)? {
                        Self::add_sector(&mut sectors, support, coefficient);
                    }
                }
                Ok(sectors)
            }
            AtomView::Mul(product) => {
                let mut sectors = BTreeMap::from([(Vec::new(), Atom::one())]);
                for term in product.iter() {
                    let factors = self.extract(term)?;
                    let mut next = BTreeMap::new();
                    for (left, a) in &sectors {
                        for (right, b) in &factors {
                            let mut support = left.clone();
                            support.extend(right.iter().cloned());
                            Self::add_sector(&mut next, support, a * b);
                        }
                    }
                    sectors = next;
                }
                Ok(sectors)
            }
            AtomView::Fun(call)
                if call.get_symbol() == GS.thermal_weight_wrapper
                    || call.get_symbol() == GS.tree_denom_wrapper =>
            {
                self.extract(call.iter().next().unwrap())
            }
            AtomView::Fun(call) if call.get_symbol() == Symbol::IF => {
                let args = call.iter().collect::<Vec<_>>();
                if args.len() != 3 || Self::contains_delta(args[0]) {
                    return Err(eyre!("unsupported conditional Fermi distribution"));
                }
                let mut sectors = BTreeMap::new();
                for (index, body) in args[1..].iter().enumerate() {
                    for (support, coefficient) in self.extract(*body)? {
                        let branches = if index == 0 {
                            [coefficient, Atom::Zero]
                        } else {
                            [Atom::Zero, coefficient]
                        };
                        Self::add_sector(
                            &mut sectors,
                            support,
                            Symbol::IF.call_args([
                                args[0].to_owned(),
                                branches[0].clone(),
                                branches[1].clone(),
                            ]),
                        );
                    }
                }
                Ok(sectors)
            }
            _ => Err(eyre!(
                "Fermi distribution occurs outside a linear distribution product: {atom}"
            )),
        }
    }

    fn expand_functions(
        &self,
        mut atom: Atom,
        builder: &ParamBuilder,
        expand_occupations: bool,
    ) -> Result<Atom> {
        let replacements = builder
            .reps
            .iter()
            .filter(|entry| {
                expand_occupations
                    || entry
                        .lhs
                        .as_fun_view()
                        .is_none_or(|call| call.get_symbol() != GS.thermal_distribution)
            })
            .map(|entry| entry.replacement())
            .collect::<Vec<_>>();
        // Every change must advance through an acyclic function dependency.
        // The bound also catches accidental recursive parameter definitions.
        for _ in 0..=replacements.len() {
            let next = atom.replace_multiple(&replacements);
            if next == atom {
                return Ok(atom);
            }
            atom = next;
        }
        Err(eyre!(
            "cyclic scalar functions in Fermi-surface localization"
        ))
    }

    fn spatial(edge: EdgeIndex, axis: usize) -> Atom {
        GS.emr_mom(edge, Atom::from(ExpandedIndex::from_iter([axis + 1])))
    }

    fn momentum(&self, edge: EdgeIndex, builder: &ParamBuilder) -> Result<[Atom; 3]> {
        Ok([
            self.expand_functions(Self::spatial(edge, 0), builder, false)?,
            self.expand_functions(Self::spatial(edge, 1), builder, false)?,
            self.expand_functions(Self::spatial(edge, 2), builder, false)?,
        ])
    }

    fn energy(&self, edge: EdgeIndex, builder: &ParamBuilder) -> Result<Atom> {
        self.expand_functions(self.graph.explicit_ose_atom(edge), builder, false)
    }

    fn norm_squared(momentum: &[Atom; 3]) -> Atom {
        Atom::add_many(momentum.iter().map(|component| component.pow(2)))
    }

    fn basis(&self, support: &[FermiDistribution]) -> Result<LoopMomentumBasis> {
        // Grouped thermal derivatives leave a tree connector unselected at
        // each CFF contraction. Their active edges therefore extend a cotree;
        // check that invariant again after generating intersection contacts.
        let edges = support.iter().map(|factor| factor.edge).collect::<Vec<_>>();
        if edges.iter().copied().collect::<BTreeSet<_>>().len() != edges.len() {
            return Err(eyre!(
                "coincident Fermi factors must be reduced before localization: {support:?}"
            ));
        }
        let basis = self.graph.lmb_with_loop_edges(edges.as_slice())?;
        if edges.iter().any(|edge| !basis.loop_edges.contains(edge)) {
            return Err(eyre!(
                "Fermi supports do not extend an independent loop basis: {support:?}"
            ));
        }
        Ok(basis)
    }

    fn sign_atom(sign: SignOrZero) -> Atom {
        Atom::num(match sign {
            SignOrZero::Plus => 1,
            SignOrZero::Minus => -1,
            SignOrZero::Zero => 0,
        })
    }

    fn reduce_derivatives(
        &self,
        sectors: Sectors,
        builder: &ParamBuilder,
    ) -> Result<(Sectors, Vec<FermiSurfaceDomain>)> {
        let mut pending = sectors;
        let mut result = BTreeMap::new();
        let mut domains = BTreeSet::new();
        while let Some((mut support, coefficient)) = pending.pop_first() {
            if coefficient.is_zero() {
                continue;
            }
            let Some(index) = support.iter().position(|factor| factor.order > 1) else {
                Self::add_sector(&mut result, support, coefficient);
                continue;
            };
            let target = support[index].edge;
            self.basis(&support)?;
            let FermiNormalFlow {
                velocity: flow,
                energy_square_rates,
                domains: constraints,
            } = self.normal_flow(&support, target, builder, &coefficient)?;
            domains.extend(constraints);
            support[index].order -= 1;

            // Differentiate energy powers through their complete directional
            // rates. Differentiating 1/E coordinate by coordinate would create
            // artificial 1/E³ terms whose cancellation loses precision at a
            // soft spectator, even though the normal flow preserves that locus.
            let mut energies = Vec::new();
            let mut opaque = coefficient.clone();
            for edge in self.graph.iter_edge_ids() {
                if self.graph.loop_momentum_basis.edge_signatures[edge]
                    .internal
                    .iter()
                    .all(|sign| sign.is_zero())
                {
                    // A constant bridge energy must not capture matching
                    // literals, including thermal order and orientation data.
                    continue;
                }
                let square = self.energy(edge, builder)?.pow(2);
                let variable = (0usize..)
                    .map(|index| {
                        symbol!(&format!("gammalooprs::fermi_surface_energy_square_{index}"))
                    })
                    .find(|&symbol| {
                        !coefficient.contains_symbol(symbol)
                            && !opaque.contains_symbol(symbol)
                            && !flow
                                .iter()
                                .flatten()
                                .any(|atom| atom.contains_symbol(symbol))
                            && !energy_square_rates
                                .values()
                                .any(|atom| atom.contains_symbol(symbol))
                    })
                    .unwrap();
                let next = opaque.replace_map(|part, _, out| {
                    if part == square.as_view() {
                        **out = Atom::var(variable);
                    }
                });
                if next != opaque {
                    energies.push((edge, square, variable));
                    opaque = next;
                }
            }
            let mut smooth_derivative = Atom::Zero;
            let mut divergence = Atom::Zero;
            for (edge, velocity) in self.graph.loop_momentum_basis.loop_edges.iter().zip(&flow) {
                for (axis, component) in velocity.iter().enumerate() {
                    let coordinate = Indeterminate::try_from(Self::spatial(*edge, axis)).unwrap();
                    smooth_derivative += component * opaque.derivative(coordinate.clone());
                    divergence += component.derivative(coordinate);
                }
            }
            for (edge, square, variable) in &energies {
                let derivative = opaque.derivative(*variable);
                if derivative.is_zero() {
                    continue;
                }
                let rate = if *edge == target {
                    2 * self.energy(target, builder)?
                } else if support.iter().any(|factor| factor.edge == *edge) {
                    continue;
                } else if let Some(rate) = energy_square_rates.get(edge) {
                    square * rate
                } else {
                    let momentum = self.momentum(*edge, builder)?;
                    let routing = &self.graph.loop_momentum_basis.edge_signatures[*edge];
                    2 * Atom::add_many(flow.iter().enumerate().map(|(index, velocity)| {
                        Self::sign_atom(routing.internal[LoopIndex(index)])
                            * Atom::add_many(velocity.iter().zip(&momentum).map(|(v, q)| v * q))
                    }))
                };
                smooth_derivative += rate * derivative;
            }
            for (_, square, variable) in &energies {
                smooth_derivative = smooth_derivative
                    .replace(Atom::var(*variable))
                    .with(square.clone());
            }
            Self::add_sector(
                &mut pending,
                support.clone(),
                -smooth_derivative - divergence * &coefficient,
            );

            let mut occupation_calls: BTreeMap<FermiDistribution, BTreeSet<Atom>> = BTreeMap::new();
            let mut calls = BTreeSet::new();
            coefficient.as_view().visitor(&mut |part| {
                if Self::distribution(part)
                    .ok()
                    .flatten()
                    .is_some_and(|factor| factor.order == 0)
                {
                    calls.insert(part.to_owned());
                }
                true
            });
            for call in calls {
                let mut factor = Self::distribution(call.as_view())?.unwrap();
                if self.graph[factor.edge].chemical_potential_atom().is_none() {
                    continue;
                }
                if !self.graph[factor.edge].is_fermion() {
                    return Err(eyre!(
                        "Fermi-surface contacts require a fermion, edge {}",
                        factor.edge
                    ));
                }
                factor.order = 1;
                occupation_calls.entry(factor).or_default().insert(call);
            }
            for (factor, calls) in occupation_calls {
                // Differentiate both thermal signs in one pass. A common
                // shift preserves A*dP inside each grouped numerator; adding
                // separate partial derivatives distributes A across P and can
                // hide the exact P_empty = 0 cancellation at higher orders.
                let shift_symbol = (0usize..)
                    .map(|index| symbol!(&format!("gammalooprs::fermi_surface_shift_{index}")))
                    .find(|&symbol| !coefficient.contains_symbol(symbol))
                    .unwrap();
                let shift = Atom::var(shift_symbol);
                let shifted = coefficient.replace_map(|part, _, out| {
                    if calls.iter().any(|call| call.as_view() == part) {
                        **out = part.to_owned() + &shift;
                    }
                });
                let derivative = shifted.derivative(shift_symbol).replace(shift).with(0);
                if derivative.is_zero() {
                    continue;
                }
                let other_q = self.momentum(factor.edge, builder)?;
                let signature = &self.graph.loop_momentum_basis.edge_signatures[factor.edge];
                let normal = Atom::add_many(flow.iter().enumerate().map(|(index, velocity)| {
                    Self::sign_atom(signature.internal[LoopIndex(index)])
                        * Atom::add_many(velocity.iter().zip(&other_q).map(|(a, b)| a * b))
                })) / self.energy(factor.edge, builder)?;
                if normal.is_zero() {
                    continue;
                }
                let mut contact = support.clone();
                contact.push(factor);
                Self::add_sector(&mut pending, contact, -normal * derivative);
            }
        }
        Ok((result, domains.into_iter().collect()))
    }

    fn localize_support(
        &self,
        support: &[FermiDistribution],
        coefficient: Atom,
        builder: &ParamBuilder,
    ) -> Result<Atom> {
        let basis = self.basis(support)?;
        let mut shifts = vec![
            [Atom::Zero, Atom::Zero, Atom::Zero];
            self.graph.loop_momentum_basis.loop_edges.len()
        ];
        let mut jacobian = Atom::one();
        let mut guards = Vec::new();
        for factor in support {
            debug_assert_eq!(factor.order, 1);
            let target_loop = LoopIndex(
                basis
                    .loop_edges
                    .iter()
                    .position(|edge| *edge == factor.edge)
                    .unwrap(),
            );
            let q = self.momentum(factor.edge, builder)?;
            let radius2 = Self::norm_squared(&q);
            let radius = radius2.sqrt();
            let mass = self.graph[factor.edge].mass_atom();
            let mu = self.graph[factor.edge].chemical_potential_atom().unwrap();
            let effective_mu = &factor.orientation * mu;
            let fermi_radius = (effective_mu.pow(2) - mass.pow(2)).sqrt();
            let scale = &fermi_radius / &radius;
            for (old_index, old_edge) in
                self.graph.loop_momentum_basis.loop_edges.iter().enumerate()
            {
                let routing =
                    Self::sign_atom(basis.edge_signatures[*old_edge].internal[target_loop]);
                for (axis, component) in q.iter().enumerate() {
                    shifts[old_index][axis] += &routing * (&scale - 1) * component;
                }
            }
            // Insert int_0^infinity h(t)dt = 1, h(t)=(1+t)^(-2).
            // t^3 / |dE(tq)/dt| times h(t) at t=k_F/|q| is this density.
            jacobian *=
                &effective_mu * fermi_radius.pow(2) / (&radius2 * (&radius + &fermi_radius).pow(2));
            let gap = effective_mu - function!(Symbol::ABS, mass);
            if gap.is_zero() {
                return Err(eyre!(
                    "Fermi-surface localization on edge {} requires a nondegenerate shell; the zero-radius threshold needs a separate distributional limit",
                    factor.edge
                ));
            }
            guards.push(gap);
        }
        let replacements = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .enumerate()
            .flat_map(|(index, edge)| {
                let shifts = &shifts[index];
                (0..3).map(move |axis| {
                    let coordinate = Self::spatial(*edge, axis);
                    (coordinate.clone(), coordinate + &shifts[axis])
                })
            })
            .collect::<BTreeMap<_, _>>();
        let coefficient = self.expand_functions(coefficient, builder, true)?;
        let localized = coefficient.replace_map(|part, _, out| {
            if let Some(value) = replacements.get(&part.to_owned()) {
                **out = value.clone();
            }
        });
        let mut value = jacobian * localized;
        for gap in guards {
            // IF is lazy: absent shells never evaluate the square root or
            // localized coefficient. Model validation diagnoses a zero-radius
            // shell; raw standalone inputs must also fail with a nonfinite
            // value there instead of silently replacing that limit by zero.
            let positive = (1 + &gap / function!(Symbol::ABS, &gap)) / 2;
            value = Symbol::IF.call_args([
                gap.clone(),
                Symbol::IF.call_args([positive, value, Atom::Zero]),
                gap.pow(-1),
            ]);
        }
        Ok(value)
    }
}

#[cfg(test)]
#[path = "fermi_surface_tests.rs"]
mod tests;
