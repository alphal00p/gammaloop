//! Independent bounded-jet references for small physical OS forests.
//! Source numerators remain factorized; only scalar denominator polynomials
//! and the Taylor coefficients demanded by a component are materialized.

use super::on_shell_tests::{SpinReduction, assert_rational_zero, scalar_component};
use super::*;
use crate::uv::approx::integrated;
use crate::{
    graph::parse::NumeratorInterface,
    initialisation::test_initialise,
    uv::{Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode},
};
use idenso::{
    shorthands::UndoShorthands,
    tensor::{SymbolicNetExt, SymbolicNetParse},
};
use linnet::half_edge::{involution::Hedge, subgraph::ModifySubSet};
use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::IndexedStructure,
};

#[derive(Clone, Copy)]
pub(crate) struct Component {
    edges: &'static [usize],
    degree: i32,
    scheme: ApproximationType,
    /// Known physical fixture boundary, independent of production's search.
    on_shell: Option<(usize, usize)>, // incoming edge, outgoing spinor hedge
}

impl Component {
    pub(crate) fn node(self, graph: &Graph) -> Result<OwnedForestNode> {
        let mut selected = graph.empty_subgraph::<SuBitGraph>();
        for edge in self.edges {
            selected.add(graph[&EdgeIndex(*edge)].1);
        }
        Ok(OwnedForestNode {
            spinney: Spinney::with_scheme(
                InternalSubGraph::cleaned_filter_optimist(selected, graph.as_ref()),
                graph,
                &graph.loop_momentum_basis,
                self.scheme,
                self.degree,
            )?
            .expect("the selected component is UV divergent"),
            topo_order: 1,
        })
    }
}

struct ForestOracle {
    graph: Graph,
    masses: Vec<Atom>,
    velocity: [Atom; 4],
}

impl ForestOracle {
    fn new(source: &str) -> Result<Self> {
        let mut graphs = Graph::import_dot(source, &crate::utils::load_generic_model("sm"))?;
        assert_eq!(graphs.len(), 1, "the forest fixture contains one graph");
        let graph = graphs.remove(0);
        let masses = graph
            .underlying
            .iter_edges()
            .map(|(_, edge, _)| graph[edge].mass_atom(&graph.model))
            .filter(|mass| !mass.is_zero())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(Self {
            graph,
            masses,
            velocity: [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero],
        })
    }

    fn root(&self) -> OwnedForestNode {
        OwnedForestNode {
            spinney: Spinney::empty(&self.graph),
            topo_order: 0,
        }
    }

    fn route(&self, atom: &Atom, node: &OwnedForestNode) -> Atom {
        // The independent local oracle contracts finite four-dimensional
        // components; production retains symbolic d for analytic integration.
        atom.replace_multiple(self.graph.uv_wrapped_replacement(
            node.subgraph(),
            node.lmb(),
            &[W_.x___],
        ))
        .replace(GS.dim)
        .with(Atom::num(4))
    }

    fn source(&self, node: &OwnedForestNode, given: &OwnedForestNode) -> Atom {
        let reduced = node.reduced_subgraph(given);
        let numerator = self
            .graph
            .numerator(&reduced, given.subgraph())
            .expect("the graph source has an admitted tensor interface")
            .to_d_dim(4)
            .expect("physical source dimension conversion preserves its interface")
            .color_simplify()
            .expect("the physical source color algebra is supported")
            .get_single_atom()
            .expect("the physical source is one factorized graph numerator");
        self.route(
            &(numerator / self.graph.denominator(&reduced, &self.graph.model, |_| 1)),
            node,
        )
    }

    fn external_edges(&self, node: &OwnedForestNode) -> Vec<EdgeIndex> {
        let owned = self.graph.paired_edges(node.subgraph());
        node.lmb()
            .ext_edges
            .iter_enumerated()
            .filter(|(external, _)| {
                owned
                    .iter()
                    .any(|edge| node.lmb().edge_signatures[*edge].external[*external].is_sign())
            })
            .map(|(_, edge)| *edge)
            .collect()
    }

    fn scale(
        &self,
        atom: &Atom,
        node: &OwnedForestNode,
        variable: Symbol,
        ultraviolet: bool,
    ) -> Atom {
        let external = self.external_edges(node);
        let scale = Atom::var(variable);
        let scaled = GS.resolve_os_momenta(atom).replace_map(|view, _, output| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() > 0
                && usize::try_from(momentum.get(0))
                    .is_ok_and(|edge| external.contains(&EdgeIndex(edge)))
            {
                **output = view.to_owned() * &scale;
            }
        });
        if ultraviolet {
            // A free expansion mass produced by a completed child has the
            // same Taylor weight as a physical mass. Terminal denominators
            // contain mUV, which retains its hard weight independently.
            scaled
                .replace(GS.m_uv_expansion)
                .with(Atom::var(GS.m_uv_expansion) * &scale)
                .replace_multiple(
                    self.masses.iter().map(|mass| {
                        Replacement::new(mass.to_pattern(), (mass * &scale).to_pattern())
                    }),
                )
        } else {
            scaled
        }
    }

    fn coefficient(atom: &Atom, variable: Symbol, degree: usize) -> Atom {
        let derivative = match degree {
            0 => atom.clone(),
            1 => atom.derivative(variable),
            2 => atom.derivative(variable).derivative(variable) / Atom::num(2),
            _ => unreachable!("the independent small-forest jet needs degrees zero through two"),
        };
        derivative.replace(variable).with(Atom::Zero)
    }

    fn multiply_jets(left: &[Atom], right: &[Atom], degree: usize) -> Vec<Atom> {
        (0..=degree)
            .map(|power| {
                (0..=power).fold(Atom::Zero, |sum, index| {
                    sum + &left[index] * &right[power - index]
                })
            })
            .collect()
    }

    /// Independent direct-soft scaling, with no production series/Taylor call.
    /// An original denominator acquires -M²+t²*mUVexp² under U; a completed
    /// terminal vacuum denominator is already hard and is not rearranged again.
    fn jet(&self, atom: &Atom, node: &OwnedForestNode, ultraviolet: bool) -> Result<Atom> {
        let degree = if ultraviolet {
            node.dod()
        } else {
            node.dod() - 1
        };
        if degree < 0 {
            return Ok(Atom::Zero);
        }
        let degree = degree as usize;
        assert!(degree <= 2);
        let variable = symbol!("os_forest_reference::jet");
        let scale = Atom::var(variable);
        let vacuum_mass = Atom::var(GS.m_uv_vacuum).pow(2);
        let terminal_mass = Atom::var(GS.m_uv_expansion).pow(2);
        let routed = self.route(atom, node);
        let mut result = Atom::Zero;
        for term in FourDTerm::from_view(routed.as_view())? {
            // Scalar differentiation acts on the explicit vector coefficients
            // of gamma matrices, including compact arguments produced by a
            // completed child. Open those arguments without distributing the
            // surrounding graph numerator or evaluating its spin algebra.
            let explicit = NumeratorInterface::export(&term.numerator)
                .0
                .parse_to_symbolic_net::<Aind>(&ParseSettings {
                    shorthand_parsing: ShorthandParsing::Expand {
                        schoonschip: SchoonschipExpansionMode::full(),
                        trace: false,
                        chain: false,
                    },
                    ..Default::default()
                })?
                .simple_execute::<()>()?;
            let numerator = self.scale(&explicit, node, variable, ultraviolet);
            let mut coefficients = (0..=degree)
                .map(|power| Self::coefficient(&numerator, variable, power))
                .collect::<Vec<_>>();
            for denominator in term.denominators {
                let physical = scalar_component(&GS.resolve_os_momenta(&denominator.full_expr))?;
                let mut polynomial = self.scale(&physical, node, variable, ultraviolet);
                let is_terminal = denominator.mass_squared == terminal_mass;
                if ultraviolet && !is_terminal {
                    polynomial += &terminal_mass * scale.pow(2) - &vacuum_mass;
                }
                let constant = Self::coefficient(&polynomial, variable, 0).expand();
                let momentum = Self::coefficient(
                    &self.scale(&denominator.momentum, node, variable, ultraviolet),
                    variable,
                    0,
                );
                let mass = if ultraviolet {
                    terminal_mass.clone()
                } else {
                    denominator.mass_squared
                };
                let leading = GS.den(
                    usize::from(denominator.source_edge),
                    momentum,
                    mass,
                    constant,
                );
                let mut inverse = vec![leading.pow(-1)];
                if degree >= 1 {
                    let linear = Self::coefficient(&polynomial, variable, 1).expand();
                    inverse.push(-&linear / leading.pow(2));
                    if degree >= 2 {
                        let quadratic = Self::coefficient(&polynomial, variable, 2).expand();
                        inverse.push(linear.pow(2) / leading.pow(3) - quadratic / leading.pow(2));
                    }
                }
                coefficients = Self::multiply_jets(&coefficients, &inverse, degree);
            }
            result += coefficients.into_iter().sum::<Atom>();
        }
        assert!(
            !result.contains_symbol(Symbol::DERIVATIVE),
            "the independent bounded jet must differentiate physical coefficients, not metadata",
        );
        Ok(result)
    }

    fn on_shell(&self, atom: &Atom, node: &OwnedForestNode, boundary: (usize, usize)) -> Atom {
        let (incoming, outgoing) = boundary;
        let mass = self.graph[EdgeIndex(incoming)].mass_atom(&self.graph.model);
        let port = self.graph.underlying[Hedge(outgoing)]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom();
        let external = self.external_edges(node);
        let routed = self.route(atom, node);
        [1, -1].into_iter().fold(Atom::Zero, |sum, sigma| {
            let point = GS.on_shell_momentum(
                EdgeIndex(incoming),
                &(Atom::num(sigma) * &mass),
                &self.velocity,
            );
            let evaluated = routed.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() >= 1
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
            });
            let dummy = idenso::bis!(4, Atom::from(Aind::new_dummy()));
            // The explicit fixed-vector basis is contracted with the
            // Minkowski metric, which supplies the spatial minus signs.
            let lorentz = spenso::mink!(4, Atom::from(Aind::new_dummy()));
            let slash = &self.velocity[0] * idenso::gamma0!(&port, &dummy)
                + self
                    .velocity
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(index, value)| {
                        value
                            * function!(GS.delta_vec, GS.cind(index), &lorentz)
                            * idenso::gamma!(&port, &dummy, &lorentz)
                    })
                    .sum::<Atom>();
            let gamma = slash * evaluated.replace(port.to_pattern()).with(dummy);
            sum + (evaluated + Atom::num(sigma) * gamma) / Atom::num(2)
        })
    }

    fn operation(&self, atom: &Atom, node: &OwnedForestNode, component: Component) -> Result<Atom> {
        let ordinary = self.jet(atom, node, true)?;
        let refinement = match component.on_shell {
            Some(boundary) => self.on_shell(atom, node, boundary),
            None => self.jet(atom, node, false)?,
        };
        let overlap = self.jet(&refinement, node, true)?;
        Ok(ordinary + refinement - overlap)
    }

    fn assert_operator_branches(
        &self,
        label: &str,
        node: &OwnedForestNode,
        component: Component,
        ctx: &UVCtx<'_>,
        given: (&OwnedForestNode, &Full4dCts),
    ) -> Result<Atom> {
        let (given_node, prefix) = given;
        let actual_source = grow(prefix, ctx, node, given_node)?;
        // A completed prefix is used here only after its own independent
        // branch certificate. Congruence then composes that certificate with
        // this independently constructed cograph, before the outer operation.
        let prefix_tensor = prefix
            .atom()
            .with_checked_expression(GS.erase_uv_momentum_provenance(prefix.atom().expression()))?;
        let source = finalize(
            NumeratorInterface::admit(self.source(node, given_node), None)?
                .try_into_indexed()?
                .multiply(prefix_tensor)?,
        )?;
        let source = NumeratorInterface::export(&source).0;
        let lmb = compatible_lmb(ctx, node, given_node, prefix)?;
        assert_rational_zero(
            &format!("{label} source"),
            &self.graph,
            node,
            &(NumeratorInterface::export(&actual_source).0 - &source),
            SpinReduction::SymbolicFirst,
        )?;
        let actual_u = t_raw(&actual_source, ctx, node, given_node, &lmb)?;
        let reference_u = self.jet(&source, node, true)?;
        assert_rational_zero(
            &format!("{label} U branch"),
            &self.graph,
            node,
            &(NumeratorInterface::export(&actual_u).0 - reference_u),
            SpinReduction::SymbolicFirst,
        )?;
        let (actual_refinement, reference_refinement) = match component.on_shell {
            Some(boundary) => {
                let [positive, negative] = self.graph.on_shell_branches(
                    &self.graph.model,
                    &actual_source,
                    node.subgraph(),
                    &lmb,
                    &self.velocity,
                )?;
                (
                    positive.try_add(&negative)?,
                    self.on_shell(&source, node, boundary),
                )
            }
            None => (
                tilde_t_raw(&actual_source, ctx, node, given_node, &lmb)?,
                self.jet(&source, node, false)?,
            ),
        };
        assert_rational_zero(
            &format!("{label} refinement branch"),
            &self.graph,
            node,
            &(NumeratorInterface::export(&actual_refinement).0 - &reference_refinement),
            SpinReduction::SymbolicFirst,
        )?;
        let actual_overlap = t_raw(
            &actual_refinement
                .with_checked_expression(GS.resolve_os_momenta(actual_refinement.expression()))?,
            ctx,
            node,
            given_node,
            &lmb,
        )?;
        let reference_overlap = self.jet(&reference_refinement, node, true)?;
        assert_rational_zero(
            &format!("{label} overlap branch"),
            &self.graph,
            node,
            &(NumeratorInterface::export(&actual_overlap).0 - reference_overlap),
            SpinReduction::SymbolicFirst,
        )?;
        let completed = finalize(
            actual_u
                .try_add(&actual_refinement)?
                .try_sub(&actual_overlap)?,
        )?;
        Ok(NumeratorInterface::export(&completed).0)
    }

    /// The four completed prefixes below the outer OS operation, in forest
    /// order: empty, gamma, Gamma, Gamma-after-gamma. Only local terms enter.
    fn three_loop_prefixes(&self) -> Result<[(OwnedForestNode, Full4dCts); 4]> {
        let child = THREE_LOOP_CHILD.node(&self.graph)?;
        let parent = THREE_LOOP_PARENT.node(&self.graph)?;
        let root = self.root();
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&self.graph, &self.graph.model, &settings);
        let seed = Full4dCts::root();
        let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
        let child_prefix = Full4dCts::from_factorized_local(&child_local);
        let parent_local = uv_limit(&seed, &ctx, &parent, &root, &parent, &root)?;
        let nested_local = uv_limit(&child_prefix, &ctx, &parent, &child, &parent, &child)?;
        Ok([
            (root, seed),
            (child, child_prefix),
            (
                THREE_LOOP_PARENT.node(&self.graph)?,
                Full4dCts::from_factorized_local(&parent_local),
            ),
            (parent, Full4dCts::from_factorized_local(&nested_local)),
        ])
    }

    fn three_loop_families(
        &self,
        prefixes: &[(OwnedForestNode, Full4dCts); 4],
    ) -> Result<[Atom; 8]> {
        let whole = THREE_LOOP_WHOLE.node(&self.graph)?;
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&self.graph, &self.graph.model, &settings);
        let mut families = std::array::from_fn(|_| Atom::Zero);
        for (index, (given, prefix)) in prefixes.iter().enumerate() {
            families[2 * index] = NumeratorInterface::export(&grow(prefix, &ctx, &whole, given)?).0;
            families[2 * index + 1] = NumeratorInterface::export(
                uv_limit(prefix, &ctx, &whole, given, &whole, given)?.atom(),
            )
            .0;
        }
        Ok(families)
    }

    fn assert_three_loop_hard_strata(&self, strata: &[[usize; 3]]) -> Result<()> {
        let whole = THREE_LOOP_WHOLE.node(&self.graph)?;
        let families = self.three_loop_families(&self.three_loop_prefixes()?)?;
        let complete = Atom::add_many(&families);
        let parameter = Atom::var(symbol!("os_forest_reference::hard"));
        for blocks in strata {
            // Zero denotes a fixed carrier. Each positive label denotes an
            // independent hard vector; equal labels retain a fixed difference.
            // These partitions also cover vanishing hard gluon momenta k-p/l-k.
            let loops = *blocks.iter().max().unwrap();
            assert!(loops > 0);
            let mut momenta = Vec::new();
            for (index, block) in blocks.iter().enumerate() {
                if *block == 0 {
                    continue;
                }
                let representative = blocks
                    .iter()
                    .position(|candidate| candidate == block)
                    .unwrap();
                let mut momentum = function!(GS.emr_mom, representative + 2) / &parameter;
                if representative != index {
                    momentum += function!(GS.emr_mom, index + 2);
                }
                momenta.push((index + 2, momentum));
            }
            let [p, k, l] = *blocks;
            let hard_fermions =
                2 * usize::from(p != 0) + 2 * usize::from(k != 0) + usize::from(l != 0);
            let hard_gluons = usize::from(p != 0) + usize::from(k != p) + usize::from(l != k);
            let bare_degree = 4 * loops as i32 - hard_fermions as i32 - 2 * hard_gluons as i32;
            self.assert_hard_limit(
                &format!("three-loop hard blocks {blocks:?}"),
                (&complete, (bare_degree >= 0).then_some(&families[0])),
                &whole,
                &momenta,
                loops,
            )?;
        }
        Ok(())
    }

    fn physical_limit(atom: &Atom) -> Atom {
        // Match the completed physical-integrand boundary, not a raw Taylor
        // operation: ownership tags/markers are metadata, and the expansion
        // bookkeeping mass is identified with the auxiliary vacuum mass.
        let physical = GS
            .resolve_os_momenta(&GS.erase_uv_momentum_provenance(atom))
            .replace(function!(GS.ct_marker, W_.a_))
            .with(Atom::one())
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum)
            .replace(GS.dim_epsilon)
            .with(Atom::Zero)
            .replace(GS.integrated_loop_scale)
            .with(Atom::one());
        assert!(!physical.contains_symbol(GS.uv_momentum_provenance));
        assert!(
            !physical.contains_symbol(Symbol::DERIVATIVE),
            "a completed physical-limit copy must not contain metadata derivatives",
        );
        physical
    }

    fn hard_coefficients(
        &self,
        label: &str,
        atom: &Atom,
        node: &OwnedForestNode,
        momenta: &[(usize, Atom)],
        hard_loops: usize,
    ) -> Result<BTreeMap<symbolica::domains::rational::Rational, Atom>> {
        let started = std::time::Instant::now();
        debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "os_hard_coefficients_start", label, hard_loops,
            input_bytes = atom.as_view().get_byte_size(),
            "Constructing exact hard-limit coefficients"
        );
        let parameter = symbol!("os_forest_reference::hard");
        let scale = Atom::var(parameter);
        let physical = self.route(&Self::physical_limit(atom), node);
        let replacements = momenta
            .iter()
            .map(|(edge, momentum)| {
                Replacement::new(
                    GS.emr_mom(EdgeIndex(*edge), W_.x___).to_pattern(),
                    GS.indexed_momentum(momentum, &[Atom::var(W_.x___)])
                        .to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        let mut coefficients = BTreeMap::new();
        for (term_index, term) in FourDTerm::from_view(physical.as_view())?
            .into_iter()
            .enumerate()
        {
            let term_started = std::time::Instant::now();
            debug_tags!(#generation, #profile, #uv, #numerator, #summary;
                stage = "os_hard_term_start", label, term_index,
                numerator_bytes = term.numerator.expression().as_view().get_byte_size(),
                denominator_count = term.denominators.len(),
                "Constructing one hard-limit denominator class"
            );
            let mut scaled = NumeratorInterface::export(&term.numerator)
                .0
                .parse_to_symbolic_net::<Aind>(&ParseSettings {
                    shorthand_parsing: ShorthandParsing::Expand {
                        schoonschip: SchoonschipExpansionMode::full(),
                        trace: false,
                        chain: false,
                    },
                    ..Default::default()
                })?
                .simple_execute::<()>()?
                .replace_multiple(&replacements)
                / scale.pow(4 * hard_loops);
            for denominator in term.denominators {
                // Only denominator polynomials are expanded. Make each hard
                // t^-2 explicit before the bounded numerator/denominator jet.
                let polynomial = scalar_component(&denominator.full_expr)?
                    .replace_multiple(&replacements)
                    .expand();
                let hard_power = if polynomial.contains_symbol(parameter) {
                    2
                } else {
                    0
                };
                let polynomial = (polynomial * scale.pow(hard_power)).expand();
                let constant = polynomial.replace(parameter).with(Atom::Zero);
                assert!(
                    !constant.is_zero() && !constant.contains_symbol(parameter),
                    "{label}: a generic hard propagator must have a nonzero quadratic leading term",
                );
                scaled *= scale.pow(hard_power)
                    / GS.den(
                        usize::from(denominator.source_edge),
                        denominator.momentum,
                        denominator.mass_squared,
                        polynomial,
                    );
            }
            let jet = scaled.series_preserving_factors(parameter, Atom::Zero.as_view(), 0, &[])?;
            for (power, coefficient) in jet.series(parameter, Atom::Zero, 0)?.terms() {
                *coefficients.entry(power).or_insert(Atom::Zero) += coefficient;
            }
            debug_tags!(#generation, #profile, #uv, #numerator, #summary;
                stage = "os_hard_term_done", label, term_index,
                elapsed_ms = term_started.elapsed().as_secs_f64() * 1000.0,
                jet_bytes = jet.as_view().get_byte_size(),
                "Constructed one hard-limit denominator class"
            );
        }
        debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "os_hard_coefficients_done", label, hard_loops,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            coefficient_count = coefficients.len(),
            coefficient_bytes = coefficients.values().map(|atom| atom.as_view().get_byte_size()).sum::<usize>(),
            "Constructed exact hard-limit coefficients"
        );
        Ok(coefficients)
    }

    fn assert_hard_limit(
        &self,
        label: &str,
        sources: (&Atom, Option<&Atom>),
        node: &OwnedForestNode,
        momenta: &[(usize, Atom)],
        hard_loops: usize,
    ) -> Result<()> {
        let started = std::time::Instant::now();
        debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "os_hard_stratum_start", label, hard_loops,
            "Certifying an exact hard-limit stratum"
        );
        let (complete, bare_control) = sources;
        for (power, coefficient) in
            self.hard_coefficients(label, complete, node, momenta, hard_loops)?
        {
            let coefficient_started = std::time::Instant::now();
            debug_tags!(#generation, #profile, #uv, #numerator, #summary;
                stage = "os_hard_coefficient_start", label, power = %power,
                coefficient_bytes = coefficient.as_view().get_byte_size(),
                "Certifying one exact hard-limit coefficient"
            );
            assert_rational_zero(
                &format!("{label}: measure-weighted inverse-hard coefficient {power}"),
                &self.graph,
                node,
                &coefficient,
                SpinReduction::SymbolicFirst,
            )?;
            debug_tags!(#generation, #profile, #uv, #numerator, #summary;
                stage = "os_hard_coefficient_done", label, power = %power,
                elapsed_ms = coefficient_started.elapsed().as_secs_f64() * 1000.0,
                "Certified one exact hard-limit coefficient"
            );
        }
        if let Some(bare) = bare_control {
            let bare = self.hard_coefficients(label, bare, node, momenta, hard_loops)?;
            let (_, leading) = bare
                .first_key_value()
                .expect("the bare graph has a divergent jet");
            let terms = FourDTerm::from_view(leading.as_view())?;
            assert_eq!(
                terms.len(),
                1,
                "the leading bare coefficient has one denominator shell"
            );
            let tensor = crate::uv::overlap_numerators::NumeratorJetCache::default().canonicalize(
                &NumeratorInterface::export(&terms[0].numerator)
                    .0
                    .undo_dots::<Aind>()?,
            )?;
            let mut nonzero = false;
            tensor.key.visitor(&mut |view| {
                let AtomView::Fun(entry) = view else {
                    return true;
                };
                if entry.get_symbol() != symbol!("gammalooprs::uv::overlap_entry") {
                    return true;
                }
                // One exact nonzero finite entry proves the control is
                // nonzero. Coordinate samples are used only for this
                // negative control, never for the UV zero certificate.
                let value = entry
                    .get(1)
                    .to_owned()
                    .replace(idenso::color_cas!(2, idenso::cof!(3)))
                    .with(Atom::num(4) / Atom::num(3))
                    .replace(idenso::color_idx!(2, idenso::cof!(3)))
                    .with(Atom::one() / Atom::num(2))
                    .replace_map(|view, _, output| match view {
                        AtomView::Var(_) => **output = Atom::one(),
                        AtomView::Fun(momentum)
                            if momentum.get_symbol() == GS.emr_mom && momentum.get_nargs() == 2 =>
                        {
                            let Ok(edge) = usize::try_from(momentum.get(0)) else {
                                return;
                            };
                            let AtomView::Fun(index) = momentum.get(1) else {
                                return;
                            };
                            let Ok(component) = usize::try_from(index.get(0)) else {
                                return;
                            };
                            **output =
                                Atom::num(1 + (edge + 2) * (component + 1) + component * component);
                        }
                        _ => {}
                    });
                nonzero |= matches!(value.as_view(), AtomView::Num(_)) && !value.is_zero();
                !nonzero
            });
            assert!(
                nonzero,
                "{label}: the bare divergent jet needs an exact nonzero control"
            );
        }
        debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "os_hard_stratum_done", label, hard_loops,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Certified an exact hard-limit stratum and its bare control"
        );
        Ok(())
    }

    fn nested_families(
        &self,
        child_spec: Component,
        parent_spec: Component,
        hard_carriers: [usize; 2],
    ) -> Result<Atom> {
        let child = child_spec.node(&self.graph)?;
        let parent = parent_spec.node(&self.graph)?;
        let root = self.root();
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&self.graph, &self.graph.model, &settings);
        let seed = Full4dCts::root();
        let bare_child = self.source(&child, &root);
        let child_reference = -self.operation(&bare_child, &child, child_spec)?;
        let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
        let certified_child =
            self.assert_operator_branches("child", &child, child_spec, &ctx, (&root, &seed))?;
        // Each branch already equals its independent bounded-jet reference.
        // Certify the one forest minus sign against their exact composition.
        assert_rational_zero(
            "nested forest child is minus its independently certified branch sum",
            &self.graph,
            &child,
            &(NumeratorInterface::export(child_local.atom()).0 + certified_child),
            SpinReduction::SymbolicFirst,
        )?;
        let completed_child = Full4dCts::from_factorized_local(&child_local);
        if child_spec.scheme == ApproximationType::OS && parent_spec.scheme == ApproximationType::IR
        {
            // The completed child depends on the parent's loop p and on its
            // own fixed rest-frame points, but not on the parent's external q.
            // Check this boundary before attaching the three-line cograph.
            let lmb = compatible_lmb(&ctx, &parent, &child, &completed_child)?;
            let actual_soft = tilde_t_raw(child_local.atom(), &ctx, &parent, &child, &lmb)?;
            let reference_soft = self.jet(&child_reference, &parent, false)?;
            assert_rational_zero(
                "the parent soft jet leaves the independent OS child fixed",
                &self.graph,
                &parent,
                &(&reference_soft - &child_reference),
                SpinReduction::SymbolicFirst,
            )?;
            assert_rational_zero(
                "parent S of the actual completed OS child matches the independent jet",
                &self.graph,
                &parent,
                &(NumeratorInterface::export(&actual_soft).0 - reference_soft),
                SpinReduction::SymbolicFirst,
            )?;
        }

        let bare_parent = self.source(&parent, &root);
        let cograph = self.source(&parent, &child);
        let nested_reference =
            -self.operation(&(&child_reference * &cograph), &parent, parent_spec)?;
        let parent_local = uv_limit(&seed, &ctx, &parent, &root, &parent, &root)?;
        let certified_parent =
            self.assert_operator_branches("parent", &parent, parent_spec, &ctx, (&root, &seed))?;
        let nested_local = uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child)?;
        // Compose the independently certified open-child tensor with an
        // independently certified cograph. Recontracting their difference
        // after joining all boundary indices can hide an already proved zero
        // inside products of finite tensor entries.
        let production_cograph = grow_factor(&ctx, &parent, &child)?;
        assert_rational_zero(
            "nested forest cograph",
            &self.graph,
            &parent,
            &(NumeratorInterface::export(&production_cograph).0 - &cograph),
            SpinReduction::SymbolicFirst,
        )?;
        let child_insertion = grow(&completed_child, &ctx, &parent, &child)?;
        assert_eq!(
            child_insertion,
            finalize(child_local.atom().clone().multiply(production_cograph)?)?,
            "the child forest family inserts the complete signed child exactly once",
        );
        assert_rational_zero(
            "nested forest parent is minus its independently certified branch sum",
            &self.graph,
            &parent,
            &(NumeratorInterface::export(parent_local.atom()).0 + certified_parent),
            SpinReduction::SymbolicFirst,
        )?;
        assert_rational_zero(
            "nested forest completed-child family",
            &self.graph,
            &parent,
            &(NumeratorInterface::export(nested_local.atom()).0 - &nested_reference),
            SpinReduction::SymbolicFirst,
        )?;
        if child_spec.scheme == ApproximationType::OS
            && parent_spec.scheme == ApproximationType::IR
            && self.graph.has_qcd_self_energy_insertion_ward_identity(
                &self.graph.model,
                parent.subgraph(),
                child.subgraph(),
            )
        {
            let completed_source = &child_reference * &cograph;
            let ordinary = self.jet(&completed_source, &parent, true)?;
            let soft = self.jet(&completed_source, &parent, false)?;
            let overlap = self.jet(&soft, &parent, true)?;
            let analytic = nested_local
                .active_sectors()
                .iter()
                .fold(Atom::Zero, |sum, sector| {
                    sum + NumeratorInterface::export(sector.integration_atom()).0
                });
            assert_rational_zero(
                "signed nested analytic input is minus U plus US of the completed child",
                &self.graph,
                &parent,
                &(analytic + ordinary - overlap),
                SpinReduction::SymbolicFirst,
            )?;
        }
        let production_source = grow(&seed, &ctx, &parent, &root)?;
        // Export only the independently certified families used by the physical
        // finite-component oracle; the native forest retains its tensor carriers.
        let signed_families = [
            NumeratorInterface::export(&production_source).0,
            NumeratorInterface::export(&child_insertion).0,
            NumeratorInterface::export(parent_local.atom()).0,
            NumeratorInterface::export(nested_local.atom()).0,
        ];
        let complete =
            &signed_families[0] + &signed_families[1] + &signed_families[2] + &signed_families[3];
        assert_eq!(
            complete,
            Atom::add_many(&signed_families),
            "the complete forest adds each certified signed family once without a second depth sign",
        );
        let parameter = Atom::var(symbol!("os_forest_reference::hard"));
        let [parent_carrier, child_carrier] = hard_carriers;
        assert!(parent.lmb().loop_edges.contains(&EdgeIndex(parent_carrier)));
        assert!(child.lmb().loop_edges.contains(&EdgeIndex(child_carrier)));
        let p = function!(GS.emr_mom, parent_carrier);
        let k = function!(GS.emr_mom, child_carrier);
        for (label, momenta, loops) in [
            ("hard child loop", vec![(child_carrier, &k / &parameter)], 1),
            (
                "hard parent loop",
                vec![(parent_carrier, &p / &parameter)],
                1,
            ),
            (
                "both generic hard loops",
                vec![
                    (parent_carrier, &p / &parameter),
                    (child_carrier, &k / &parameter),
                ],
                2,
            ),
            (
                "correlated hard k=p+r",
                vec![
                    (parent_carrier, &p / &parameter),
                    (child_carrier, &p / &parameter + &k),
                ],
                1,
            ),
        ] {
            let divergent_control = label == "hard child loop"
                || label == "both generic hard loops"
                || (label == "correlated hard k=p+r"
                    && parent_spec.scheme == ApproximationType::IR);
            self.assert_hard_limit(
                label,
                (&complete, divergent_control.then_some(&bare_parent)),
                &parent,
                &momenta,
                loops,
            )?;
        }
        // Raw U/O/S identities above retain independent bookkeeping masses;
        // the completed physical limit uses the common finalization boundary.
        if parent_spec.scheme == ApproximationType::IR {
            assert_rational_zero(
                "the complete IR-parent forest has no constant or linear external soft jet",
                &self.graph,
                &parent,
                &self.jet(&Self::physical_limit(&complete), &parent, false)?,
                SpinReduction::SymbolicFirst,
            )?;
        }
        if child_spec.scheme == ApproximationType::IR {
            assert_rational_zero(
                "the IR gluon child remainder has no constant or linear soft jet",
                &self.graph,
                &child,
                &self.jet(
                    &Self::physical_limit(
                        &(bare_child + NumeratorInterface::export(child_local.atom()).0),
                    ),
                    &child,
                    false,
                )?,
                SpinReduction::SymbolicFirst,
            )?;
        }
        // The bare source, signed child insertion, signed parent and signed
        // completed-child parent have separate exact certificates above.
        // Linearity composes those certificates into the full four-family
        // identity without re-expanding a difference of complete forests.
        Ok(complete)
    }
}

const QUARK_CHILD: Component = Component {
    edges: &[3, 4],
    degree: 1,
    scheme: ApproximationType::OS,
    on_shell: Some((2, 8)),
};
const QUARK_PARENT: Component = Component {
    edges: &[2, 3, 4, 5, 6],
    degree: 1,
    scheme: ApproximationType::OS,
    on_shell: Some((0, 1)),
};
const GLUON_CHILD: Component = Component {
    edges: &[3, 4],
    degree: 2,
    scheme: ApproximationType::IR,
    on_shell: None,
};
const GLUON_PARENT: Component = Component {
    edges: &[2, 3, 4, 5, 6],
    degree: 2,
    scheme: ApproximationType::IR,
    on_shell: None,
};

const THREE_LOOP_CHILD: Component = Component {
    edges: &[4, 5],
    degree: 1,
    scheme: ApproximationType::OS,
    on_shell: Some((3, 10)),
};
const THREE_LOOP_PARENT: Component = Component {
    edges: &[3, 4, 5, 6, 7],
    degree: 1,
    scheme: ApproximationType::OS,
    on_shell: Some((2, 14)),
};
const THREE_LOOP_WHOLE: Component = Component {
    edges: &[2, 3, 4, 5, 6, 7, 8, 9],
    degree: 1,
    scheme: ApproximationType::OS,
    on_shell: Some((0, 1)),
};

pub(crate) const FIG10_LOGARITHMIC_BOX: Component = Component {
    edges: &[2, 3, 5, 6],
    degree: 0,
    scheme: ApproximationType::MUV,
    on_shell: None,
};

const DISJOINT_LEFT: Component = Component {
    edges: &[2, 3],
    degree: 1,
    scheme: ApproximationType::OS,
    on_shell: Some((0, 6)),
};
const DISJOINT_RIGHT: Component = Component {
    edges: &[5, 6],
    degree: 1,
    scheme: ApproximationType::OS,
    on_shell: Some((4, 1)),
};

pub(crate) fn nested_quark_fixture(
    soft_child: bool,
) -> Result<(Graph, OwnedForestNode, OwnedForestNode, OwnedForestNode)> {
    let source = if soft_child {
        include_str!(
            "../../../../../tests/resources/graphs/local_soft_gluon_in_top_self_energy.dot"
        )
    } else {
        include_str!("../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot")
    };
    let oracle = ForestOracle::new(source)?;
    let parent = QUARK_PARENT.node(&oracle.graph)?;
    let child = (if soft_child { GLUON_CHILD } else { QUARK_CHILD }).node(&oracle.graph)?;
    let root = oracle.root();
    Ok((oracle.graph, parent, child, root))
}

pub(crate) fn three_loop_quark_fixture()
-> Result<(Graph, OwnedForestNode, OwnedForestNode, OwnedForestNode)> {
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"
    ))?;
    let mut parent = THREE_LOOP_PARENT.node(&oracle.graph)?;
    parent.topo_order = 2;
    let child = THREE_LOOP_CHILD.node(&oracle.graph)?;
    let root = oracle.root();
    Ok((oracle.graph, parent, child, root))
}

pub(crate) fn disjoint_quark_fixture()
-> Result<(Graph, OwnedForestNode, OwnedForestNode, OwnedForestNode)> {
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_disjoint_os_top_self_energies.dot"
    ))?;
    let left = DISJOINT_LEFT.node(&oracle.graph)?;
    let right = DISJOINT_RIGHT.node(&oracle.graph)?;
    let root = oracle.root();
    Ok((oracle.graph, left, right, root))
}

pub(crate) fn nested_gluon_fixture()
-> Result<(Graph, OwnedForestNode, OwnedForestNode, OwnedForestNode)> {
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
    ))?;
    let parent = GLUON_PARENT.node(&oracle.graph)?;
    let child = QUARK_CHILD.node(&oracle.graph)?;
    let root = oracle.root();
    Ok((oracle.graph, parent, child, root))
}

#[test]
fn nested_child_scalar_soft_residues_distinguish_physical_and_os_points() -> Result<()> {
    test_initialise()?;
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
    ))?;
    let child = QUARK_CHILD.node(&oracle.graph)?;
    let root = oracle.root();
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&oracle.graph, &oracle.graph.model, &settings);
    let seed = Full4dCts::root();
    let source = grow(&seed, &ctx, &child, &root)?;
    let lmb = compatible_lmb(&ctx, &child, &root, &seed)?;
    let [positive, negative] = oracle.graph.on_shell_branches(
        &oracle.graph.model,
        &source,
        child.subgraph(),
        &lmb,
        &settings.on_shell_velocity()?,
    )?;
    let mass = oracle.graph[EdgeIndex(3)].mass_atom(&oracle.graph.model);
    let radial_symbol = symbol!("os_forest_reference::soft_radius");
    let radius = Atom::var(radial_symbol);
    let energy_symbol = symbol!("os_forest_reference::child_energy");
    let energy = Atom::var(energy_symbol);
    // Unit rays and m=12 give exact nonzero soft coefficients. The moving
    // physical boundary is p=(15,9,0,0), while each OS branch fixes the whole
    // component boundary at (+/-12,0,0,0). These are distinct soft strata.
    for (label, input, boundary_energy, centre, ray, shell_energy) in [
        ("positive OS point", positive, 12, [0, 0, 0], [1, 0, 0], 12),
        ("negative OS point", negative, -12, [0, 0, 0], [0, 1, 0], 12),
        (
            "moving physical point, parallel",
            source.clone(),
            15,
            [9, 0, 0],
            [1, 0, 0],
            15,
        ),
        (
            "moving physical point, transverse",
            source,
            15,
            [9, 0, 0],
            [0, 1, 0],
            15,
        ),
    ] {
        let spatial = std::array::from_fn::<_, 3, _>(|index| {
            Atom::num(centre[index]) + Atom::num(ray[index]) * &radius
        });
        let input = oracle.route(
            &GS.resolve_os_momenta(
                &GS.erase_uv_momentum_provenance(&NumeratorInterface::export(&input).0),
            ),
            &child,
        );
        let terms = FourDTerm::from_view(input.as_view())?;
        assert_eq!(
            terms.len(),
            1,
            "each source or OS branch has one denominator product"
        );
        assert_eq!(terms[0].denominators.len(), 2);
        let energy_squared =
            Atom::num(144) + spatial.iter().map(|value| value.pow(2)).sum::<Atom>();
        let boundary = Atom::num(boundary_energy);
        for denominator in &terms[0].denominators {
            let mut polynomial = scalar_component(&denominator.full_expr)?
                .replace(mass.to_pattern())
                .with(12)
                .replace(GS.emr_mom(EdgeIndex(3), GS.cind(0)))
                .with(energy.to_pattern())
                .replace(GS.emr_mom(EdgeIndex(2), GS.cind(0)))
                .with(boundary.to_pattern());
            for index in 1..=3 {
                polynomial = polynomial
                    .replace(GS.emr_mom(EdgeIndex(3), GS.cind(index)))
                    .with(spatial[index - 1].to_pattern())
                    .replace(GS.emr_mom(EdgeIndex(2), GS.cind(index)))
                    .with(centre[index - 1]);
            }
            let expected = match usize::from(denominator.source_edge) {
                3 => energy.pow(2) - &energy_squared,
                4 => (&energy - &boundary).pow(2) - radius.pow(2),
                _ => unreachable!("the physical child has only quark e3 and gluon e4"),
            };
            assert!(
                (polynomial - expected).expand().is_zero(),
                "{label}: exact child pole polynomial"
            );
        }
        // The scalar residue oracle differentiates the explicit linear vector
        // coefficients, including slashes retained compactly by native rules.
        let numerator = NumeratorInterface::export(&terms[0].numerator)
            .0
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::full(),
                    trace: false,
                    chain: false,
                },
                ..Default::default()
            })?
            .simple_execute::<()>()?
            .replace(mass.to_pattern())
            .with(12)
            .replace_map(|view, _, output| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() != 2 {
                    return;
                }
                let Ok(owner) = usize::try_from(momentum.get(0)) else {
                    return;
                };
                if owner != 3 {
                    return;
                }
                let index = momentum.get(1);
                **output = function!(GS.delta_vec, GS.cind(0), index) * &energy
                    + (1..=3)
                        .map(|component| {
                            function!(GS.delta_vec, GS.cind(component), index)
                                * &spatial[component - 1]
                        })
                        .sum::<Atom>();
            });
        assert!(
            numerator
                .derivative(energy_symbol)
                .derivative(energy_symbol)
                .is_zero()
        );
        assert!(!numerator.derivative(energy_symbol).is_zero());

        // Clockwise Feynman residues at z=E and z=p0+rho, separately for
        // the two coefficients of the certified affine energy numerator.
        // Only these small scalar contour weights are Taylor expanded.
        let pole = energy_squared.sqrt();
        let quark_weight =
            -Atom::one() / (Atom::num(2) * &pole * ((&boundary - &pole).pow(2) - radius.pow(2)));
        let gluon_weight =
            -Atom::one() / (Atom::num(2) * &radius * ((&boundary + &radius).pow(2) - pole.pow(2)));
        let residues = [
            &quark_weight + &gluon_weight,
            &pole * quark_weight + (&boundary + &radius) * gluon_weight,
        ];
        let ray_product: i32 = centre
            .iter()
            .zip(ray)
            .map(|(component, direction)| component * direction)
            .sum();
        let leading = Atom::one() / Atom::num(4 * (shell_energy + ray_product));
        for (coefficient, residue) in residues.into_iter().enumerate() {
            let weighted = (radius.pow(2) * residue).series(radial_symbol, Atom::Zero, 0)?;
            assert_eq!(
                weighted.get_trailing_exponent(),
                0,
                "{label}: scalar kernel must be rho^-2 integrable"
            );
            assert_eq!(
                weighted.coefficient(0.into()),
                Some(if coefficient == 0 {
                    leading.clone()
                } else {
                    &boundary * &leading
                }),
                "{label}: the leading radial weight must remain nonzero"
            );
        }
    }
    // These are child scalar-kernel certificates. Complete nested CFF sums
    // still require their own radial checks after all parent residues.
    Ok(())
}

#[test]
fn parent_os_promotes_its_inherited_child_soft_carrier_to_a_physical_source() -> Result<()> {
    test_initialise()?;
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
    ))?;
    let parent = QUARK_PARENT.node(&oracle.graph)?;
    let child = QUARK_CHILD.node(&oracle.graph)?;
    let root = oracle.root();
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&oracle.graph, &oracle.graph.model, &settings);
    let child_local = uv_limit(&Full4dCts::root(), &ctx, &child, &root, &child, &root)?;
    let completed_child = Full4dCts::from_factorized_local(&child_local);
    let source = grow(&completed_child, &ctx, &parent, &child)?;
    let lmb = compatible_lmb(&ctx, &parent, &child, &completed_child)?;
    let owner_two = |atom: &Atom| {
        let mut records = BTreeSet::new();
        atom.visitor(&mut |view| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() > 0
                && let Some((owner, role, payload)) =
                    GS.uv_momentum_provenance_data(momentum.get(0))
                && owner == EdgeIndex(2)
            {
                records.insert((
                    role,
                    payload,
                    momentum
                        .iter()
                        .skip(1)
                        .map(|index| index.to_owned())
                        .collect::<Vec<_>>(),
                ));
                return false;
            }
            true
        });
        records
    };
    let inherited = owner_two(source.expression())
        .into_iter()
        .filter(|(role, _, _)| *role == UvMomentumProvenanceRole::DenominatorDerivedSoft)
        .collect::<Vec<_>>();
    assert!(
        !inherited.is_empty(),
        "the physical child must carry its soft crown owner e2"
    );
    let branches = oracle.graph.on_shell_branches(
        &oracle.graph.model,
        &source,
        parent.subgraph(),
        &lmb,
        &settings.on_shell_velocity()?,
    )?;
    for branch in &branches {
        let records = owner_two(branch.expression());
        assert!(
            records
                .iter()
                .all(|(role, _, _)| { *role != UvMomentumProvenanceRole::DenominatorDerivedSoft })
        );
        for (_, payload, indices) in &inherited {
            assert!(
                records.contains(&(
                    UvMomentumProvenanceRole::PhysicalSourceFixed,
                    payload.clone(),
                    indices.clone(),
                )),
                "parent O must retain the promoted carrier's owner, payload and component indices"
            );
        }
    }
    let reference = oracle.on_shell(
        &GS.erase_uv_momentum_provenance(&NumeratorInterface::export(&source).0),
        &parent,
        QUARK_PARENT.on_shell.unwrap(),
    );
    assert_rational_zero(
        "parent O provenance promotion leaves the independent owned momentum operation unchanged",
        &oracle.graph,
        &parent,
        &(branches
            .iter()
            .map(|branch| NumeratorInterface::export(branch).0)
            .sum::<Atom>()
            - reference),
        SpinReduction::SymbolicFirst,
    )
}

#[test]
fn os_quark_in_ir_gluon_forest_matches_independent_component_jets() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
    ))?
    .nested_families(QUARK_CHILD, GLUON_PARENT, [2, 3])
    .map(|_| ())
}

#[test]
fn os_quark_in_os_quark_forest_matches_independent_component_jets() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
    ))?
    .nested_families(QUARK_CHILD, QUARK_PARENT, [2, 3])
    .map(|_| ())
}

#[test]
fn ir_gluon_in_os_quark_forest_matches_independent_component_jets() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_soft_gluon_in_top_self_energy.dot"
    ))?
    .nested_families(GLUON_CHILD, QUARK_PARENT, [2, 3])
    .map(|_| ())
}

#[test]
fn proper_two_loop_os_forest_in_three_loop_graph_keeps_owned_points_and_uv_limits() -> Result<()> {
    test_initialise()?;
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"
    ))?;
    let child = THREE_LOOP_CHILD.node(&oracle.graph)?;
    let parent = THREE_LOOP_PARENT.node(&oracle.graph)?;
    let whole = THREE_LOOP_WHOLE.node(&oracle.graph)?;
    let root = oracle.root();
    assert_eq!(oracle.graph.n_loops(whole.subgraph()), 3);
    assert_eq!(oracle.graph.n_loops(parent.subgraph()), 2);
    assert_eq!(oracle.graph.n_loops(child.subgraph()), 1);
    assert_eq!(oracle.external_edges(&parent), [EdgeIndex(2)]);
    assert_eq!(oracle.external_edges(&child), [EdgeIndex(3)]);

    // The same bounded four-family certificate now acts on a proper component:
    // k/l become hard while the surrounding p loop remains formal and live.
    let complete = oracle.nested_families(THREE_LOOP_CHILD, THREE_LOOP_PARENT, [3, 4])?;
    let mut point_owners = BTreeSet::new();
    complete.visitor(&mut |view| {
        if let AtomView::Fun(point) = view
            && point.get_symbol() == GS.os_mom
            && let Ok(owner) = usize::try_from(point.get(0))
        {
            point_owners.insert(owner);
            assert!(owner == 2 || owner == 3);
            let mass = oracle.graph[EdgeIndex(owner)].mass_atom(&oracle.graph.model);
            assert!(point.get(1) == mass.as_view() || point.get(1) == (-mass).as_view());
        }
        true
    });
    assert_eq!(point_owners, BTreeSet::from([2, 3]));

    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&oracle.graph, &oracle.graph.model, &settings);
    let cograph = oracle.source(&whole, &parent);
    let production_cograph = grow_factor(&ctx, &whole, &parent)?;
    assert_rational_zero(
        "three-loop spectator cograph",
        &oracle.graph,
        &whole,
        &(NumeratorInterface::export(&production_cograph).0 - &cograph),
        SpinReduction::SymbolicFirst,
    )?;
    let cograph_terms = FourDTerm::from_view(cograph.as_view())?;
    assert_eq!(cograph_terms.len(), 1);
    assert_eq!(
        cograph_terms[0]
            .denominators
            .iter()
            .map(|den| den.source_edge)
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([EdgeIndex(2), EdgeIndex(8), EdgeIndex(9)]),
    );
    assert!(!cograph.contains_symbol(GS.os_mom));
    let mut live_outer_momentum = false;
    cograph.visitor(&mut |view| {
        if let AtomView::Fun(momentum) = view
            && momentum.get_symbol() == GS.emr_mom
            && let Ok(owner) = usize::try_from(momentum.get(0))
        {
            assert!(
                owner != 3 && owner != 4,
                "the cograph cannot carry either proper-component hard loop"
            );
            live_outer_momentum |= owner == 2;
        }
        true
    });
    assert!(live_outer_momentum);

    // Check the production attachment boundary on the doubly nested signed
    // family, without distributing a full three-loop numerator. The certified
    // cograph is independent of k/l, so multiplication preserves all four
    // component UV identities above. No outer-p UV subtraction is claimed.
    let seed = Full4dCts::root();
    let local_child = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
    let completed_child = Full4dCts::from_factorized_local(&local_child);
    let nested = uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child)?;
    let completed_parent = Full4dCts::from_factorized_local(&nested);
    assert_eq!(
        grow(&completed_parent, &ctx, &whole, &parent)?,
        finalize(nested.atom().clone().multiply(production_cograph)?)?,
        "the completed proper parent attaches once to the untouched three-loop cograph",
    );
    Ok(())
}

#[test]
fn full_three_loop_os_forests_match_independent_component_jets() -> Result<()> {
    test_initialise()?;
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"
    ))?;
    let whole = THREE_LOOP_WHOLE.node(&oracle.graph)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        renormalization_prescription: crate::uv::RenormalizationPrescriptionSettings {
            massive_power_divergent: ApproximationType::OS,
            ..Default::default()
        },
        ..Default::default()
    };
    let spinneys = oracle.graph.classified_spinneys(
        &oracle.graph.full_filter(),
        &oracle.graph.model,
        &settings,
        &oracle.graph.loop_momentum_basis,
    )?;
    assert!(
        spinneys.iter().all(|spinney| spinney.n_components() <= 1),
        "the nested chain must not have an omitted disconnected forest component"
    );
    let connected = spinneys
        .iter()
        .filter(|spinney| spinney.n_components() == 1)
        .map(|spinney| {
            assert_eq!(spinney.dod, 1);
            assert_eq!(spinney.renormalization_scheme, ApproximationType::OS);
            oracle
                .graph
                .paired_edges(spinney.filter())
                .into_iter()
                .map(usize::from)
                .collect::<BTreeSet<_>>()
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(
        connected,
        [THREE_LOOP_CHILD, THREE_LOOP_PARENT, THREE_LOOP_WHOLE]
            .into_iter()
            .map(|component| component.edges.iter().copied().collect::<BTreeSet<_>>())
            .collect::<BTreeSet<_>>(),
        "the production classifier must yield exactly the three components of the eight forests"
    );
    let prefixes = oracle.three_loop_prefixes()?;
    let ctx = UVCtx::new(&oracle.graph, &oracle.graph.model, &settings);
    // Prove the prefixes in dependency order before reusing their actual
    // factors in the outer independent jet. In particular Gamma-after-gamma
    // acts on the completed child, including its one forest minus sign.
    for (label, output, input, specification) in [
        ("three-loop primitive gamma", 1, 0, THREE_LOOP_CHILD),
        ("three-loop proper Gamma", 2, 0, THREE_LOOP_PARENT),
        (
            "three-loop proper Gamma after gamma",
            3,
            1,
            THREE_LOOP_PARENT,
        ),
    ] {
        let (node, prefix) = &prefixes[output];
        let certified = oracle.assert_operator_branches(
            label,
            node,
            specification,
            &ctx,
            (&prefixes[input].0, &prefixes[input].1),
        )?;
        assert_rational_zero(
            &format!("{label}: one signed completed prefix"),
            &oracle.graph,
            node,
            &(NumeratorInterface::export(prefix.atom()).0 + certified),
            SpinReduction::SymbolicFirst,
        )?;
    }
    let families = oracle.three_loop_families(&prefixes)?;
    for (index, (label, (given, prefix))) in [
        "G",
        "G after gamma (middle omitted)",
        "G after Gamma",
        "G after Gamma after gamma",
    ]
    .into_iter()
    .zip(&prefixes)
    .enumerate()
    {
        let cograph = grow_factor(&ctx, &whole, given)?;
        assert_rational_zero(
            &format!("{label}: independent live cograph"),
            &oracle.graph,
            &whole,
            &(NumeratorInterface::export(&cograph).0 - oracle.source(&whole, given)),
            SpinReduction::SymbolicFirst,
        )?;
        assert_eq!(
            families[2 * index],
            NumeratorInterface::export(&finalize(prefix.atom().clone().multiply(cograph)?)?).0
        );
        let certified = oracle.assert_operator_branches(
            label,
            &whole,
            THREE_LOOP_WHOLE,
            &ctx,
            (given, prefix),
        )?;
        assert_rational_zero(
            &format!("{label}: one outer forest minus"),
            &oracle.graph,
            &whole,
            &(&families[2 * index + 1] + certified),
            SpinReduction::SymbolicFirst,
        )?;
    }
    // The four pairs are the eight ordered forests of gamma < Gamma < G.
    // Each outer branch and each source factor has its own exact certificate;
    // associativity/linearity compose them without a giant tensor difference.
    let paired_sum = families
        .as_chunks::<2>()
        .0
        .iter()
        .fold(Atom::Zero, |sum, pair| sum + &pair[0] + &pair[1]);
    assert_eq!(Atom::add_many(&families), paired_sum);
    Ok(())
}

#[test]
fn three_loop_open_spinor_hard_coefficient_has_exact_finite_certificate() -> Result<()> {
    test_initialise()?;
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"
    ))?;
    let whole = THREE_LOOP_WHOLE.node(&oracle.graph)?;
    let families = oracle.three_loop_families(&oracle.three_loop_prefixes()?)?;
    let parameter = Atom::var(symbol!("os_forest_reference::hard"));
    let coefficients = oracle.hard_coefficients(
        "three-loop innermost hard logarithmic coefficient",
        &Atom::add_many(&families),
        &whole,
        &[(4, function!(GS.emr_mom, 4) / parameter)],
        1,
    )?;
    // This complete coefficient includes the open-spinor denominator class
    // whose canonical Clifford reduction stalled despite its compact input.
    // Certify all classes and all finite components; no coordinate sampling
    // or special-case denominator index enters the proof.
    let logarithmic = coefficients
        .get(&0.into())
        .expect("the completed forest retains a logarithmic coefficient to certify");
    assert_rational_zero(
        "three-loop innermost hard logarithmic coefficient",
        &oracle.graph,
        &whole,
        logarithmic,
        SpinReduction::SymbolicFirst,
    )
}

#[test]
fn full_three_loop_os_forest_removes_generic_hard_limits() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"
    ))?
    .assert_three_loop_hard_strata(&[
        [0, 0, 1],
        [0, 1, 0],
        [1, 0, 0],
        [0, 1, 2],
        [1, 0, 2],
        [1, 2, 0],
        [1, 2, 3],
    ])
}

#[test]
fn full_three_loop_os_forest_removes_correlated_hard_limits() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"
    ))?
    .assert_three_loop_hard_strata(&[
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
        [1, 1, 2],
        [1, 2, 1],
        [1, 2, 2],
    ])
}

#[test]
fn fig10_logarithmic_box_forest_families_cancel_at_the_parent_soft_point() -> Result<()> {
    test_initialise()?;
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
    ))?;
    let child = FIG10_LOGARITHMIC_BOX.node(&oracle.graph)?;
    let parent = GLUON_PARENT.node(&oracle.graph)?;
    let root = oracle.root();
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&oracle.graph, &oracle.graph.model, &settings);
    let seed = Full4dCts::root();
    let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
    let cograph = oracle.source(&parent, &child);
    let denominators = FourDTerm::from_view(cograph.as_view())?;
    assert_eq!(denominators.len(), 1);
    assert_eq!(denominators[0].denominators.len(), 1);
    assert_eq!(denominators[0].denominators[0].source_edge, EdgeIndex(4));
    assert!(denominators[0].denominators[0].mass_squared.is_zero());

    let completed_child = Full4dCts::from_factorized_local(&child_local);
    let completed = child_local
        .atom()
        .clone()
        .multiply(NumeratorInterface::admit(cograph.clone(), None)?.try_into_indexed()?)?;
    let completed_source = NumeratorInterface::export(&completed).0;
    let soft_parameter = symbol!("os_forest_reference::box_soft");
    assert!(
        !oracle
            .scale(
                &oracle.route(&completed_source, &parent),
                &parent,
                soft_parameter,
                false
            )
            .contains_symbol(soft_parameter),
        "the complete logarithmic-box insertion is independent of the parent external momentum",
    );
    let reference_soft = oracle.jet(&completed_source, &parent, false)?;
    assert_rational_zero(
        "the parent S0 reproduces the complete logarithmic-box insertion",
        &oracle.graph,
        &parent,
        &(&reference_soft - &completed_source),
        SpinReduction::SymbolicFirst,
    )?;
    let lmb = compatible_lmb(&ctx, &parent, &child, &completed_child)?;
    let actual_soft = tilde_t_raw(&completed, &ctx, &parent, &child, &lmb)?;
    assert_rational_zero(
        "the logarithmic-box parent soft branch matches its independent jet",
        &oracle.graph,
        &parent,
        &(NumeratorInterface::export(&actual_soft).0 - &reference_soft),
        SpinReduction::SymbolicFirst,
    )?;
    let ordinary = t_raw(&completed, &ctx, &parent, &child, &lmb)?;
    let overlap = t_raw(&actual_soft, &ctx, &parent, &child, &lmb)?;
    assert_rational_zero(
        "the logarithmic-box parent U and US branches agree",
        &oracle.graph,
        &parent,
        &NumeratorInterface::export(&ordinary.try_sub(&overlap)?).0,
        SpinReduction::SymbolicFirst,
    )?;
    let nested = uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child)?;
    assert_rational_zero(
        "the extra logarithmic-box child and nested forest families cancel",
        &oracle.graph,
        &parent,
        &(completed_source + NumeratorInterface::export(nested.atom()).0),
        SpinReduction::SymbolicFirst,
    )
}

#[test]
fn disjoint_os_quark_components_keep_the_live_bridge_momentum() -> Result<()> {
    test_initialise()?;
    let oracle = ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_disjoint_os_top_self_energies.dot"
    ))?;
    let left_spec = DISJOINT_LEFT;
    let right_spec = DISJOINT_RIGHT;
    let left = left_spec.node(&oracle.graph)?;
    let right = right_spec.node(&oracle.graph)?;
    let root = oracle.root();
    let union_spec = Component {
        edges: &[2, 3, 5, 6],
        degree: 2,
        scheme: ApproximationType::MUV,
        on_shell: None,
    };
    let union = union_spec.node(&oracle.graph)?;
    let full_spec = Component {
        edges: &[2, 3, 4, 5, 6],
        degree: 0,
        scheme: ApproximationType::MUV,
        on_shell: None,
    };
    let full = full_spec.node(&oracle.graph)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&oracle.graph, &oracle.graph.model, &settings);
    let seed = Full4dCts::root();
    let left_local = uv_limit(&seed, &ctx, &left, &root, &left, &root)?;
    let right_local = uv_limit(&seed, &ctx, &right, &root, &right, &root)?;
    let product = Local4dCts::from_full_product([
        Full4dCts::from_factorized_local(&left_local),
        Full4dCts::from_factorized_local(&right_local),
    ])?;
    let left_reference = -oracle.operation(&oracle.source(&left, &root), &left, left_spec)?;
    let right_reference = -oracle.operation(&oracle.source(&right, &root), &right, right_spec)?;
    let bridge = oracle.source(&full, &union);
    for (label, component, actual, reference) in [
        (
            "left disjoint OS component",
            &left,
            left_local.atom(),
            &left_reference,
        ),
        (
            "right disjoint OS component",
            &right,
            right_local.atom(),
            &right_reference,
        ),
    ] {
        assert_rational_zero(
            label,
            &oracle.graph,
            component,
            &(NumeratorInterface::export(actual).0 - reference),
            SpinReduction::SymbolicFirst,
        )?;
    }
    // Each open component has its independent exact tensor certificate
    // above. Their unchanged product therefore supplies the product oracle
    // without distributing either numerator over the other component.
    assert_eq!(
        product.atom(),
        &left_local
            .atom()
            .clone()
            .multiply(right_local.atom().clone())?,
        "the disjoint forest is the product of both complete signed components",
    );
    assert_rational_zero(
        "disjoint forest live bridge",
        &oracle.graph,
        &full,
        &(NumeratorInterface::export(&grow_factor(&ctx, &full, &union)?).0 - &bridge),
        SpinReduction::SymbolicFirst,
    )?;
    // Disjoint forests factor into two independently signed one-loop
    // remainders and their live bridge. Check the completed physical product
    // in each hard subspace, not merely the isolated counterterm factors.
    let left_source = oracle.source(&left, &root);
    let right_source = oracle.source(&right, &root);
    let bare = &left_source * &bridge * &right_source;
    let complete = (left_source + NumeratorInterface::export(left_local.atom()).0)
        * &bridge
        * (right_source + NumeratorInterface::export(right_local.atom()).0);
    let parameter = Atom::var(symbol!("os_forest_reference::hard"));
    let left_momentum = function!(GS.emr_mom, 2);
    let right_momentum = function!(GS.emr_mom, 5);
    for (label, momenta, loops, divergent) in [
        (
            "disjoint left hard",
            vec![(2, &left_momentum / &parameter)],
            1,
            true,
        ),
        (
            "disjoint right hard",
            vec![(5, &right_momentum / &parameter)],
            1,
            true,
        ),
        (
            "disjoint both hard",
            vec![
                (2, &left_momentum / &parameter),
                (5, &right_momentum / &parameter),
            ],
            2,
            true,
        ),
        (
            "disjoint correlated hard",
            vec![
                (2, &left_momentum / &parameter),
                (5, &left_momentum / &parameter + &right_momentum),
            ],
            1,
            false,
        ),
    ] {
        oracle.assert_hard_limit(
            label,
            (&complete, divergent.then_some(&bare)),
            &full,
            &momenta,
            loops,
        )?;
    }
    let bridge_terms = FourDTerm::from_view(bridge.as_view())?;
    assert_eq!(bridge_terms.len(), 1);
    assert_eq!(bridge_terms[0].denominators.len(), 1);
    let denominator = &bridge_terms[0].denominators[0];
    assert_eq!(denominator.source_edge, EdgeIndex(4));
    assert!(!denominator.full_expr.contains_symbol(GS.os_mom));
    assert!(denominator.full_expr.contains_symbol(GS.emr_mom));
    let mut point_owners = BTreeSet::new();
    let _ = product.atom().expression().replace_map(|view, _, _| {
        if let AtomView::Fun(point) = view
            && point.get_symbol() == GS.os_mom
            && let Ok(edge) = usize::try_from(point.get(0))
        {
            point_owners.insert(edge);
        }
    });
    assert_eq!(point_owners, BTreeSet::from([0, 4]));
    Ok(())
}

impl ForestOracle {
    fn rainbow_soft_coordinates(
        &self,
        input: &SymbolicTensor<IndexedStructure>,
        external: &[Atom; 4],
        radius: &Atom,
    ) -> SymbolicTensor<IndexedStructure> {
        let bound = input
            .expression()
            .replace_multiple(
                self.masses
                    .iter()
                    .map(|mass| Replacement::new(mass.to_pattern(), Atom::num(12).to_pattern())),
            )
            .replace_map(|view, _, output| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() != 2 {
                    return;
                }
                let Ok(owner) = usize::try_from(momentum.get(0)) else {
                    return;
                };
                if owner > 2 {
                    return;
                }
                let index = momentum.get(1);
                let fixed = external
                    .iter()
                    .enumerate()
                    .map(|(component, value)| {
                        value * function!(GS.delta_vec, GS.cind(component), index)
                    })
                    .sum::<Atom>();
                **output = if owner == 2 {
                    fixed + radius * view.to_owned()
                } else {
                    fixed
                };
            });
        let bound = input
            .with_checked_expression(bound)
            .expect("the fixed soft coordinates retain their native tensor interface");
        finalize(bound).expect("the fixed soft coordinates retain their native tensor interface")
    }

    /// The massive rainbow uses r=e2 and l=e3. Its outer gluon is r-q;
    /// approach it with r=q+rho*K while leaving the inner loop completely
    /// symbolic. Only the required Laurent coefficients are constructed.
    fn rainbow_soft_coefficients(
        &self,
        atom: &Atom,
        node: &OwnedForestNode,
        external: &[Atom; 4],
        through_power: i64,
    ) -> Result<BTreeMap<symbolica::domains::rational::Rational, Atom>> {
        let variable = symbol!("os_forest_reference::rainbow_soft");
        let radius = Atom::var(variable);
        let physical = self.route(&Self::physical_limit(atom), node);
        let mut coefficients = BTreeMap::new();
        for term in FourDTerm::from_view(physical.as_view())? {
            let numerator = SymbolicTensor::from(term.numerator)
                .undo_dots()?
                .try_into_indexed()?;
            let scaled = self.rainbow_soft_coordinates(&numerator, external, &radius);
            // This is the explicit independent finite-component jet boundary.
            let mut scaled = NumeratorInterface::export(&scaled).0;
            // The finite jet differentiates vector coefficients, not a compact
            // gamma function's vector metadata. Reuse the same selective native
            // materialization as the independent component jets above.
            scaled = scaled
                .parse_to_symbolic_net::<Aind>(&ParseSettings {
                    shorthand_parsing: ShorthandParsing::Expand {
                        schoonschip: SchoonschipExpansionMode::full(),
                        trace: false,
                        chain: false,
                    },
                    ..Default::default()
                })?
                .simple_execute::<()>()?;
            for denominator in term.denominators {
                let polynomial = self.rainbow_soft_coordinates(
                    &NumeratorInterface::admit(denominator.full_expr.clone(), None)?
                        .try_into_indexed()?,
                    external,
                    &radius,
                );
                let polynomial =
                    scalar_component(&NumeratorInterface::export(&polynomial).0)?.expand();
                let leading = polynomial
                    .series(variable, Atom::Zero, 2)?
                    .get_trailing_exponent();
                assert!(
                    (symbolica::domains::rational::Rational::from(0)
                        ..=symbolica::domains::rational::Rational::from(2))
                        .contains(&leading),
                    "a soft denominator is at most quadratic"
                );
                let polynomial = (polynomial / radius.pow(leading.clone())).expand();
                assert!(!polynomial.replace(variable).with(Atom::Zero).is_zero());
                scaled /= radius.pow(leading)
                    * GS.den(
                        usize::from(denominator.source_edge),
                        denominator.momentum,
                        denominator.mass_squared,
                        polynomial,
                    );
            }
            // rho^-4 is logarithmic with the four-dimensional radial measure;
            // rho^-3 is the first absolutely integrable candidate.
            let jet = scaled.series_preserving_factors(
                variable,
                Atom::Zero.as_view(),
                through_power,
                &[],
            )?;
            assert!(
                !jet.contains_symbol(Symbol::DERIVATIVE),
                "the rainbow soft jet must contain explicit coefficient derivatives"
            );
            for (power, coefficient) in jet.series(variable, Atom::Zero, through_power)?.terms() {
                *coefficients.entry(power).or_insert(Atom::Zero) += coefficient;
            }
        }
        Ok(coefficients)
    }

    fn assert_rainbow_soft_nonzero(&self, atom: &Atom, label: &str) -> Result<()> {
        // A single exact rational point certifies non-identity; it never
        // supplies a zero certificate. K=(2,1,0,0), l=(7,2,3,1) are away from
        // every soft-coefficient pole in these massive fixtures.
        let sampled = atom
            .replace_map(|view, _, output| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() != 2 {
                    return;
                }
                let Ok(owner) = usize::try_from(momentum.get(0)) else {
                    return;
                };
                let AtomView::Fun(index) = momentum.get(1) else {
                    return;
                };
                let Ok(component) = usize::try_from(index.get(0)) else {
                    return;
                };
                if component >= 4 {
                    return;
                }
                let point = if owner == 2 {
                    [2, 1, 0, 0]
                } else {
                    [7, 2, 3, 1]
                };
                **output = Atom::num(point[component]);
            })
            .replace(function!(symbol!(
                "os_forest_reference::frozen_child_localizer"
            )))
            .with(Atom::num(20) / (Atom::var(GS.pi).pow(2) * Atom::num(414).pow(2)))
            .replace(GS.mu_r_sq)
            .with(400)
            .replace(GS.m_uv_vacuum)
            .with(20)
            .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
            .with(W_.d_);
        // Keep the small physical matrix product factorized. Symbolic gamma
        // collection at a nonaxial fixed point manufactures large scalar sums.
        let tensor = crate::uv::overlap_numerators::NumeratorJetCache::default()
            .canonicalize(&sampled.undo_dots::<Aind>()?)?;
        let mut nonzero = false;
        tensor.key.visitor(&mut |view| {
            let AtomView::Fun(entry) = view else {
                return true;
            };
            if entry.get_symbol() != symbol!("gammalooprs::uv::overlap_entry") {
                return true;
            }
            let value = entry
                .get(1)
                .to_owned()
                .replace(idenso::color_cas!(2, idenso::cof!(3)))
                .with(Atom::num(4) / Atom::num(3))
                .replace(idenso::color_idx!(2, idenso::cof!(3)))
                .with(Atom::one() / Atom::num(2))
                .replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() == 2
                        && let Ok(owner) = usize::try_from(momentum.get(0))
                        && let AtomView::Fun(index) = momentum.get(1)
                        && let Ok(component) = usize::try_from(index.get(0))
                        && component < 4
                    {
                        let point = if owner == 2 {
                            [2, 1, 0, 0]
                        } else {
                            [7, 2, 3, 1]
                        };
                        **output = Atom::num(point[component]);
                    } else if let AtomView::Var(variable) = view
                        && variable.get_symbol() != GS.pi
                    {
                        **output = Atom::one();
                    }
                });
            // At this exact rational point, only an exactly nonzero scalar
            // factorization certifies the control. Zero iterations take no
            // numerical samples; an unresolved sum remains inconclusive.
            nonzero |= value
                .collect_compact_factors()
                .check_zero_with(
                    &symbolica::atom::ZeroTestSettings::default()
                        .finite_field_samples(0)
                        .numerical_samples(0)
                        .exact_fallback(false),
                )
                .is_false();
            !nonzero
        });
        assert!(
            nonzero,
            "{label}: an exact nonzero tensor component is required"
        );
        Ok(())
    }

    fn assert_boosted_rainbow_soft_limit(&mut self, velocity: [&str; 4]) -> Result<()> {
        let external_velocity = UVgenerationSettings {
            on_shell_subtraction_velocity: velocity.map(str::to_owned),
            ..Default::default()
        }
        .on_shell_velocity()?;
        let external = external_velocity
            .clone()
            .map(|component| Atom::num(12) * component);
        let child = QUARK_CHILD.node(&self.graph)?;
        let parent = QUARK_PARENT.node(&self.graph)?;
        let root = self.root();
        let seed = Full4dCts::root();
        let matched_boost = velocity != ["1", "0", "0", "0"];
        for (label, subtraction_velocity) in [("rest", ["1", "0", "0", "0"]), ("matched", velocity)]
            .into_iter()
            .filter(|(label, _)| matched_boost || *label == "rest")
        {
            let settings = UVgenerationSettings {
                generate_integrated: true,
                on_shell_subtraction_velocity: subtraction_velocity.map(str::to_owned),
                ..Default::default()
            };
            self.velocity = settings.on_shell_velocity()?;
            let ctx = UVCtx::new(&self.graph, &self.graph.model, &settings);
            let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
            let vakint_settings = vakint::VakintSettings {
                number_of_terms_in_epsilon_expansion: 3,
                evaluation_order: vakint::EvaluationOrder(vec![
                    vakint::EvaluationMethod::AlphaLoop(vakint::AlphaLoopOptions {
                        susbstitute_masters: false,
                    }),
                ]),
                ..settings.vakint.true_settings()
            };
            let child_integrated = integrated::Integrated::new(
                crate::utils::vakint()?,
                &vakint_settings,
            )
            .run(&child_local, &ctx, &child, &root, &child, &root)?;
            let completed_child = Full4dCts::recursion_input(
                &child_local,
                &child_integrated,
                ApproximationType::OS,
                false,
                child.lmb(),
            )?;
            let parent_local = uv_limit(&seed, &ctx, &parent, &root, &parent, &root)?;
            let nested_local = uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child)?;
            let bare = grow(&seed, &ctx, &parent, &root)?;
            let local_child = Full4dCts::from_factorized_local(&child_local);
            let insertion = grow(&local_child, &ctx, &parent, &child)?;
            let finite_insertion =
                grow(&completed_child, &ctx, &parent, &child)?.try_sub(&insertion)?;
            let frozen = function!(symbol!("os_forest_reference::frozen_child_localizer"));
            assert!(nested_local.recursive_completion().is_empty());
            assert_eq!(
                nested_local
                    .active_sectors()
                    .iter()
                    .filter(|sector| !sector.frozen_lmbs().is_empty())
                    .count(),
                1,
                "the complete parent operates on exactly one finite child prefix"
            );
            let nested = nested_local.active_sectors().iter().try_fold(
                bare.with_checked_expression(Atom::Zero)?,
                |sum, sector| -> Result<_> {
                    let weighted = if sector.frozen_lmbs().is_empty() {
                        sector.atom.clone()
                    } else {
                        assert_eq!(sector.frozen_lmbs(), std::slice::from_ref(child.lmb()));
                        sector
                            .atom
                            .with_checked_expression(&frozen * sector.atom.expression())?
                    };
                    Ok(sum.try_add(&weighted)?)
                },
            )?;
            let frozen_finite_insertion = finite_insertion
                .with_checked_expression(&frozen * finite_insertion.expression())?;
            let complete = bare
                .try_add(&insertion)?
                .try_add(parent_local.atom())?
                .try_add(&nested)?
                .try_add(&frozen_finite_insertion)?;

            // Independent source and bounded U/O/UO jets certify all four
            // local forest families. The analytic child is carried separately
            // with its frozen localization owner through the actual recursion.
            let child_reference =
                -self.operation(&self.source(&child, &root), &child, QUARK_CHILD)?;
            assert_rational_zero(
                "boosted rainbow child operation",
                &self.graph,
                &child,
                &(NumeratorInterface::export(child_local.atom()).0 - child_reference),
                SpinReduction::SymbolicFirst,
            )?;
            let parent_operator = self.assert_operator_branches(
                "boosted rainbow parent",
                &parent,
                QUARK_PARENT,
                &ctx,
                (&root, &seed),
            )?;
            assert_rational_zero(
                "rainbow parent has exactly one forest minus sign",
                &self.graph,
                &parent,
                &(NumeratorInterface::export(parent_local.atom()).0 + parent_operator),
                SpinReduction::SymbolicFirst,
            )?;
            let nested_operator = self.assert_operator_branches(
                "boosted rainbow nested parent",
                &parent,
                QUARK_PARENT,
                &ctx,
                (&child, &local_child),
            )?;

            let nested_local_family = nested_local
                .active_sectors()
                .iter()
                .filter(|sector| sector.frozen_lmbs().is_empty())
                .try_fold(
                    bare.with_checked_expression(Atom::Zero)?,
                    |sum, sector| -> Result<_> { Ok(sum.try_add(&sector.atom)?) },
                )?;
            assert_rational_zero(
                "rainbow nested parent has exactly one additional forest minus sign",
                &self.graph,
                &parent,
                &(NumeratorInterface::export(&nested_local_family).0 + nested_operator),
                SpinReduction::SymbolicFirst,
            )?;

            // Independent leading-pole formula: the cograph has two massive
            // quark lines r and one gluon r-q. At r=q+rho*K their leading
            // polynomials are 2*rho*q.K, 2*rho*q.K and rho^2*K^2. This
            // constructs gamma.(qslash+m).R_child.(qslash+m).gamma without
            // using the Laurent extraction or any production Taylor routine.
            let cograph = self.route(
                &Self::physical_limit(&self.source(&parent, &child)),
                &parent,
            );
            let cograph_terms = FourDTerm::from_view(cograph.as_view())?;
            assert_eq!(cograph_terms.len(), 1);
            assert_eq!(
                cograph_terms[0]
                    .denominators
                    .iter()
                    .map(|den| usize::from(den.source_edge))
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([2, 5, 6])
            );
            let cograph_term = cograph_terms.into_iter().next().unwrap();
            let cograph_numerator = SymbolicTensor::from(cograph_term.numerator)
                .undo_dots()?
                .try_into_indexed()?;
            let cograph_numerator =
                self.rainbow_soft_coordinates(&cograph_numerator, &external, &Atom::Zero);
            let child_finite = child_integrated.physical_finite_counterterm()?;
            let frozen_child_finite =
                child_finite.with_checked_expression(&frozen * child_finite.expression())?;
            let child_completed = NumeratorInterface::admit(self.source(&child, &root), None)?
                .try_into_indexed()?
                .try_add(child_local.atom())?
                .try_add(&frozen_child_finite)?;
            let child_completed = child_completed.with_checked_expression(
                self.route(&Self::physical_limit(child_completed.expression()), &parent),
            )?;
            let child_at_point =
                self.rainbow_soft_coordinates(&child_completed, &external, &Atom::Zero);
            let hard = function!(GS.emr_mom, 2);
            let k_squared = (0..4)
                .map(|index| {
                    Atom::num(if index == 0 { 1 } else { -1 })
                        * GS.emr_mom(EdgeIndex(2), GS.cind(index)).pow(2)
                })
                .sum::<Atom>();
            let twice_q_dot_k = Atom::num(2)
                * external
                    .iter()
                    .enumerate()
                    .map(|(index, value)| {
                        Atom::num(if index == 0 { 1 } else { -1 })
                            * value
                            * GS.emr_mom(EdgeIndex(2), GS.cind(index))
                    })
                    .sum::<Atom>();
            let expected_source_coefficient = NumeratorInterface::export(&cograph_numerator).0
                * NumeratorInterface::export(&child_at_point).0
                / (GS.den(6, &hard, Atom::Zero, k_squared)
                    * GS.den(2, &hard, Atom::num(144), &twice_q_dot_k)
                    * GS.den(5, &hard, Atom::num(144), twice_q_dot_k));
            let source = bare
                .try_add(&insertion)?
                .try_add(&frozen_finite_insertion)?;
            let source_coefficients = self.rainbow_soft_coefficients(
                &NumeratorInterface::export(&source).0,
                &parent,
                &external,
                -3,
            )?;
            assert_rational_zero(
                "rainbow leading pole equals independent doubled-quark cograph formula",
                &self.graph,
                &parent,
                &(source_coefficients.get(&(-4).into()).unwrap() - &expected_source_coefficient),
                SpinReduction::FiniteComponents,
            )?;
            if label == "matched" || !matched_boost {
                assert_rational_zero(
                    "matched completed child removes the cograph logarithmic coefficient",
                    &self.graph,
                    &parent,
                    &expected_source_coefficient,
                    SpinReduction::FiniteComponents,
                )?;
            }

            // Overall analytic addbacks have no live loop momentum in their
            // coefficient. Their actual normalized spatial kernel is regular
            // on this stratum, so an arbitrary overall coefficient cannot
            // contribute to either negative power tested here. The child
            // localizer above stays frozen during the parent operation.
            let overall_localizer = GS
                .localizing_integrand(parent.lmb())
                .replace(function!(GS.localizing_integrand, W_.a_))
                .with(W_.a_);
            for power in self
                .rainbow_soft_coefficients(&overall_localizer, &parent, &external, 0)?
                .into_keys()
            {
                assert!(
                    power >= 0,
                    "the overall analytic localizer is regular, not merely free of logarithmic soft poles"
                );
            }

            let coefficients = self.rainbow_soft_coefficients(
                &NumeratorInterface::export(&complete).0,
                &parent,
                &external,
                -3,
            )?;
            for (power, coefficient) in &coefficients {
                assert!(*power >= -4, "no worse than the doubled-fermion soft pole");
                if *power <= -4 {
                    if matched_boost && label == "rest" {
                        self.assert_rainbow_soft_nonzero(
                            coefficient,
                            "boosted on-shell rainbow with rest-frame subtraction",
                        )?;
                    } else {
                        assert_rational_zero(
                            "matched complete rainbow logarithmic soft coefficient",
                            &self.graph,
                            &parent,
                            coefficient,
                            SpinReduction::FiniteComponents,
                        )?;
                    }
                }
            }
            if label == "matched" || !matched_boost {
                self.assert_rainbow_soft_nonzero(
                    coefficients
                        .get(&(-3).into())
                        .expect("the first integrable soft power is retained"),
                    "matched complete rainbow has a nonzero rho^-3 coefficient",
                )?;
            }
            let logarithmic = coefficients
                .get(&(-4).into())
                .expect("the soft jet retains its logarithmic coefficient");
            if matched_boost && label == "rest" {
                // At a different fixed subtraction point, parent U and both
                // parent O branches lack this outer-gluon soft denominator.
                // The complete rho^-4 coefficient must therefore equal the
                // actual cograph insertion of the completed child exactly.
                assert_rational_zero(
                    "rest prescription retains precisely the completed-child soft coefficient",
                    &self.graph,
                    &parent,
                    &(logarithmic - &expected_source_coefficient),
                    SpinReduction::FiniteComponents,
                )?;
            }
            debug_tags!(#uv, #local, #soft, #summary;
                stage = "boosted_rainbow_soft_certificate", label,
                external_velocity = ?external_velocity,
                logarithmic_coefficient_zero = !(matched_boost && label == "rest"),
                coefficient_powers = ?coefficients.keys().collect::<Vec<_>>(),
                "Certified the complete rainbow outer-soft Laurent coefficients"
            );
        }
        Ok(())
    }
}

#[test]
fn complete_two_loop_rainbow_soft_limit_at_rest_matches_fixed_point() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
    ))?
    .assert_boosted_rainbow_soft_limit(["1", "0", "0", "0"])
}

#[test]
fn complete_two_loop_rainbow_soft_limit_at_axial_boost_compares_fixed_points() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
    ))?
    .assert_boosted_rainbow_soft_limit(["5/4", "3/4", "0", "0"])
}

#[test]
fn complete_two_loop_rainbow_soft_limit_at_nonaxial_boost_compares_fixed_points() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
    ))?
    .assert_boosted_rainbow_soft_limit(["5/4", "1/4", "1/2", "1/2"])
}
