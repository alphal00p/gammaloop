//! Independent bounded-jet references for small physical OS forests.
//! Source numerators remain factorized; only scalar denominator polynomials
//! and the Taylor coefficients demanded by a component are materialized.

use super::on_shell_tests::{assert_rational_zero, scalar_component};
use super::*;
use crate::{
    graph::parse::IntoGraph,
    initialisation::test_initialise,
    uv::{Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode},
};
use idenso::shorthands::UndoShorthands;
use linnet::half_edge::{involution::Hedge, subgraph::ModifySubSet};

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
}

impl ForestOracle {
    fn new(source: &str) -> Result<Self> {
        let graph: Graph = source.into_graph(&crate::utils::load_generic_model("sm"))?;
        let masses = graph
            .underlying
            .iter_edges()
            .map(|(_, edge, _)| graph[edge].mass_atom())
            .filter(|mass| !mass.is_zero())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(Self { graph, masses })
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
            .to_d_dim(4)
            .color_simplify()
            .get_single_atom()
            .expect("the physical source is one factorized graph numerator");
        self.route(&(numerator / self.graph.denominator(&reduced, |_| 1)), node)
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
            let numerator = self.scale(
                &term.numerator.undo_dots::<Aind>()?,
                node,
                variable,
                ultraviolet,
            );
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
        let mass = self.graph[EdgeIndex(incoming)].mass_atom();
        let port = self.graph.underlying[Hedge(outgoing)]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom();
        let external = self.external_edges(node);
        let routed = self.route(atom, node);
        [1, -1].into_iter().fold(Atom::Zero, |sum, sigma| {
            let point = function!(GS.os_mom, incoming, Atom::num(sigma) * &mass);
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
            let gamma =
                idenso::gamma0!(&dummy, &port) * evaluated.replace(port.to_pattern()).with(dummy);
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

    fn assert_primitive_branches(
        &self,
        label: &str,
        node: &OwnedForestNode,
        component: Component,
        ctx: &UVCtx<'_>,
        root: &OwnedForestNode,
    ) -> Result<Atom> {
        let seed = Full4dCts::root();
        let actual_source = grow(&seed, ctx, node, root)?;
        let source = self.source(node, root);
        let lmb = compatible_lmb(ctx, node, root, &seed)?;
        assert_rational_zero(
            &format!("{label} source"),
            &self.graph,
            node,
            &(&actual_source - &source),
        )?;
        let actual_u = t_raw(&actual_source, ctx, node, root, &lmb)?;
        let reference_u = self.jet(&source, node, true)?;
        assert_rational_zero(
            &format!("{label} U branch"),
            &self.graph,
            node,
            &(&actual_u - reference_u),
        )?;
        let (actual_refinement, reference_refinement) = match component.on_shell {
            Some(boundary) => {
                let [positive, negative] =
                    self.graph
                        .on_shell_branches(&actual_source, node.subgraph(), &lmb)?;
                (positive + negative, self.on_shell(&source, node, boundary))
            }
            None => (
                tilde_t_raw(&actual_source, ctx, node, root, &lmb)?,
                self.jet(&source, node, false)?,
            ),
        };
        assert_rational_zero(
            &format!("{label} refinement branch"),
            &self.graph,
            node,
            &(&actual_refinement - &reference_refinement),
        )?;
        let actual_overlap = t_raw(
            &GS.resolve_os_momenta(&actual_refinement),
            ctx,
            node,
            root,
            &lmb,
        )?;
        let reference_overlap = self.jet(&reference_refinement, node, true)?;
        assert_rational_zero(
            &format!("{label} overlap branch"),
            &self.graph,
            node,
            &(&actual_overlap - reference_overlap),
        )?;
        Ok(finalize(actual_u + actual_refinement - actual_overlap))
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
            .with(Atom::Zero);
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
        for term in FourDTerm::from_view(physical.as_view())? {
            let mut scaled = term
                .numerator
                .undo_dots::<Aind>()?
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
        }
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
        let (complete, bare_control) = sources;
        for (power, coefficient) in
            self.hard_coefficients(label, complete, node, momenta, hard_loops)?
        {
            assert_rational_zero(
                &format!("{label}: measure-weighted inverse-hard coefficient {power}"),
                &self.graph,
                node,
                &coefficient,
            )?;
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
            let tensor = crate::uv::overlap_numerators::NumeratorJetCache::default()
                .canonicalize(&terms[0].numerator.undo_dots::<Aind>()?)?;
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
        Ok(())
    }

    fn nested_families(&self, child_spec: Component, parent_spec: Component) -> Result<()> {
        let child = child_spec.node(&self.graph)?;
        let parent = parent_spec.node(&self.graph)?;
        let root = self.root();
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&self.graph, &settings);
        let seed = Full4dCts::root();
        let bare_child = self.source(&child, &root);
        let child_reference = -self.operation(&bare_child, &child, child_spec)?;
        let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
        let certified_child =
            self.assert_primitive_branches("child", &child, child_spec, &ctx, &root)?;
        // Each branch already equals its independent bounded-jet reference.
        // Certify the one forest minus sign against their exact composition.
        assert_rational_zero(
            "nested forest child is minus its independently certified branch sum",
            &self.graph,
            &child,
            &(child_local.atom() + certified_child),
        )?;
        let completed_child = Full4dCts::from_factorized_local(&child_local);
        if child_spec.scheme == ApproximationType::OS && parent_spec.scheme == ApproximationType::IR
        {
            // The completed child depends on the parent's loop p and on its
            // own fixed rest-frame points, but not on the parent's external q.
            // Check this boundary before attaching the three-line cograph.
            let lmb = compatible_lmb(&ctx, &parent, &child, &completed_child)?;
            let actual_soft = tilde_t_raw(&child_local.atom(), &ctx, &parent, &child, &lmb)?;
            let reference_soft = self.jet(&child_reference, &parent, false)?;
            assert_rational_zero(
                "the parent soft jet leaves the independent OS child fixed",
                &self.graph,
                &parent,
                &(&reference_soft - &child_reference),
            )?;
            assert_rational_zero(
                "parent S of the actual completed OS child matches the independent jet",
                &self.graph,
                &parent,
                &(actual_soft - reference_soft),
            )?;
        }

        let bare_parent = self.source(&parent, &root);
        let cograph = self.source(&parent, &child);
        let nested_reference =
            -self.operation(&(&child_reference * &cograph), &parent, parent_spec)?;
        let parent_local = uv_limit(&seed, &ctx, &parent, &root, &parent, &root)?;
        let certified_parent =
            self.assert_primitive_branches("parent", &parent, parent_spec, &ctx, &root)?;
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
            &(&production_cograph - &cograph),
        )?;
        let child_insertion = grow(&completed_child, &ctx, &parent, &child)?;
        assert_eq!(
            child_insertion,
            (child_local.atom() * production_cograph).simplify_metrics(),
            "the child forest family inserts the complete signed child exactly once",
        );
        assert_rational_zero(
            "nested forest parent is minus its independently certified branch sum",
            &self.graph,
            &parent,
            &(parent_local.atom() + certified_parent),
        )?;
        assert_rational_zero(
            "nested forest completed-child family",
            &self.graph,
            &parent,
            &(nested_local.atom() - &nested_reference),
        )?;
        if child_spec.scheme == ApproximationType::OS
            && parent_spec.scheme == ApproximationType::IR
            && self
                .graph
                .has_qcd_self_energy_insertion_ward_identity(parent.subgraph(), child.subgraph())
        {
            let completed_source = &child_reference * &cograph;
            let ordinary = self.jet(&completed_source, &parent, true)?;
            let soft = self.jet(&completed_source, &parent, false)?;
            let overlap = self.jet(&soft, &parent, true)?;
            let analytic = nested_local
                .active_sectors()
                .iter()
                .fold(Atom::Zero, |sum, sector| sum + sector.integration_atom());
            assert_rational_zero(
                "signed nested analytic input is minus U plus US of the completed child",
                &self.graph,
                &parent,
                &(analytic + ordinary - overlap),
            )?;
        }
        let production_source = grow(&seed, &ctx, &parent, &root)?;
        let complete =
            &production_source + &child_insertion + parent_local.atom() + nested_local.atom();
        let signed_families = [
            production_source,
            child_insertion,
            parent_local.atom().clone(),
            nested_local.atom().clone(),
        ];
        assert_eq!(
            complete,
            Atom::add_many(&signed_families),
            "the complete forest adds each certified signed family once without a second depth sign",
        );
        let parameter = Atom::var(symbol!("os_forest_reference::hard"));
        let p = function!(GS.emr_mom, 2);
        let k = function!(GS.emr_mom, 3);
        for (label, momenta, loops) in [
            ("hard child loop", vec![(3, &k / &parameter)], 1),
            ("hard parent loop", vec![(2, &p / &parameter)], 1),
            (
                "both generic hard loops",
                vec![(2, &p / &parameter), (3, &k / &parameter)],
                2,
            ),
            (
                "correlated hard k=p+r",
                vec![(2, &p / &parameter), (3, &p / &parameter + &k)],
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
            )?;
        }
        if child_spec.scheme == ApproximationType::IR {
            assert_rational_zero(
                "the IR gluon child remainder has no constant or linear soft jet",
                &self.graph,
                &child,
                &self.jet(
                    &Self::physical_limit(&(bare_child + child_local.atom())),
                    &child,
                    false,
                )?,
            )?;
        }
        // The bare source, signed child insertion, signed parent and signed
        // completed-child parent have separate exact certificates above.
        // Linearity composes those certificates into the full four-family
        // identity without re-expanding a difference of complete forests.
        Ok(())
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
    let ctx = UVCtx::new(&oracle.graph, &settings);
    let seed = Full4dCts::root();
    let source = grow(&seed, &ctx, &child, &root)?;
    let lmb = compatible_lmb(&ctx, &child, &root, &seed)?;
    let [positive, negative] = oracle
        .graph
        .on_shell_branches(&source, child.subgraph(), &lmb)?;
    let mass = oracle.graph[EdgeIndex(3)].mass_atom();
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
            &GS.resolve_os_momenta(&GS.erase_uv_momentum_provenance(&input)),
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
        let numerator = terms[0]
            .numerator
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
    let ctx = UVCtx::new(&oracle.graph, &settings);
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
    let inherited = owner_two(&source)
        .into_iter()
        .filter(|(role, _, _)| *role == UvMomentumProvenanceRole::DenominatorDerivedSoft)
        .collect::<Vec<_>>();
    assert!(
        !inherited.is_empty(),
        "the physical child must carry its soft crown owner e2"
    );
    let branches = oracle
        .graph
        .on_shell_branches(&source, parent.subgraph(), &lmb)?;
    for branch in &branches {
        let records = owner_two(branch);
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
        &GS.erase_uv_momentum_provenance(&source),
        &parent,
        QUARK_PARENT.on_shell.unwrap(),
    );
    assert_rational_zero(
        "parent O provenance promotion leaves the independent owned momentum operation unchanged",
        &oracle.graph,
        &parent,
        &(branches.into_iter().sum::<Atom>() - reference),
    )
}

#[test]
fn os_quark_in_ir_gluon_forest_matches_independent_component_jets() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
    ))?
    .nested_families(QUARK_CHILD, GLUON_PARENT)
}

#[test]
fn os_quark_in_os_quark_forest_matches_independent_component_jets() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
    ))?
    .nested_families(QUARK_CHILD, QUARK_PARENT)
}

#[test]
fn ir_gluon_in_os_quark_forest_matches_independent_component_jets() -> Result<()> {
    test_initialise()?;
    ForestOracle::new(include_str!(
        "../../../../../tests/resources/graphs/local_soft_gluon_in_top_self_energy.dot"
    ))?
    .nested_families(GLUON_CHILD, QUARK_PARENT)
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
    let ctx = UVCtx::new(&oracle.graph, &settings);
    let seed = Full4dCts::root();
    let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
    let cograph = oracle.source(&parent, &child);
    let denominators = FourDTerm::from_view(cograph.as_view())?;
    assert_eq!(denominators.len(), 1);
    assert_eq!(denominators[0].denominators.len(), 1);
    assert_eq!(denominators[0].denominators[0].source_edge, EdgeIndex(4));
    assert!(denominators[0].denominators[0].mass_squared.is_zero());

    let completed_child = Full4dCts::from_factorized_local(&child_local);
    let completed = child_local.atom() * &cograph;
    let soft_parameter = symbol!("os_forest_reference::box_soft");
    assert!(
        !oracle
            .scale(
                &oracle.route(&completed, &parent),
                &parent,
                soft_parameter,
                false
            )
            .contains_symbol(soft_parameter),
        "the complete logarithmic-box insertion is independent of the parent external momentum",
    );
    let reference_soft = oracle.jet(&completed, &parent, false)?;
    assert_rational_zero(
        "the parent S0 reproduces the complete logarithmic-box insertion",
        &oracle.graph,
        &parent,
        &(&reference_soft - &completed),
    )?;
    let lmb = compatible_lmb(&ctx, &parent, &child, &completed_child)?;
    let actual_soft = tilde_t_raw(&completed, &ctx, &parent, &child, &lmb)?;
    assert_rational_zero(
        "the logarithmic-box parent soft branch matches its independent jet",
        &oracle.graph,
        &parent,
        &(&actual_soft - &reference_soft),
    )?;
    let ordinary = t_raw(&completed, &ctx, &parent, &child, &lmb)?;
    let overlap = t_raw(&actual_soft, &ctx, &parent, &child, &lmb)?;
    assert_rational_zero(
        "the logarithmic-box parent U and US branches agree",
        &oracle.graph,
        &parent,
        &(ordinary - overlap),
    )?;
    let nested = uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child)?;
    assert_rational_zero(
        "the extra logarithmic-box child and nested forest families cancel",
        &oracle.graph,
        &parent,
        &(completed + nested.atom()),
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
    let ctx = UVCtx::new(&oracle.graph, &settings);
    let seed = Full4dCts::root();
    let left_local = uv_limit(&seed, &ctx, &left, &root, &left, &root)?;
    let right_local = uv_limit(&seed, &ctx, &right, &root, &right, &root)?;
    let product = Local4dCts::from_full_product([
        Full4dCts::from_factorized_local(&left_local),
        Full4dCts::from_factorized_local(&right_local),
    ]);
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
        assert_rational_zero(label, &oracle.graph, component, &(actual - reference))?;
    }
    // Each open component has its independent exact tensor certificate
    // above. Their unchanged product therefore supplies the product oracle
    // without distributing either numerator over the other component.
    assert_eq!(
        product.atom(),
        &(left_local.atom() * right_local.atom()),
        "the disjoint forest is the product of both complete signed components",
    );
    assert_rational_zero(
        "disjoint forest live bridge",
        &oracle.graph,
        &full,
        &(grow_factor(&ctx, &full, &union)? - &bridge),
    )?;
    // Disjoint forests factor into two independently signed one-loop
    // remainders and their live bridge. Check the completed physical product
    // in each hard subspace, not merely the isolated counterterm factors.
    let left_source = oracle.source(&left, &root);
    let right_source = oracle.source(&right, &root);
    let bare = &left_source * &bridge * &right_source;
    let complete =
        (left_source + left_local.atom()) * &bridge * (right_source + right_local.atom());
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
    let _ = product.atom().replace_map(|view, _, _| {
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
