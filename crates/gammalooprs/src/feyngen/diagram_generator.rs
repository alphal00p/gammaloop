use idenso::IndexTooling;
use idenso::color::{CS, ColorSimplifier};
use idenso::dirac::{AGS, GammaSimplifier};

use linnet::permutation::Permutation;

use spenso::network::library::LibraryTensor;
use spenso::network::library::symbolic::{ExplicitKey, TensorLibrary};
use spenso::network::parsing::{ParseSettings, ShorthandParsing};
use spenso::network::{MinResultRank, Sequential};

// use spenso::network::Network;

// use symbolica_utils::AtomPrintExt;
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::{PermutedStructure, TensorStructure};
use spenso::tensors::data::DataTensor;
use spenso::tensors::parametric::ParamTensor;
use spenso_hep_lib::{gamma_data_weyl, gamma5_weyl_data, proj_m_data_weyl, proj_p_data_weyl};
use std::collections::{HashSet, VecDeque};

use std::ops::Deref;
use symbolica::atom::AtomView;
use symbolica::coefficient::Coefficient;
use symbolica::domains::algebraic_number::AlgebraicExtension;
use symbolica::domains::finite_field::PrimeIteratorU64;
use symbolica::domains::float::Complex as SymbolicaComplex;
use symbolica::function;
use symbolica::id::Replacement;

use tracing::{event_enabled, instrument};

use ahash::AHashMap;
use ahash::AHashSet;
use ahash::HashMap;
use smartstring::SmartString;
use symbolica::atom::AtomCore;
use symbolica::{parse, symbol};
use tracing::debug;

use super::FeynGenError;
use super::NumeratorAwareGraphGroupingOption;
use crate::graph::parse::string_utils::ToOrderedSimple;
use crate::graph::{Graph, LMBext};
use crate::is_interrupted;
use crate::model::ArcParticle;
use crate::model::ArcVertexRule;
use crate::numerator::ParamParsingNet;
use crate::numerator::aind::Aind;
use crate::numerator::graph::ReversibleEdge;
use crate::processes::ProcessDefinition;
use crate::settings::GlobalSettings;
use crate::utils::symbolica_ext::{COMPLEXRATPOLYFIELD, LOGPRINTOPTS, Q_I};
use crate::utils::{GS, PARAM_FUN_LIB, W_};
use crate::{feyngen::GenerationType, model::Model};
use eyre::eyre;
use symbolica_utils::{AtomFloatExt, PrimeGenerate};

pub(crate) type NumeratorSample = (
    Vec<Replacement>,
    TensorLibrary<ParamTensor<ExplicitKey<Aind>>, Aind>,
);
type RationalPoly = symbolica::poly::polynomial::MultivariatePolynomial<
    AlgebraicExtension<
        symbolica::domains::rational::FractionField<symbolica::domains::integer::IntegerRing>,
    >,
>;
type SampleEvaluationsAsPolynomial = (Vec<RationalPoly>, Vec<bool>);

use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeData, Flow, Orientation};
use symbolica::{atom::Atom, graph::Graph as SymbolicaGraph};

const ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL: bool = true;
const EXPAND_NUMERICAL_SAMPLES_BEFORE_COMPARISON: bool = false;
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeColorWithVertexRule {
    pub external_tag: i32,
    pub vertex_rule: ArcVertexRule,
}

impl NodeColorWithVertexRule {
    pub fn from_particles<B: AsRef<str>, A: IntoIterator<Item = B>>(
        iter: A,
        model: &Model,
    ) -> Self {
        let mut particles: Vec<ArcParticle> = iter
            .into_iter()
            .map(|p| model.get_particle(p.as_ref()))
            .collect();
        particles.sort_by_key(|p| p.0.pdg_code);
        let vertex_rule = model.particle_set_to_vertex_rules_map[&particles][0].clone();
        Self {
            external_tag: 0,
            vertex_rule,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default, Copy)]
pub struct NodeColorWithoutVertexRule {
    pub external_tag: i32,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct EdgeColor {
    pub pdg: isize,
}

impl std::fmt::Display for EdgeColor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pdg)
    }
}

impl EdgeColor {
    pub(crate) fn from_particle(particle: ArcParticle) -> Self {
        Self {
            pdg: particle.0.pdg_code,
        }
    }
}

fn polyrat_to_atom(
    polyrat: &symbolica::domains::rational::Fraction<
        symbolica::poly::polynomial::PolynomialRing<
            AlgebraicExtension<
                symbolica::domains::rational::FractionField<
                    symbolica::domains::integer::IntegerRing,
                >,
            >,
            u16,
        >,
    >,
) -> Atom {
    let polyrat_numerator = polyrat.numerator();
    let num = polyrat_numerator.to_expression_with_coeff_map(|_x, a, b| {
        let a_poly = a.clone().into_poly();
        *b = Atom::num(match a_poly.coefficients.len() {
            0 => Coefficient::Complex(SymbolicaComplex::new(0.into(), 0.into())),
            1 => {
                if a_poly.exponents[0] == 1 {
                    Coefficient::Complex(SymbolicaComplex::new(
                        0.into(),
                        a_poly.coefficients[0].clone(),
                    ))
                } else {
                    Coefficient::Complex(SymbolicaComplex::new(
                        a_poly.coefficients[0].clone(),
                        0.into(),
                    ))
                }
            }
            2 => Coefficient::Complex(SymbolicaComplex::new(
                a_poly.coefficients[0].clone(),
                a_poly.coefficients[1].clone(),
            )),
            _ => unreachable!(),
        });
    });
    let polyrat_denominator = polyrat.denominator();
    let den = polyrat_denominator.to_expression_with_coeff_map(|_x, a, b| {
        let a_poly = a.clone().into_poly();

        *b = Atom::num(match a_poly.coefficients.len() {
            0 => Coefficient::Complex(SymbolicaComplex::new(0.into(), 0.into())),
            1 => {
                if a_poly.exponents[0] == 1 {
                    Coefficient::Complex(SymbolicaComplex::new(
                        0.into(),
                        a_poly.coefficients[0].clone(),
                    ))
                } else {
                    Coefficient::Complex(SymbolicaComplex::new(
                        a_poly.coefficients[0].clone(),
                        0.into(),
                    ))
                }
            }
            2 => Coefficient::Complex(SymbolicaComplex::new(
                a_poly.coefficients[0].clone(),
                a_poly.coefficients[1].clone(),
            )),
            _ => unreachable!(),
        });
    });
    num / den
}

impl std::fmt::Display for NodeColorWithoutVertexRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.external_tag)
    }
}

impl std::fmt::Display for NodeColorWithVertexRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({}|{})",
            self.external_tag,
            self.vertex_rule.0.name.clone()
        )
    }
}

pub(crate) fn follow_chain(
    current_node: usize,
    vetos: &mut Vec<bool>,
    adj_map: &HashMap<usize, Vec<(usize, usize)>>,
    one_step_only: bool,
) -> Result<(Option<usize>, usize), FeynGenError> {
    if let Some(connected_edges) = adj_map.get(&current_node) {
        let targets = connected_edges
            .iter()
            .filter_map(|(i_e, next_node)| {
                if !vetos[*i_e] {
                    Some((i_e, next_node))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if targets.is_empty() {
            Ok((None, current_node))
        } else if targets.len() == 1 {
            let (next_edge, next_node) = targets.first().unwrap();
            vetos[**next_edge] = true;
            if one_step_only {
                Ok((Some(**next_edge), **next_node))
            } else {
                follow_chain(**next_node, vetos, adj_map, one_step_only)
            }
        } else {
            Ok((None, *targets.first().unwrap().1))
            // return Err(FeynGenError::GenericError(
            //     "GammaLoop does not support four-fermion vertices yet".to_string(),
            // ));
        }
    } else {
        Ok((None, current_node))
    }
}

pub fn evaluate_overall_factor(factor: AtomView) -> Atom {
    let mut res = factor.to_owned();
    for header in [
        "AutG",
        "CouplingsMultiplicity",
        "InternalFermionLoopSign",
        "ExternalFermionOrderingSign",
        "AntiFermionSpinSumSign",
        "NumeratorIndependentSymmetryGrouping",
    ] {
        res = res
            .replace(function!(symbol!(header), Atom::var(symbol!("x_"))).to_pattern())
            .with(Atom::var(symbol!("x_")).to_pattern());
    }
    res = res
        .replace(
            function!(
                symbol!("NumeratorDependentGrouping"),
                Atom::var(symbol!("GraphId_")),
                Atom::var(symbol!("ratio_")),
                Atom::var(symbol!("GraphSymmetryFactor_"))
            )
            .to_pattern(),
        )
        .with(
            (Atom::var(symbol!("ratio_")) * Atom::var(symbol!("GraphSymmetryFactor_")))
                .to_pattern(),
        );
    res.expand()
}

pub fn evaluate_sign_origin(factor: AtomView) -> Atom {
    let mut res = factor.to_owned();
    for header in [
        "AutG",
        "CouplingsMultiplicity",
        "InternalFermionLoopSign",
        "ExternalFermionOrderingSign",
        "AntiFermionSpinSumSign",
        "NumeratorIndependentSymmetryGrouping",
    ] {
        res = res
            .replace(function!(symbol!(header), Atom::var(symbol!("x_"))).to_pattern())
            .with(Atom::var(symbol!("x_")).to_pattern());
    }
    res = res
        .replace(function!(
            symbol!("NumeratorDependentGrouping"),
            Atom::var(symbol!("GraphId_")),
            Atom::var(symbol!("ratio_")),
            Atom::var(symbol!("GraphSymmetryFactor_"))
        ))
        .with(function!(
            symbol!("Group"),
            Atom::var(symbol!("GraphId_")),
            Atom::var(symbol!("ratio_")),
            Atom::var(symbol!("GraphSymmetryFactor_"))
        ));
    res.expand()
}

impl ProcessDefinition {
    pub(crate) fn sample_lib(
        &self,
        sample_iterator: &mut PrimeIteratorU64,
        add_model_params: bool,
        symmetric_polarizations: bool,
        model: &Model,
    ) -> (
        Vec<Replacement>,
        TensorLibrary<ParamTensor<ExplicitKey<Aind>>, Aind>,
    ) {
        // let ff = Zp64::new(PrimeIteratorU64::new(100000000).next().unwrap());

        let mut weyl = TensorLibrary::new();
        weyl.update_ids();

        let gamma_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            gamma_data_weyl(AGS.gamma_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        // println!("permutation{}", gamma_key.rep_permutation);
        weyl.insert_explicit(gamma_key);

        let gamma5_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            gamma5_weyl_data(AGS.gamma5_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        weyl.insert_explicit(gamma5_key);

        let projm_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            proj_m_data_weyl(AGS.projm_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        weyl.insert_explicit(projm_key);

        let projp_key = PermutedStructure::identity(ParamTensor::composite(DataTensor::Sparse(
            proj_p_data_weyl(AGS.projp_strct::<Aind>(4), Atom::num(1), Atom::num(0))
                .map_data(|a| a.re + Atom::i() * a.im),
        )));
        weyl.insert_explicit(projp_key);
        let mut lib = weyl;

        for i in 0..self.loop_count_range.1 {
            let key = ExplicitKey::from_iter(
                [Minkowski {}.new_rep(4)],
                GS.loop_mom,
                Some(vec![Atom::num(i)]),
            );

            //debug!("lib_loop:{}", key.clone().permute_with_metric());
            let key = ParamTensor::from_dense(
                key.structure,
                (0..4)
                    .map(|_| Atom::prime_generate_rat(sample_iterator))
                    .collect(),
            )
            .unwrap();

            lib.insert_explicit(PermutedStructure::identity(key));
        }

        let mut pol_vals = vec![];
        for (i, pdg) in self.initial_pdgs.iter().enumerate() {
            let additional_args = Some(vec![Atom::num(i)]);
            let key = ExplicitKey::from_iter(
                [Minkowski {}.new_rep(4)],
                GS.external_mom,
                additional_args.clone(),
            );

            //debug!("lib_ext:{}", key.clone().permute_with_metric());

            let key = ParamTensor::from_dense(
                key.structure,
                (0..4)
                    .map(|_| Atom::prime_generate_rat(sample_iterator))
                    .collect(),
            )
            .unwrap();

            lib.insert_explicit(PermutedStructure::identity(key));

            let p = model.get_particle_from_pdg(*pdg as isize);

            let structure = p.spin_reps();
            let global_name = EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Sink);

            let len = structure.size().unwrap();

            let key = PermutedStructure::identity(ExplicitKey {
                structure,
                global_name,
                additional_args,
            });

            pol_vals.push(
                (0..len)
                    .map(|_| Atom::prime_generate_rat_complex(sample_iterator))
                    .collect::<Vec<_>>(),
            );

            //debug!("lib_pol:{}", key.clone().permute_with_metric());
            let key =
                ParamTensor::from_dense(key.structure, pol_vals.last().unwrap().clone()).unwrap();

            lib.insert_explicit(PermutedStructure::identity(key));
        }
        match self.generation_type {
            GenerationType::Amplitude => {
                let ext_shift = self.initial_pdgs.len();
                for (i, pdg) in self.final_pdgs_lists[0].iter().enumerate() {
                    let additional_args = Some(vec![Atom::num(i + ext_shift)]);
                    let key = ExplicitKey::from_iter(
                        [Minkowski {}.new_rep(4)],
                        GS.external_mom,
                        additional_args.clone(),
                    );

                    let key = ParamTensor::from_dense(
                        key.structure,
                        (0..4)
                            .map(|_| Atom::prime_generate_rat(sample_iterator))
                            .collect(),
                    )
                    .unwrap();

                    lib.insert_explicit(PermutedStructure::identity(key));

                    let p = model.get_particle_from_pdg(*pdg as isize);

                    let structure = p.spin_reps();
                    let global_name = EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Sink);

                    let len = structure.size().unwrap();

                    let key = ExplicitKey {
                        structure,
                        global_name,
                        additional_args,
                    };

                    let key = ParamTensor::from_dense(
                        key,
                        (0..len)
                            .map(|_| Atom::prime_generate_rat_complex(sample_iterator))
                            .collect(),
                    )
                    .unwrap();

                    lib.insert_explicit(PermutedStructure::identity(key));
                }
            }
            GenerationType::CrossSection => {
                for (i, pdg) in self.initial_pdgs.iter().enumerate() {
                    let additional_args = Some(vec![Atom::num(i)]);

                    let p = model.get_particle_from_pdg(*pdg as isize);

                    let structure = p.spin_reps();
                    let global_name =
                        EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Source);

                    // let len = structure.size().unwrap();

                    let key = PermutedStructure::identity(ExplicitKey {
                        structure,
                        global_name,
                        additional_args,
                    });

                    //debug!("lib_pol:{}", key.clone().permute_with_metric());
                    let key = ParamTensor::from_dense(
                        key.structure,
                        pol_vals[i]
                            .iter()
                            .map(|a| {
                                if symmetric_polarizations {
                                    a.clone()
                                } else {
                                    Atom::prime_generate_rat_complex(sample_iterator)
                                }
                            })
                            .collect(),
                    )
                    .unwrap();

                    lib.insert_explicit(PermutedStructure::identity(key));
                }
            }
        }
        let mut reps = vec![];

        if add_model_params {
            for param in model.parameters.values().filter(|p| p.value.is_some()) {
                if param.value.is_some() {
                    let name: Atom = param.name.into();
                    reps.push(Replacement::new(
                        name.to_pattern(),
                        Atom::prime_generate_rat_complex(sample_iterator),
                    ));
                }
            }
        }
        // for r in &reps {
        //     info!("Model replacements: {}", r);
        // }
        (reps, lib)
    }

    #[instrument(skip_all)]

    pub(crate) fn unresolved_cut_content(&self, model: &Model) -> (usize, AHashSet<ArcParticle>) {
        if let Some(p) = self.cross_section_filters.get_perturbative_orders() {
            let mut unresolved = AHashSet::new();
            for k in p.keys() {
                let k: SmartString<_> = k.clone().into();
                if let Some(p) = model.unresolved_particles.get(&k) {
                    unresolved = unresolved.union(p).cloned().collect();
                }
            }
            (p.values().sum(), unresolved)
        } else {
            (0, AHashSet::new())
        }
    }

    pub(crate) fn canonize_external_momenta_assignment(
        &self,
        model: &Model,
        node_colors_for_external_symmetrization: &HashMap<i32, i32>,
        graph: &mut SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    ) {
        let mut initial_pdgs = self
            .initial_pdgs
            .iter()
            .enumerate()
            .map(|(i_n, pdg)| (pdg, i_n + 1))
            .collect::<Vec<_>>();
        let mut final_pdgs = if matches!(self.generation_type, GenerationType::CrossSection) {
            self.initial_pdgs
                .iter()
                .enumerate()
                .map(|(i_n, pdg)| (pdg, i_n + 1 + initial_pdgs.len()))
                .collect::<Vec<_>>()
        } else {
            self.final_pdgs_lists[0]
                .iter()
                .enumerate()
                .map(|(i_n, pdg)| (pdg, i_n + 1 + initial_pdgs.len()))
                .collect::<Vec<_>>()
        };

        let mut all_pdgs = initial_pdgs.clone();
        all_pdgs.extend(final_pdgs.clone());

        let mut new_node_data = vec![];
        let mut new_edge_data = vec![];

        // We do two passes, first assigning "specified externals" only (i.e. with tags > 0) and finally the remaining non specified ones (tags < 0)
        // The three pass steps are:
        // 0: distribute all externals with forced assignment, i.e. external_tag > 0
        // 1: distribute all externals with forced assignment up to left-right symmetry, i.e. external_tag < -2000
        // 2: distribute symmetrized externals, i.e. external_tag > -2000 && external_tag < 0

        for pass_steps in [0, 1, 2] {
            for (i_e, e) in graph.edges().iter().enumerate() {
                // All edges supposed to be incoming at this stage
                assert!(graph.nodes()[e.vertices.1].data.external_tag == 0);
                assert!(graph.nodes()[e.vertices.0].data.external_tag >= 0);
                let p = model.get_particle_from_pdg(e.data.pdg);
                let external_tag = graph.nodes()[e.vertices.0].data.external_tag;
                let symmetrized_external_tag = node_colors_for_external_symmetrization
                    .get(&external_tag)
                    .copied()
                    .unwrap_or(external_tag);
                let is_initial_state = external_tag <= self.initial_pdgs.len() as i32;
                let container = if is_initial_state {
                    &mut initial_pdgs
                } else {
                    &mut final_pdgs
                };
                if pass_steps == 0 && symmetrized_external_tag > 0 {
                    if self.symmetrize_left_right_states {
                        let matched_external_pos = all_pdgs
                            .iter()
                            .position(|(_pdg, i_ext)| (*i_ext as i32) == symmetrized_external_tag)
                            .unwrap();
                        all_pdgs.remove(matched_external_pos);
                    } else {
                        let matched_external_pos = container
                            .iter()
                            .position(|(_pdg, i_ext)| (*i_ext as i32) == symmetrized_external_tag)
                            .unwrap();
                        container.remove(matched_external_pos);
                    };
                } else if pass_steps == 1 && symmetrized_external_tag < -2000 {
                    // Node colors below -2000 indicate a left-right symmetrization of the external legs for a forward-scattering diagrams
                    // without symmetrizing the initial-states.
                    let external_leg_position = (-symmetrized_external_tag) % 1000;
                    // Try and find this either in initial or final states
                    let (matched_position, is_initial_match) = if let Some(matched_initial_pos) =
                        all_pdgs
                            .iter()
                            .position(|(_pdg, i_ext)| (*i_ext as i32) == external_leg_position)
                    {
                        (matched_initial_pos, true)
                    } else if let Some(matched_final_pos) =
                        all_pdgs.iter().position(|(_pdg, i_ext)| {
                            (*i_ext as i32)
                                == external_leg_position + (self.initial_pdgs.len() as i32)
                        })
                    {
                        (matched_final_pos, false)
                    } else {
                        unreachable!(
                            "Logical mistake in feyngen: external legs in canonicalized graphs should always be matchable."
                        )
                    };

                    let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                    let new_external_tag = all_pdgs[matched_position].1 as i32;
                    // If we swapped initial and final state assignment, then we must also flip the pdg code of the corresponding half-edges
                    if is_initial_state != is_initial_match {
                        let mut e_data = e.data;
                        e_data.pdg = model
                            .get_particle_from_pdg(e_data.pdg)
                            .0
                            .get_anti_particle(model)
                            .0
                            .pdg_code;
                        new_edge_data.push((i_e, e_data));
                    }
                    new_data.external_tag = new_external_tag;
                    new_node_data.push((e.vertices.0, new_data));
                    all_pdgs.remove(matched_position);
                } else if pass_steps == 2 && (-2000..0).contains(&symmetrized_external_tag) {
                    let pdg_code = if is_initial_state {
                        p.0.pdg_code
                    } else {
                        p.0.get_anti_particle(model).0.pdg_code
                    };
                    if self.symmetrize_left_right_states {
                        let matched_external_pos: usize = all_pdgs
                            .iter()
                            .position(|(pdg, _i_ext)| **pdg == pdg_code as i64)
                            .unwrap();
                        let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                        let new_external_tag = all_pdgs[matched_external_pos].1 as i32;
                        // If we swapped initial and final state assignment, then we must also flip the pdg code of the corresponding half-edges
                        if (new_external_tag > initial_pdgs.len() as i32
                            && new_data.external_tag <= initial_pdgs.len() as i32)
                            || (new_external_tag <= initial_pdgs.len() as i32
                                && new_data.external_tag > initial_pdgs.len() as i32)
                        {
                            let mut e_data = e.data;
                            e_data.pdg = model
                                .get_particle_from_pdg(e_data.pdg)
                                .0
                                .get_anti_particle(model)
                                .0
                                .pdg_code;
                            new_edge_data.push((i_e, e_data));
                        }
                        new_data.external_tag = new_external_tag;
                        new_node_data.push((e.vertices.0, new_data));
                        all_pdgs.remove(matched_external_pos);
                    } else {
                        let matched_external_pos = container
                            .iter()
                            .position(|(pdg, _i_ext)| **pdg == pdg_code as i64)
                            .unwrap();
                        let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                        new_data.external_tag = container[matched_external_pos].1 as i32;
                        new_node_data.push((e.vertices.0, new_data));
                        container.remove(matched_external_pos);
                    }
                }
            }
        }

        for (node_pos, new_data) in new_node_data {
            graph.set_node_data(node_pos, new_data);
        }
        for (edge_pos, new_data) in new_edge_data {
            graph.set_edge_data(edge_pos, new_data);
        }
    }

    /// This function canonizes the edge and vertex ordering of a graph based on the skeletton graph with only propagator mass as edge color.
    /// This is useful to then allow for further grouping of isomorphic graphs, incl numerator.
    #[allow(clippy::type_complexity)]
    pub(crate) fn canonicalize_edge_and_vertex_ordering(
        &self,
        model: &Model,
        input_graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        node_colors_for_external_symmetrization: &HashMap<i32, i32>,
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        // This option contains: (self.symmetrize_initial_states, self.symmetrize_left_right_states) if any is true, else None.
        manually_canonalize_initial_states_cross_section_ordering: Option<(bool, bool)>,
    ) -> Result<
        (
            SymbolicaGraph<NodeColorWithoutVertexRule, std::string::String>,
            SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        ),
        FeynGenError,
    > {
        // println!("INPUT GRAPH:\n{}", input_graph.to_dot());
        let mut canonized_skelettons = vec![];
        let external_node_positions = input_graph
            .nodes()
            .iter()
            .enumerate()
            .filter_map(|(i_n, n)| {
                if n.data.external_tag > 0 {
                    Some((n.data.external_tag, i_n))
                } else {
                    None
                }
            })
            .collect::<HashMap<_, _>>();
        // Contains the two tuple (are_left_and_right_swapped, remapping)
        let mut external_remappings = vec![];
        if let Some((symmetrize_initial_states, symmetrize_left_right_states)) =
            manually_canonalize_initial_states_cross_section_ordering
        {
            let external_ids = (1..=self.initial_pdgs.len()).collect::<Vec<_>>();
            let mut initial_states_orders = vec![];
            if symmetrize_initial_states {
                'permutation_loop: for permutation in external_ids
                    .iter()
                    .cloned()
                    .permutations(external_ids.len())
                {
                    if permutation
                        .iter()
                        .enumerate()
                        .any(|(src, trgt)| self.initial_pdgs[src] != self.initial_pdgs[*trgt - 1])
                    {
                        continue 'permutation_loop;
                    }
                    initial_states_orders.push(permutation);
                }
            } else {
                initial_states_orders.push(external_ids);
            };
            for initial_state_order in initial_states_orders {
                let mut remapping = AHashMap::<usize, usize>::default();
                for ext_id in 1..=self.initial_pdgs.len() {
                    remapping.insert(
                        *external_node_positions.get(&(ext_id as i32)).unwrap(),
                        *external_node_positions
                            .get(&((initial_state_order[ext_id - 1]) as i32))
                            .unwrap(),
                    );
                    remapping.insert(
                        *external_node_positions
                            .get(&((ext_id + self.initial_pdgs.len()) as i32))
                            .unwrap(),
                        *external_node_positions
                            .get(
                                &((initial_state_order[ext_id - 1] + self.initial_pdgs.len())
                                    as i32),
                            )
                            .unwrap(),
                    );
                }
                external_remappings.push((false, remapping));
                if symmetrize_left_right_states {
                    let mut remapping = AHashMap::<usize, usize>::default();
                    for ext_id in 1..=self.initial_pdgs.len() {
                        remapping.insert(
                            *external_node_positions.get(&(ext_id as i32)).unwrap(),
                            *external_node_positions
                                .get(
                                    &((initial_state_order[ext_id - 1] + self.initial_pdgs.len())
                                        as i32),
                                )
                                .unwrap(),
                        );
                        remapping.insert(
                            *external_node_positions
                                .get(&((ext_id + self.initial_pdgs.len()) as i32))
                                .unwrap(),
                            *external_node_positions
                                .get(&(initial_state_order[ext_id - 1] as i32))
                                .unwrap(),
                        );
                    }
                    external_remappings.push((true, remapping));
                }
            }
        } else {
            external_remappings.push((false, AHashMap::<usize, usize>::default()));
        }

        for (are_left_and_right_swapped, remapping) in external_remappings {
            // Make sure to canonize the edge ordering based on the skeletton graph with only propagator mass as edge color
            let mut skeletton_graph = SymbolicaGraph::new();
            for node in input_graph.nodes() {
                skeletton_graph.add_node(NodeColorWithoutVertexRule {
                    external_tag: *node_colors_for_external_symmetrization
                        .get(&node.data.external_tag)
                        .unwrap_or(&node.data.external_tag),
                });
            }
            let color_according_to_mass: bool = if let Some(grouping_options) =
                numerator_aware_isomorphism_grouping.get_options()
            {
                grouping_options.differentiate_particle_masses_only
            } else {
                true
            };
            for edge in input_graph.edges() {
                // We must maintain the directionality of the external edges to make sure that antifermions are not mapped to fermions
                // when grouping is done based on the skeletton graph and only the mass of the propagators
                // Also external edges must always retain their PDG since exchanging the momentum of e.g. a massless electron and muon
                // would not be a valid isomorphism w.r.t the process definition
                let is_edge_external = [edge.vertices.0, edge.vertices.1]
                    .iter()
                    .any(|e| input_graph.nodes()[*e].data.external_tag != 0);
                let mut remapped_edge_vertices: (usize, usize) = (
                    *remapping.get(&edge.vertices.0).unwrap_or(&edge.vertices.0),
                    *remapping.get(&edge.vertices.1).unwrap_or(&edge.vertices.1),
                );
                // Force all externals incoming in the canonicalization
                let particle = if input_graph.nodes()[edge.vertices.1].data.external_tag > 0 {
                    remapped_edge_vertices = (remapped_edge_vertices.1, remapped_edge_vertices.0);
                    model
                        .get_particle_from_pdg(edge.data.pdg)
                        .0
                        .get_anti_particle(model)
                } else {
                    model.get_particle_from_pdg(edge.data.pdg)
                };

                let keep_direction = is_edge_external;

                // // We want to forget the orientation of fermions (except for externals) to capture furry, but we need the canonization to be sensitive to charged vector bosons when it is connected to fermion lines
                // if !keep_direction && edge.directed && particle.is_vector() {
                //     let src_ferm = input_graph.node(edge.vertices.0).edges.iter().any(|a| {
                //         model
                //             .get_particle_from_pdg(input_graph.edge(*a).data.pdg)
                //             .is_fermion()
                //     });
                //     let sink_ferm = input_graph.node(edge.vertices.1).edges.iter().any(|a| {
                //         model
                //             .get_particle_from_pdg(input_graph.edge(*a).data.pdg)
                //             .is_fermion()
                //     });

                //     if src_ferm || sink_ferm {
                //         keep_direction = true;
                //     }
                // }

                // WARNING: It is important to note that the colouring choice below dictates what representative diagram (i.e. "sorted_g") will be used
                // for performing numerical comparisons for grouping isomorphic graphs. If we do not include the spin in the colouring, we may incorrectly
                // sort two isomorphic graphs with and interchange of massless quarks and gluons which will prevent their grouping (final result still correct, but more diagrams).
                // This is avoided by including the spin in the color, which will prevent massless quark and gluons from ever being interchanged.
                // Of course, even then, it could be that we pick two sorted representatives that are not isomorphic, even though there exist a different isomorphic skeletton graph
                // which would have given a match when doing the numerical comparison. In the SM, this should never happen, and for BSM we can afford to lose some grouping.
                // The only way to avoid this would be to numerically test all isomorphic permutations of the skeletton graph, which is prohibitively slow.
                skeletton_graph
                    .add_edge(
                        remapped_edge_vertices.0,
                        remapped_edge_vertices.1,
                        keep_direction,
                        if color_according_to_mass && !is_edge_external {
                            format!("{} | {}", particle.0.mass.0, particle.0.spin)
                        } else {
                            particle.0.name.to_string()
                        },
                    )
                    .unwrap();
            }
            canonized_skelettons.push((
                (are_left_and_right_swapped, remapping),
                skeletton_graph.canonize(),
            ));
        }
        canonized_skelettons.sort_by(|a, b| {
            (a.1.graph.nodes(), a.1.graph.edges()).cmp(&(b.1.graph.nodes(), b.1.graph.edges()))
        });
        let ((was_left_and_right_swapped, selected_external_remapping), canonized_skeletton) =
            canonized_skelettons.first().unwrap();

        let mut can_graph_node_pos_to_input_graph_node_pos: Vec<usize> =
            vec![0; input_graph.nodes().len()];
        for (input_graph_node_position, node_order) in
            canonized_skeletton.vertex_map.iter().enumerate()
        {
            can_graph_node_pos_to_input_graph_node_pos[*node_order] = *selected_external_remapping
                .get(&input_graph_node_position)
                .unwrap_or(&input_graph_node_position);
        }
        let mut input_graph_node_pos_to_can_graph_node_pos: Vec<usize> =
            vec![0; input_graph.nodes().len()];
        for (can_graph_node_pos, input_graph_node_pos) in can_graph_node_pos_to_input_graph_node_pos
            .iter()
            .enumerate()
        {
            input_graph_node_pos_to_can_graph_node_pos[*input_graph_node_pos] = can_graph_node_pos;
        }

        // Sort nodes according to the canonized skeleton graph
        // This will also ensure that EMR variables line up
        let mut sorted_g: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor> =
            SymbolicaGraph::new();
        for &node_order in can_graph_node_pos_to_input_graph_node_pos.iter() {
            let node_data = input_graph.nodes()[*selected_external_remapping
                .get(&node_order)
                .unwrap_or(&node_order)]
            .data
            .clone();
            sorted_g.add_node(node_data);
        }

        let input_graph_nodes = input_graph.nodes();
        let mut reordered_edges = input_graph
            .edges()
            .iter()
            .map(|e| {
                // Set all external edges incoming and internal edges going from lower to higher node order
                let is_external = input_graph_nodes[e.vertices.0].data.external_tag > 0
                    || input_graph_nodes[e.vertices.1].data.external_tag > 0;
                let is_external_outgoing = input_graph_nodes[e.vertices.1].data.external_tag > 0;
                let is_flipped = is_external_outgoing
                    || (!is_external
                        && input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                            > input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]);
                let mut particle = if is_flipped {
                    model
                        .get_particle_from_pdg(e.data.pdg)
                        .0
                        .get_anti_particle(model)
                } else {
                    model.get_particle_from_pdg(e.data.pdg)
                };
                // Apply the switch from particle to anti-particle (CP symmetry) if the canonicalization swapped initial and final states.
                // IMPORTANT: we must here to apply CP not just to the external particles since the assignment of fermion flow matters, for
                // internal edges connecting these external fermions.
                // If one would fix those edges (necessary e.g. for e- d > e- d DIS process) then one could add '&& is_external' below, which
                // would allow to capture additional groupings involving the charged. W bosons.
                if *was_left_and_right_swapped {
                    particle = particle.0.get_anti_particle(model);
                }
                (
                    (e, is_flipped),
                    // This key will serve to give a unique ordering of the edges
                    (
                        if !is_flipped {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                        } else {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]
                        },
                        if !is_flipped {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]
                        } else {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                        },
                        particle.0.pdg_code,
                    ),
                )
            })
            .collect::<Vec<_>>();
        reordered_edges.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        for &((_e, _is_flipped), (v_in, v_out, pdg)) in reordered_edges.iter() {
            // We must canonize the fermion flow as well, so we force the fermion flow to always go from the lower to the higher node order,
            // and adjust the edge colour (particle vs anti-particle) accordingly.
            // This is necessary to ensure that the meaning of the edge momentum representation (which always aligns the monentum with the fermion flow)
            // is preserved
            sorted_g
                .add_edge(
                    v_in,
                    v_out,
                    // Now that the edge orientation is canonized we set all edges as directed because
                    // it is semantic even for bosons as it dictates the flow of the momentum
                    true,
                    EdgeColor::from_particle(model.get_particle_from_pdg(pdg)),
                )
                .unwrap();
        }

        // In order for the external assignment momenta to line up between the canonized versions of these graphs
        self.canonize_external_momenta_assignment(
            model,
            node_colors_for_external_symmetrization,
            &mut sorted_g,
        );

        Ok((canonized_skeletton.graph.to_owned(), sorted_g))
    }

    // This function normalize fermion and ghost flows and makes sure that all virtual charged particles without flows (like the W boson) are particles
    pub(crate) fn normalize_flows(
        &self,
        graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        model: &Model,
    ) -> Result<(SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>, bool), FeynGenError> {
        let mut adj_map: HashMap<usize, Vec<(usize, usize)>> = HashMap::default();
        for (i_e, e) in graph.edges().iter().enumerate() {
            // Build an adjacency list including only fermions
            let p = model.get_particle_from_pdg(e.data.pdg);
            if !(p.0.is_fermion() || p.0.is_ghost()) {
                continue;
            }
            adj_map
                .entry(e.vertices.0)
                .or_default()
                .push((i_e, e.vertices.1));
            if e.vertices.0 != e.vertices.1 {
                adj_map
                    .entry(e.vertices.1)
                    .or_default()
                    .push((i_e, e.vertices.0));
            }
        }

        let mut vetoed_edges: Vec<bool> = graph
            .edges()
            .iter()
            .map(|e| {
                let p = model.get_particle_from_pdg(e.data.pdg);
                !(p.0.is_fermion() || p.0.is_ghost())
            })
            .collect();
        let mut new_edges: AHashMap<usize, (usize, usize, bool, EdgeColor)> = AHashMap::default();
        let mut normalized_graph: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor> =
            SymbolicaGraph::new();
        for n in graph.nodes().iter() {
            normalized_graph.add_node(n.data.clone());
        }

        let graph_edges = graph.edges();
        // Pairing of the external fermion flows. The keys are (sorted) two-tuple of the PDGs of the external fermions (assuming all incoming)
        // and the values are the external leg ids of the fermions connected and in the order of the sorted key.
        #[allow(clippy::type_complexity)]
        let mut external_fermion_flow_pairings: AHashMap<
            (ArcParticle, ArcParticle),
            Vec<(usize, usize)>,
        > = AHashMap::default();
        // First fix flows connected to external only and after that fix all internal fermion/ghost flows
        let mut external_tags_ordering = (1..=self.initial_pdgs.len()).collect::<Vec<_>>();
        if self.generation_type == GenerationType::CrossSection {
            external_tags_ordering
                .extend((1..=self.initial_pdgs.len()).map(|i| i + self.initial_pdgs.len()));
        } else {
            external_tags_ordering
                .extend((1..=self.final_pdgs_lists[0].len()).map(|i| i + self.initial_pdgs.len()));
        }
        external_tags_ordering.push(0);
        for external_tag_to_consider in external_tags_ordering {
            for (i_e, e) in graph_edges.iter().enumerate() {
                let edge_tag = if graph.nodes()[e.vertices.0].data.external_tag != 0 {
                    graph.nodes()[e.vertices.0].data.external_tag
                } else {
                    graph.nodes()[e.vertices.1].data.external_tag
                };
                assert!(edge_tag >= 0);

                let is_a_virtual_edge = edge_tag == 0;
                if external_tag_to_consider != (edge_tag as usize) {
                    continue;
                }
                if vetoed_edges[i_e] {
                    if !new_edges.contains_key(&i_e) {
                        new_edges.insert(i_e, (e.vertices.0, e.vertices.1, e.directed, e.data));
                    }
                    continue;
                }
                let starting_particle = model.get_particle_from_pdg(e.data.pdg);
                let mut is_starting_antiparticle = starting_particle.0.is_antiparticle();
                let starting_vertices = if is_a_virtual_edge && is_starting_antiparticle {
                    // Force all virtual closed fermion loops to be particles
                    is_starting_antiparticle = false;
                    if !new_edges.contains_key(&i_e) {
                        new_edges.insert(
                            i_e,
                            (
                                e.vertices.1,
                                e.vertices.0,
                                e.directed,
                                EdgeColor::from_particle(
                                    model
                                        .get_particle_from_pdg(e.data.pdg)
                                        .0
                                        .get_anti_particle(model),
                                ),
                            ),
                        );
                        (e.vertices.1, e.vertices.0)
                    } else {
                        continue;
                    }
                } else if !new_edges.contains_key(&i_e) {
                    new_edges.insert(i_e, (e.vertices.0, e.vertices.1, e.directed, e.data));
                    (e.vertices.0, e.vertices.1)
                } else {
                    continue;
                };

                vetoed_edges[i_e] = true;
                let mut connected_leg_ids: AHashSet<usize> = AHashSet::new();
                for read_to_the_right in [true, false] {
                    let mut previous_node_position = if read_to_the_right {
                        starting_vertices.1
                    } else {
                        starting_vertices.0
                    };
                    // println!(
                    //     "Starting reading chain from edge {}->{}, from {}",
                    //     starting_vertices.0, starting_vertices.1, previous_node_position
                    // );
                    if graph.nodes()[previous_node_position].data.external_tag > 0 {
                        connected_leg_ids.insert(
                            graph.nodes()[previous_node_position].data.external_tag as usize,
                        );
                    }

                    'fermion_chain: loop {
                        let (next_fermion_edge_position, next_fermion_node_position) =
                            follow_chain(
                                previous_node_position,
                                &mut vetoed_edges,
                                &adj_map,
                                true,
                            )?;
                        if graph.nodes()[next_fermion_node_position].data.external_tag > 0 {
                            connected_leg_ids.insert(
                                graph.nodes()[next_fermion_node_position].data.external_tag
                                    as usize,
                            );
                        }
                        // println!(
                        //     "Next edge: {}",
                        //     if let Some(nfep) = next_fermion_edge_position {
                        //         format!(
                        //             "{}, {} -> {}",
                        //             graph_edges[nfep].vertices.0, graph_edges[nfep].vertices.1, nfep
                        //         )
                        //     } else {
                        //         "None".to_string()
                        //     }
                        // );
                        if let Some(nfep) = next_fermion_edge_position {
                            let this_has_same_orientation_starting = (read_to_the_right
                                && graph_edges[nfep].vertices.1 == next_fermion_node_position)
                                || (!read_to_the_right
                                    && graph_edges[nfep].vertices.0 == next_fermion_node_position);
                            if !new_edges.contains_key(&nfep) {
                                if this_has_same_orientation_starting {
                                    new_edges.insert(
                                        nfep,
                                        (
                                            graph_edges[nfep].vertices.0,
                                            graph_edges[nfep].vertices.1,
                                            graph_edges[nfep].directed,
                                            graph_edges[nfep].data,
                                        ),
                                    );
                                } else {
                                    let this_edge_particle =
                                        model.get_particle_from_pdg(graph_edges[nfep].data.pdg);
                                    new_edges.insert(
                                        nfep,
                                        (
                                            graph_edges[nfep].vertices.1,
                                            graph_edges[nfep].vertices.0,
                                            graph_edges[nfep].directed,
                                            // Force fermion to be the same species as the one starting the chain
                                            if this_edge_particle.0.is_antiparticle()
                                                != is_starting_antiparticle
                                            {
                                                EdgeColor::from_particle(
                                                    this_edge_particle.0.get_anti_particle(model),
                                                )
                                            } else {
                                                EdgeColor::from_particle(this_edge_particle)
                                            },
                                        ),
                                    );
                                }
                            }
                            previous_node_position = next_fermion_node_position;
                        } else {
                            break 'fermion_chain;
                        }
                    }
                }
                if external_tag_to_consider > 0 && starting_particle.0.is_fermion() {
                    if connected_leg_ids.len() != 2 {
                        return Err(FeynGenError::GenericError(
                            "External fermion flow must have exactly two legs".to_string(),
                        ));
                    }
                    let connected_leg_ids_vec =
                        connected_leg_ids.iter().copied().collect::<Vec<_>>();
                    let connected_leg_pdgs = (0..=1)
                        .map(|i| {
                            if connected_leg_ids_vec[i] <= self.initial_pdgs.len() {
                                model.get_particle_from_pdg(
                                    self.initial_pdgs[connected_leg_ids_vec[i] - 1] as isize,
                                )
                            } else {
                                // Assign line types assuming all incoming particles
                                let right_side_pdgs =
                                    if self.generation_type == GenerationType::CrossSection {
                                        &self.initial_pdgs
                                    } else {
                                        &self.final_pdgs_lists[0]
                                    };
                                model
                                    .get_particle_from_pdg(
                                        right_side_pdgs
                                            [connected_leg_ids_vec[i] - self.initial_pdgs.len() - 1]
                                            as isize,
                                    )
                                    .0
                                    .get_anti_particle(model)
                            }
                        })
                        .collect::<Vec<_>>();
                    if connected_leg_pdgs
                        .iter()
                        .any(|particle| !particle.0.is_fermion())
                    {
                        return Err(FeynGenError::GenericError(
                            "External fermion flow must connect two fermions".to_string(),
                        ));
                    }
                    if connected_leg_pdgs[0].0.is_antiparticle()
                        && !connected_leg_pdgs[1].0.is_antiparticle()
                    {
                        external_fermion_flow_pairings
                            .entry((connected_leg_pdgs[0].clone(), connected_leg_pdgs[1].clone()))
                            .or_default()
                            .push((connected_leg_ids_vec[0], connected_leg_ids_vec[1]));
                    } else if connected_leg_pdgs[1].0.is_antiparticle()
                        && !connected_leg_pdgs[0].0.is_antiparticle()
                    {
                        external_fermion_flow_pairings
                            .entry((connected_leg_pdgs[1].clone(), connected_leg_pdgs[0].clone()))
                            .or_default()
                            .push((connected_leg_ids_vec[1], connected_leg_ids_vec[0]));
                    } else {
                        return Err(FeynGenError::GenericError(
                            "External fermion flow must connect a fermion and an anti-fermion. GammaLoop has no support for Majorana particles yet.".to_string(),
                        ));
                    }
                }
            }
        }

        let mut external_fermion_flow_sign = 1;
        let mut concatenated_lines = vec![];
        for (_line_type, lines) in external_fermion_flow_pairings.iter().sorted() {
            for (a, b) in lines.iter() {
                concatenated_lines.push(*a);
                concatenated_lines.push(*b);
            }
        }
        if Permutation::sort(&concatenated_lines)
            .transpositions()
            .len()
            % 2
            == 1
        {
            external_fermion_flow_sign *= -1;
        }

        // Finally make sure that all virtual non-self-antiparticles that are not fermions or ghost (like the W-boson) are set to particles
        for (i_e, e) in graph_edges.iter().enumerate() {
            let is_a_virtual_edge = graph.nodes()[e.vertices.0].data.external_tag == 0
                && graph.nodes()[e.vertices.1].data.external_tag == 0;
            let this_edge_particle = model.get_particle_from_pdg(e.data.pdg);
            if is_a_virtual_edge
                && this_edge_particle.0.is_antiparticle()
                && !(this_edge_particle.0.is_fermion() || this_edge_particle.0.is_ghost())
            {
                new_edges.insert(
                    i_e,
                    (
                        e.vertices.1,
                        e.vertices.0,
                        e.directed,
                        EdgeColor::from_particle(this_edge_particle.0.get_anti_particle(model)),
                    ),
                );
            }
        }

        let mut sorted_new_edges = new_edges.iter().collect::<Vec<_>>();
        sorted_new_edges.sort_by_key(|(i_e_0, _)| *i_e_0);

        for (v0, v1, directed, edge_data) in sorted_new_edges.iter().map(|(_, v)| v) {
            normalized_graph
                .add_edge(*v0, *v1, *directed, *edge_data)
                .map_err(|e| FeynGenError::GenericError(e.to_string()))?;
        }
        assert!(normalized_graph.edges().len() == graph.edges().len());

        Ok((normalized_graph, external_fermion_flow_sign == -1))
    }

    #[instrument(skip_all)]
    pub(crate) fn compare_numerator_tensors(
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        numerator_a: &ProcessedNumeratorForComparison,
        numerator_b: &ProcessedNumeratorForComparison,
    ) -> Option<Atom> {
        if numerator_a.canonized_numerator.is_none() && numerator_b.canonized_numerator.is_some()
            || numerator_a.canonized_numerator.is_some()
                && numerator_b.canonized_numerator.is_none()
        {
            panic!(
                "Inconsistent state: one sample has canonalized numerator while the other does not."
            );
        }

        if numerator_a.sample_evaluations.len() != numerator_b.sample_evaluations.len() {
            panic!(
                "Inconsistent state: the two samples have different number of numerical evaluations."
            );
        }

        match numerator_aware_isomorphism_grouping {
            NumeratorAwareGraphGroupingOption::NoGrouping => None,
            NumeratorAwareGraphGroupingOption::OnlyDetectZeroes => None,
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToSign(_) => {
                numerator_a.compare_with_sign_only(numerator_b)
            }
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(_) => {
                numerator_a.compare_with_scalar_rescaling(numerator_b)
            }
        }
    }

    pub(crate) fn substitute_color_factors(expr: AtomView) -> Atom {
        let replacements = vec![
            (Atom::var(CS.nc), Atom::num(3)),
            (Atom::var(CS.tr), parse!("1/2")),
            // (parse!("CA"), parse!("3")),
            // (parse!("CF"), parse!("4/3")),
        ];
        let mut res = expr.to_owned();
        for (src, trgt) in replacements {
            res = res.replace(src.to_pattern()).with(trgt.to_pattern());
        }
        res
    }
}

#[derive(Clone)]
pub(crate) struct ProcessedNumeratorForComparison {
    diagram_id: usize,
    /// Canonized numerator expression for symbolic comparison.
    ///
    /// This is `None` when:
    /// - `--compare-canonized-numerator` flag is disabled (test_canonized_numerator = false)
    /// - Numerator grouping is disabled (NoGrouping or OnlyDetectZeroes variants)
    /// - No grouping options are available (numerator_aware_isomorphism_grouping.get_options() returns None)
    ///
    /// This is `Some(canonized_atom)` when:
    /// - `--compare-canonized-numerator` flag is enabled (test_canonized_numerator = true)
    /// - The numerator has been processed through loop momentum replacements
    /// - Optionally processed through color factor substitution (if --fully-numerical-substitution-when-comparing-numerators)
    /// - Optionally processed through factor collection (for GroupIdenticalGraphUpToScalarRescaling mode)
    canonized_numerator: Option<Atom>,

    /// Numerical evaluations of the numerator at sample points for numerical comparison.
    ///
    /// This is empty (`vec![]`) when:
    /// - `--number-of-samples-for-numerator-comparisons 0` is specified
    /// - Numerator grouping is disabled (NoGrouping or OnlyDetectZeroes variants)
    /// - No grouping options are available (numerator_aware_isomorphism_grouping.get_options() returns None)
    /// - Neither test_canonized_numerator nor number_of_numerical_samples > 0
    ///
    /// This contains evaluated `Atom`s when:
    /// - `--number-of-samples-for-numerator-comparisons N` where N > 0
    /// - Each element represents the numerator evaluated at a different sample point
    /// - Sample points are generated using a prime-based iterator seeded by numerical_sample_seed
    /// - Evaluations include color expansion, tensor network execution, and color factor substitution
    /// - May have common factors collected depending on the grouping mode and ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL setting
    sample_evaluations: Vec<Atom>,
    sample_evaluations_as_polynomial: Vec<
        symbolica::poly::polynomial::MultivariatePolynomial<
            AlgebraicExtension<
                symbolica::domains::rational::FractionField<
                    symbolica::domains::integer::IntegerRing,
                >,
            >,
        >,
    >,
    sample_evaluations_are_zero: Vec<bool>,
}

impl ProcessedNumeratorForComparison {
    pub(crate) fn evaluated_samples_are_zero(&self) -> bool {
        !self.sample_evaluations.is_empty()
            && self
                .sample_evaluations_are_zero
                .iter()
                .all(|&is_zero| is_zero)
    }

    /// Compare two numerators allowing only sign changes (±1 ratios).
    #[instrument(skip_all, fields(self_id = %self.diagram_id, other_id = %other.diagram_id))]
    fn compare_with_sign_only(&self, other: &ProcessedNumeratorForComparison) -> Option<Atom> {
        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            "Starting sign-only comparison between diagrams"
        );
        fn analyze_diff_and_sum(a: AtomView, b: AtomView) -> Option<Atom> {
            if (a - b).expand().is_zero() {
                return Some(Atom::num(1));
            }
            if (a + b).expand().is_zero() {
                return Some(Atom::num(-1));
            }
            debug!(a = %a.floatify(13).to_canonical_string(),b=%b.floatify(13).to_canonical_string(),"compared but no luck");
            None
        }

        // Try canonized numerator comparison first
        if let (Some(canonized_num_a), Some(canonized_num_b)) = (
            self.canonized_numerator.as_ref(),
            other.canonized_numerator.as_ref(),
        ) {
            debug!("Attempting symbolic comparison using canonized numerators");
            debug!(
                numerator = %canonized_num_a.to_canonical_string(),
                numerator_diagram_id = %self.diagram_id,
                "Canonized numerator for comparison"
            );
            debug!(
                denominator = %canonized_num_b.to_canonical_string(),
                denominator_diagram_id = %other.diagram_id,
                "Canonized denominator for comparison"
            );
            if let Some(ratio) =
                analyze_diff_and_sum(canonized_num_a.as_view(), canonized_num_b.as_view())
            {
                debug!(
                    ratio = %ratio.to_canonical_string(),
                    method = "canonized_numerators",
                    "Successfully matched diagrams using canonized numerators"
                );
                return Some(ratio);
            } else {
                debug!(
                    comparison_type = "canonized_numerator",
                    result = "failed",
                    reason = "expressions_not_identical_up_to_sign",
                    "Canonized numerator comparison failed"
                );
            }
        } else {
            debug!(
                comparison_type = "canonized_numerator",
                result = "skipped",
                reason = "missing_numerators",
                "Skipping canonized numerator comparison"
            );
        }

        // Fall back to sample evaluations
        if !self.sample_evaluations.is_empty() {
            debug!(
                comparison_type = "numerical_samples",
                sample_count = %self.sample_evaluations.len(),
                "Attempting numerical comparison using sample evaluations"
            );
            let ratios = self
                .sample_evaluations
                .iter()
                .zip(other.sample_evaluations.iter())
                .enumerate()
                .map(|(idx, (a, b))| {
                    debug!(
                        sample_idx = %idx,
                        numerator = %a.to_canonical_string(),
                        numerator_diagram_id = %self.diagram_id,
                        "Sample numerator evaluation"
                    );
                    debug!(
                        sample_idx = %idx,
                        denominator = %b.to_canonical_string(),
                        denominator_diagram_id = %other.diagram_id,
                        "Sample denominator evaluation"
                    );
                    let ratio = analyze_diff_and_sum(a.as_view(), b.as_view());
                    if let Some(a) = &ratio {
                        debug!(
                            sample_idx = %idx,
                            ratio = %a.to_canonical_string(),
                            "Sign-only comparison result for sample"
                        );
                    } else {
                        debug!(
                            sample_idx = %idx,
                            result = "none_ratio",
                            reason = "expressions_not_identical_up_to_sign",
                            "Sign-only comparison result for sample is None"
                        );
                    }

                    ratio
                })
                .collect::<HashSet<_>>();

            debug!(
                unique_ratios_found = %ratios.len(),
                "Found unique ratios from sample evaluations"
            );

            if ratios.len() == 1 {
                if let Some(ratio) = ratios.iter().next().unwrap().to_owned() {
                    debug!(
                        ratio = %ratio.to_canonical_string(),
                        method = "numerical_evaluation",
                        "Successfully matched diagrams using numerical evaluation"
                    );
                    return Some(ratio);
                } else {
                    debug!(
                        result = "none_ratio",
                        reason = "likely_zeros",
                        "Sample evaluations yielded None ratio"
                    );
                }
            } else {
                debug!(
                    result = "inconsistent_ratios",
                    unique_ratios_count = %ratios.len(),
                    "Sample evaluations yielded inconsistent ratios - cannot group diagrams"
                );
            }
        } else {
            debug!(
                comparison_type = "numerical_samples",
                result = "skipped",
                reason = "no_sample_evaluations",
                "Skipping numerical comparison"
            );
        }

        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            final_result = "no_grouping_possible",
            "No grouping possible between diagrams"
        );
        None
    }

    /// Compare two numerators allowing arbitrary scalar rescaling.
    #[instrument(skip_all, fields(self_id = %self.diagram_id, other_id = %other.diagram_id))]
    pub(crate) fn compare_with_scalar_rescaling(
        &self,
        other: &ProcessedNumeratorForComparison,
    ) -> Option<Atom> {
        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            "Starting scalar rescaling comparison between diagrams"
        );
        fn analyze_ratios(ratios: &HashSet<Option<Atom>>) -> Option<Atom> {
            if ratios.len() > 1 {
                None
            } else {
                let ratio = ratios.iter().next().unwrap().to_owned();
                if let Some(r) = ratio {
                    for head in LibraryRep::all_self_duals()
                        .chain(LibraryRep::all_inline_metrics())
                        .chain(LibraryRep::all_dualizables())
                        .map(|a| a.to_symbolic([W_.a__]).to_pattern())
                    {
                        if r.pattern_match(&head, None, None).next().is_some() {
                            return None;
                        }
                    }
                    Some(r)
                } else {
                    None
                }
            }
        }

        // Try canonized numerator comparison first
        if let (Some(canonized_num_a), Some(canonized_num_b)) = (
            self.canonized_numerator.as_ref(),
            other.canonized_numerator.as_ref(),
        ) {
            debug!(
                numerator = %canonized_num_a.to_canonical_string(),
                numerator_diagram_id = %self.diagram_id,
                denominator = %canonized_num_b.to_canonical_string(),
                denominator_diagram_id = %other.diagram_id,
                "Attempting symbolic comparison using canonized numerators",
            );
            let mut ratios = HashSet::<Option<Atom>>::default();
            let r = if canonized_num_a == canonized_num_b {
                debug!("Canonized numerators are identical");
                Some(Atom::num(1))
            } else if *canonized_num_a == canonized_num_b * Atom::num(-1) {
                debug!("Canonized numerators differ by sign only");
                Some(Atom::num(-1))
            } else if canonized_num_b.is_zero() {
                debug!("Reference canonized numerator is zero - cannot compute ratio");
                None
            } else {
                let ratio = canonized_num_a / canonized_num_b;
                debug!(
                    computed_ratio = %ratio.to_canonical_string(),
                    "Computed canonized numerator ratio"
                );
                Some(ratio)
            };
            ratios.insert(r);
            if let Some(ratio) = analyze_ratios(&ratios) {
                debug!(
                    ratio = %ratio.to_canonical_string(),
                    method = "canonized_numerators",
                    "Successfully matched diagrams using canonized numerators"
                );
                return Some(ratio);
            } else {
                debug!(
                    comparison_type = "canonized_numerator",
                    result = "rejected",
                    reason = "problematic_patterns",
                    rejection_stage = "analyze_ratios",
                    "Canonized numerator ratio was rejected"
                );
            }
        } else {
            debug!(
                comparison_type = "canonized_numerator",
                result = "skipped",
                reason = "missing_numerators",
                "Skipping canonized numerator comparison"
            );
        }

        // Fall back to sample evaluations
        if !self.sample_evaluations.is_empty() {
            debug!(
                comparison_type = "numerical_samples",
                sample_count = %self.sample_evaluations.len(),
                "Attempting numerical comparison using sample evaluations"
            );
            let evaluations_a = &self.sample_evaluations;
            let evaluations_b = &other.sample_evaluations;
            if evaluations_a.is_empty() {
                debug!(
                    result = "cannot_proceed",
                    reason = "empty_self_evaluations",
                    "Self sample evaluations are empty"
                );
                return None;
            }

            let ratios = evaluations_a
                .iter()
                .zip(evaluations_b.iter())
                .enumerate()
                .map(|(idx, (a, b))| {
                    debug!(
                        sample_id= %idx,
                        numerator = %a.to_canonical_string(),
                        numerator_diagram_id = %self.diagram_id,
                        denominator = %b.to_canonical_string(),
                        denominator_diagram_id = %other.diagram_id,
                        "Sample evaluation A"
                    );
                    if a == b {
                        Some(Atom::num(1))
                    } else if *a == b * Atom::num(-1) {
                        Some(Atom::num(-1))
                    } else if b.is_zero() || a.is_zero() {
                        debug!(
                            sample_idx = %idx,
                            "Skipping sample due to zero value"
                        );
                        None
                    } else {
                        let ratio = if ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL {
                            // let a_poly = a.to_polynomial(&Q_I.clone(), None);
                            // let b_poly = b.to_polynomial(&Q_I.clone(), None);
                            let a_poly = &self.sample_evaluations_as_polynomial[idx];
                            let b_poly = &other.sample_evaluations_as_polynomial[idx];
                            if a_poly.is_zero() || b_poly.is_zero() {
                                debug!(
                                    sample_idx = %idx,
                                    "Skipping sample due to zero value after expansion"
                                );
                                None
                            } else {
                                let element = COMPLEXRATPOLYFIELD.to_element(a_poly.clone(), b_poly.clone(), true);
                                Some(polyrat_to_atom(&element))
                            }
                            // let element = COMPLEXRATPOLYFIELD.to_element(a.to_polynomial(&Q_I, None), b.to_polynomial(&Q_I, None), true);
                            // Some(polyrat_to_atom(&element))
                        } else {
                            Some((a / b).cancel())
                        };
                        debug!(
                            sample_idx = %idx,
                            computed_ratio = %ratio.clone().map(|r| r.to_ordered_simple().to_string()).unwrap_or("None".into()),
                            "Computed ratio for sample"
                        );
                        ratio
                    }
                })
                .collect::<HashSet<_>>();

            debug!(
                unique_ratios_found = %ratios.len(),
                "Found unique ratios from sample evaluations"
            );
            if event_enabled!(tracing::Level::DEBUG, parent: &Span::current()) {
                for ((rat, a), b) in ratios.iter().zip(evaluations_a).zip(evaluations_b) {
                    debug!(
                        self_diagram_id = %self.diagram_id,
                        other_diagram_id = %other.diagram_id,
                        ratio_value = ?rat.as_ref().map(|ra| ra.floatify(13).to_ordered_simple()).unwrap_or("None".into()),
                        numerator_value = %a.floatify(13).to_ordered_simple(),
                        denominator_value = %b.floatify(13).to_ordered_simple(),
                        "Detailed sample evaluation ratio information"
                    );
                }
            }
            if let Some(ratio) = analyze_ratios(&ratios) {
                debug!(
                    ratio = %ratio.to_canonical_string(),
                    method = "numerical_evaluation",
                    "Successfully matched diagrams using numerical evaluation"
                );
                return Some(ratio);
            } else {
                debug!(
                    comparison_type = "numerical_samples",
                    result = "rejected",
                    reason = "inconsistent_or_problematic_patterns",
                    rejection_stage = "analyze_ratios",
                    "Sample evaluation ratios were rejected"
                );
            }
        } else {
            debug!(
                comparison_type = "numerical_samples",
                result = "skipped",
                reason = "no_sample_evaluations",
                "Skipping numerical comparison"
            );
        }

        debug!(
            self_diagram_id = %self.diagram_id,
            other_diagram_id = %other.diagram_id,
            final_result = "no_grouping_possible",
            "No grouping possible between diagrams"
        );
        None
    }
    #[instrument(skip_all,fields(graph_name = %graph.name))]
    pub(crate) fn from_numerator_symbolic_expression(
        diagram_id: usize,
        graph: &Graph,
        mut numerator: Atom,
        samples: &[NumeratorSample],
        settings: &GlobalSettings,
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
    ) -> Result<Self, FeynGenError> {
        let default_processed_data = ProcessedNumeratorForComparison {
            diagram_id,
            canonized_numerator: None,
            sample_evaluations: vec![],
            sample_evaluations_as_polynomial: vec![],
            sample_evaluations_are_zero: vec![],
        };
        let res = if let Some(group_options) = numerator_aware_isomorphism_grouping.get_options() {
            if group_options.test_canonized_numerator
                || group_options.number_of_numerical_samples > 0
            {
                let lmb_reps = graph.integrand_replacement(
                    &graph.full_filter(),
                    &graph.loop_momentum_basis,
                    &[W_.x___],
                );

                numerator = numerator.replace_multiple(&lmb_reps);

                debug!(numerator=%numerator.to_ordered_simple(),diagram_id=%diagram_id,debug_dot=%graph.debug_dot(),"Initial Numerator");

                let canonized_numerator = if group_options.test_canonized_numerator {
                    let mut canonized_numerator_to_consider = numerator.canonize(Aind::Dummy);
                    if group_options.fully_numerical_substitution_when_comparing_numerators {
                        canonized_numerator_to_consider =
                            ProcessDefinition::substitute_color_factors(
                                canonized_numerator_to_consider.as_atom_view(),
                            )
                    };
                    // IMPORTANT: we must make sure to collect all common coefficients first
                    // with `collect_factors()` to ensure that common factors get simplified when looking at ratios.
                    if matches!(
                        numerator_aware_isomorphism_grouping,
                        NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(
                            _
                        )
                    ) {
                        canonized_numerator_to_consider =
                            canonized_numerator_to_consider.collect_factors();
                    }
                    // numerator = canonized_numerator_to_consider.clone();
                    Some(canonized_numerator_to_consider)
                } else {
                    None
                };

                let mut numerators = vec![numerator];
                let sample_parse_settings = ParseSettings {
                    shorthand_parsing: ShorthandParsing::expand_schoonschip_only(),
                    ..Default::default()
                };

                if settings
                    .generation
                    .feyngen
                    .gamma_simplification_closure_check
                {
                    debug!(numerator = %numerators[0].to_ordered_simple(),"Gamma Simplifying");
                    let gamma_simplified_numerator = numerators[0].simplify_gamma();
                    crate::debug_tags!(#generation, #profile, #graph, #inspect, #dump;
                        stage = "feyngen_closure_check_after_simplify_gamma",
                        log.after_gamma = gamma_simplified_numerator,
                        "Feyngen closure check after gamma simplification"
                    );
                    numerators.push(gamma_simplified_numerator);
                    debug!("Done Simplifying");
                }

                let mut sample_evals = VecDeque::new();
                for numerator in numerators {
                    if is_interrupted() {
                        return Err(FeynGenError::Interrupted);
                    }
                    let expanded = numerator.expand_color();

                    let sample_evaluations = samples
                    .iter()
                    .enumerate()
                    .map(|(sample_index, (reps, lib))| {
                        if is_interrupted() {
                            return Err(FeynGenError::Interrupted);
                        }
                        let mut sample_evaluation = expanded
                            .iter()
                            .enumerate()
                            .try_fold(Atom::Zero, |acc, (term_index, (c, l))| {
                                if is_interrupted() {
                                    return Err(FeynGenError::Interrupted);
                                }
                                let context = || {
                                    format!(
                                        "while evaluating numerator sample for graph '{}' (diagram id {}, sample {}, expanded color term {}, color_len={}, lorentz_len={})",
                                        graph.name,
                                        diagram_id,
                                        sample_index,
                                        term_index,
                                        c.to_plain_string().len(),
                                        l.to_plain_string().len()
                                    )
                                };
                                debug!("Sample evaluation inputs c:{c},l:{l}");
                                let mut net = ParamParsingNet::try_from_view(
                                    l.as_view(),
                                    lib,
                                    &sample_parse_settings,
                                )
                                .map_err(|source| {
                                    FeynGenError::Eyre(eyre!(source).wrap_err(context()))
                                })?;
                                net.store
                                    .scalar
                                    .iter_mut()
                                    .for_each(|a| *a = a.replace_multiple(reps));

                                // debug!(net=?net.dot_pretty());
                                net.execute::<Sequential, MinResultRank, _, _, _>(
                                    lib,
                                    PARAM_FUN_LIB.deref(),
                                )
                                .map_err(|source| {
                                    FeynGenError::Eyre(eyre!(source).wrap_err(format!(
                                        "{}; failed tensor-network execution for net:\n{}",
                                        context(),
                                        net.dot_pretty()
                                    )))
                                })?;

                                // let c = ProcessDefinition::substitute_color_factors(c.as_view())
                                //     .expand();

                                let scalar = net
                                    .result_scalar();

                                if scalar.is_err(){
                                    debug!("Failed scalar for {} for graph:{}",net.dot_pretty(),graph.debug_dot());
                                }

                                let scalar:Atom = scalar
                                    .map_err(|source| {
                                        FeynGenError::Eyre(
                                            eyre!(source).wrap_err(format!(
                                                "{}; expected scalar from net:\n{}\nfor graph:\n{}",
                                                context(),
                                                net.dot_pretty(),
                                                graph.debug_dot()
                                            )),
                                        )
                                    })?
                                    .into();

                                // println!("Trying to canonize:{c}");
                                let canonized_color = c.canonize::<Aind>(Aind::Dummy);
                                debug!("canonizing \n{c}\n gives\n{canonized_color}");
                                let a = ProcessDefinition::substitute_color_factors(
                                    (canonized_color * scalar).as_view(),
                                );

                                debug!(evaluated=%a.printer(LOGPRINTOPTS.clone()),"evaluated{}",a.floatify(13).printer(LOGPRINTOPTS.clone()));
                                Ok::<_, FeynGenError>(acc + a)
                            })?;

                        if EXPAND_NUMERICAL_SAMPLES_BEFORE_COMPARISON {
                            sample_evaluation = sample_evaluation.expand();
                        }
                        // TODO: optimize the above by instead directly storing "a.to_polynomial(&Q_I.clone(), None)" in the sample record, and not the expanded atom.
                        // Only do that for the if-branch below though.
                        let sample_evaluation = if ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL || !matches!(numerator_aware_isomorphism_grouping,NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(_))
                        {
                            sample_evaluation
                        } else {
                            // When not looking at the ratio of samples as a rational polynomial, we must make sure to collect all common coefficients first
                            // with `collect_factors()` to ensure that common factors get simplified when looking at ratios.
                            sample_evaluation.collect_factors()
                        };
                        Ok(sample_evaluation)
                    })
                    .collect::<Result<Vec<_>, FeynGenError>>()?;

                    sample_evals.push_back(sample_evaluations);
                }

                let sample_evaluations = sample_evals.pop_front().unwrap();
                let sample_gamma_simplified_evaluations =
                    sample_evals.pop_front().unwrap_or(vec![]);

                if settings
                    .generation
                    .feyngen
                    .gamma_simplification_closure_check
                {
                    for (i, (a, b)) in sample_evaluations
                        .iter()
                        .zip(sample_gamma_simplified_evaluations.iter())
                        .enumerate()
                    {
                        if !(a - b).expand().is_zero() {
                            return Err(FeynGenError::Eyre(eyre!(
                                "Gamma simplification closure check failed for diagram ID {} on sample evaluation #{}. Numerator evaluation before gamma simplification: {}. After gamma simplification: {}.",
                                diagram_id,
                                i,
                                a.to_ordered_simple(),
                                b.to_ordered_simple()
                            )));
                        }
                    }
                }
                let samples_evals_as_polynomial: SampleEvaluationsAsPolynomial =
                    if ANALYZE_RATIO_AS_RATIONAL_POLYNOMIAL
                        || !matches!(
                        numerator_aware_isomorphism_grouping,
                        NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(
                            _
                        )
                    ) {
                        let se_as_poly = sample_evaluations
                            .iter()
                            .map(|a| a.to_polynomial(&Q_I.clone(), None))
                            .collect::<Vec<_>>();
                        let se_are_zero =
                            se_as_poly.iter().map(|p| p.is_zero()).collect::<Vec<_>>();
                        (se_as_poly, se_are_zero)
                    } else {
                        (
                            vec![],
                            sample_evaluations
                                .iter()
                                .map(|a| a.expand().is_zero())
                                .collect::<Vec<_>>(),
                        )
                    };
                ProcessedNumeratorForComparison {
                    diagram_id,
                    canonized_numerator,
                    sample_evaluations,
                    sample_evaluations_as_polynomial: samples_evals_as_polynomial.0,
                    sample_evaluations_are_zero: samples_evals_as_polynomial.1,
                }
            } else {
                default_processed_data
            }
        } else {
            default_processed_data
        };
        Ok(res)
    }
}

// #[test]
// pub(crate) fn symbolica_symm_factors_bug() {
//     let external_edges_for_generation: Vec<(usize, (Option<bool>, usize))> = vec![];

//     let vertex_signatures_for_generation_a = vec![
//         vec![(None, 21), (None, 21), (None, 21)],
//         vec![(None, 21), (None, 21), (None, 21), (None, 21)],
//         vec![(Some(true), 1), (Some(false), 1), (None, 21)],
//         vec![(Some(true), 6), (Some(false), 6), (None, 21)],
//         vec![(Some(false), 9000005), (Some(true), 9000005), (None, 21)],
//         vec![(Some(true), 2), (Some(false), 2), (None, 21)],
//     ];

//     let settings = GenerationSettings::new()
//         .max_loops(5)
//         .max_bridges(0)
//         .allow_self_loops(true);

//     //let mut graphs_a = SymbolicaGraph::generate(
//     //    &external_edges_for_generation,
//     //    vertex_signatures_for_generation_a.as_slice(),
//     //    None,
//     //    Some(5),
//     //    Some(0),
//     //    true,
//     //);

//     let mut graphs_a = SymbolicaGraph::generate(
//         &external_edges_for_generation,
//         vertex_signatures_for_generation_a.as_slice(),
//         &settings,
//     )
//     .unwrap();

//     graphs_a.retain(|g, _| g.num_loops() >= 5);

//     let mut tot_symm_fact_graphs_a = Atom::num(0);
//     for (_g, s) in graphs_a.iter() {
//         tot_symm_fact_graphs_a = tot_symm_fact_graphs_a + Atom::num(1) / Atom::num(s.clone());
//     }
//     println!("tot_symm_fact_graphs_A = {}", tot_symm_fact_graphs_a);

//     let vertex_signatures_for_generation_b = vec![
//         vec![(Some(true), 1), (Some(false), 1), (None, 21)],
//         vec![(None, 21), (None, 21), (None, 21), (None, 21)],
//         vec![(Some(true), 6), (Some(false), 6), (None, 21)],
//         vec![(None, 21), (None, 21), (None, 21)],
//         vec![(Some(false), 9000005), (Some(true), 9000005), (None, 21)],
//         vec![(Some(true), 2), (Some(false), 2), (None, 21)],
//     ];
//     let mut graphs_b = SymbolicaGraph::generate(
//         &external_edges_for_generation,
//         vertex_signatures_for_generation_b.as_slice(),
//         &settings,
//     )
//     .unwrap();

//     graphs_b.retain(|g, _| g.num_loops() >= 5);

//     let mut tot_symm_fact_graphs_b = Atom::num(0);
//     for (_g, s) in graphs_b.iter() {
//         tot_symm_fact_graphs_b = tot_symm_fact_graphs_b + Atom::num(1) / Atom::num(s.clone());
//     }
//     println!("tot_symm_fact_graphs_B = {}", tot_symm_fact_graphs_b);

//     assert!(tot_symm_fact_graphs_a == tot_symm_fact_graphs_b);
// }
