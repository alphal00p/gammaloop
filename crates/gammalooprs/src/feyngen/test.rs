use std::sync::Arc;

use ahash::HashMap;
use ahash::HashMapExt;
use ahash::HashSet;
use idenso::IndexTooling;
use idenso::gamma::AGS;
use idenso::metric::PermuteWithMetric;
use insta::assert_snapshot;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeData;
use linnet::half_edge::involution::Flow;
use linnet::half_edge::involution::Orientation;
use spenso::network::library::LibraryTensor;
use spenso::network::library::symbolic::ExplicitKey;
use spenso::network::library::symbolic::TensorLibrary;
use spenso::structure::PermutedStructure;
use spenso::structure::representation::Minkowski;
use spenso::structure::representation::RepName;
use spenso::tensors::data::DataTensor;
use spenso::tensors::parametric::ParamTensor;
use spenso_hep_lib::gamma_data_weyl;
use spenso_hep_lib::gamma5_weyl_data;
use spenso_hep_lib::proj_m_data_weyl;
use spenso_hep_lib::proj_p_data_weyl;
use symbolica::atom::Atom;
use symbolica::atom::AtomCore;
use symbolica::coefficient::Coefficient;
use symbolica::graph::Graph as SymbolicaGraph;
use symbolica::id::Replacement;
use symbolica::parse_lit;
use tracing::debug;

use crate::dot;
use crate::feyngen::diagram_generator::EdgeColor;
use crate::graph::FeynmanGraph;
use crate::initialisation::test_initialise;
use crate::numerator::aind::Aind;
use crate::numerator::graph::ReversibleEdge;
use crate::settings::GlobalSettings;
use crate::settings::global::Parallelisation;

use crate::utils::GS;
use crate::uv::UltravioletGraph;
// use crate::graph::BareGraph;
use super::GenerationType;
use super::GraphGroupingOptions;
use super::NumeratorAwareGraphGroupingOption;
use super::SewedFilterOptions;
use super::SnailFilterOptions;
use super::TadpolesFilterOptions;
use super::diagram_generator::NodeColorWithVertexRule;
use super::diagram_generator::ProcessedNumeratorForComparison;
use super::{FeynGenFilter, FeynGenFilters};
use crate::graph::{Graph, parse::IntoGraph};
use crate::model::ArcVertexRule;
use crate::model::ColorStructure;
use crate::model::Model;
use crate::model::VertexRule;
use crate::numerator::GlobalPrefactor;
use crate::processes::ProcessDefinition;
use crate::utils::load_generic_model;

fn manual_lib<C: Into<Coefficient>>(
    loop_momenta: Vec<Vec<C>>,
    pol_v: Vec<(isize, Vec<C>, Vec<C>)>,
    pol_out: Vec<(isize, Vec<C>)>,
    model: &Model,
) -> TensorLibrary<ParamTensor<ExplicitKey<Aind>>, Aind> {
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

    let mut lib = weyl; //hep_lib(F(1.), F(0.));
    for (i, v) in loop_momenta.into_iter().enumerate() {
        let key = ExplicitKey::from_iter(
            [Minkowski {}.new_rep(4)],
            GS.loop_mom,
            Some(vec![Atom::num(i)]),
        );

        debug!("lib_loop:{}", key.clone().permute_with_metric());
        let key =
            ParamTensor::from_dense(key.structure, v.into_iter().map(|n| Atom::num(n)).collect())
                .unwrap();

        lib.insert_explicit(PermutedStructure::identity(key));
    }

    for (i, (pdg, pol, ext_mom)) in pol_v.into_iter().enumerate() {
        let additional_args = Some(vec![Atom::num(i)]);
        let key = ExplicitKey::from_iter(
            [Minkowski {}.new_rep(4)],
            GS.external_mom,
            additional_args.clone(),
        );

        debug!("lib_ext:{}", key.clone().permute_with_metric());

        let key =
            ParamTensor::from_dense(key.structure, ext_mom.into_iter().map(Atom::num).collect())
                .unwrap();

        lib.insert_explicit(PermutedStructure::identity(key));

        let p = model.get_particle_from_pdg(pdg);

        let structure = p.spin_reps();
        let global_name = EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Sink);
        let key = PermutedStructure::identity(ExplicitKey {
            structure,
            global_name,
            additional_args,
        });

        debug!("lib_pol:{}", key.clone().permute_with_metric());
        let key = ParamTensor::from_dense(key.structure, pol.into_iter().map(Atom::num).collect())
            .unwrap();

        lib.insert_explicit(PermutedStructure::identity(key));
    }

    for (i, (pdg, pol)) in pol_out.into_iter().enumerate() {
        let additional_args = Some(vec![Atom::num(i)]);

        let p = model.get_particle_from_pdg(pdg);

        let structure = p.spin_reps();
        let global_name = EdgeData::new(p, Orientation::Default).pol_symbol(Flow::Source);

        let key = PermutedStructure::identity(ExplicitKey {
            structure,
            global_name,
            additional_args,
        });
        debug!("lib_pol:{}", key.clone().permute_with_metric());
        let key = ParamTensor::from_dense(key.structure, pol.into_iter().map(Atom::num).collect())
            .unwrap();

        lib.insert_explicit(PermutedStructure::identity(key));
    }
    lib
}

#[test]
fn gl_11_vs_gl_12() {
    test_initialise().unwrap();

    let _gl_12:Graph = dot!(
        digraph GL12{
            num = "1";
        overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
        // projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
        projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵ(0,spenso::mink(4,hedge(0)))";
        // projector = "spenso::g(spenso::mink(4,hedge(1)),spenso::mink(4,hedge(0)))";
        2:1-> 3:0 [id=0  is_cut=0   particle="a"];
        1:8-> 2:9 [id=1 dir=back particle="d~"];
        0:2-> 1:3 [id=2 dir=back lmb_id=0 particle="d~"];
        0:4-> 1:5 [id=3 dir=none  particle="Z"];
        3:6-> 0:7 [id=4 dir=back lmb_id=1 particle="d~"];
        2:10-> 3:11 [id=5 dir=back particle="d~"];
        }
    ).unwrap();

    let _gl_11:Graph = dot!(
        digraph GL11{
            num = "1";
        overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
        // projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
        projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵ(0,spenso::mink(4,hedge(0)))";
        // projector = "spenso::g(spenso::mink(4,hedge(1)),spenso::mink(4,hedge(0)))";
        2:1-> 3:0 [id=0  is_cut=0   particle="a"];
        1:8-> 2:9 [id=1 particle="d"];
        0:2-> 1:3 [id=2 lmb_id=0 particle="d"];
        0:4-> 1:5 [id=3 dir=none  particle="Z"];
        3:6-> 0:7 [id=4 lmb_id=1 particle="d"];
        2:10-> 3:11 [id=5 particle="d"];
        }
    ).unwrap();

    let gl_12:Graph = dot!(
        digraph GL12{
            num = "1";
        overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
        // projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
        projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵ(0,spenso::mink(4,hedge(0)))";
        // projector = "spenso::g(spenso::mink(4,hedge(1)),spenso::mink(4,hedge(0)))";
        0[dod=0 int_id=V_79 num_old="(UFO::GC_58*(-2*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projp(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(2)))+spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(2))))+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(2))))*spenso::g(spenso::dind(spenso::cof(3,hedge(6))),spenso::cof(3,hedge(2)))"];
        1[dod=0 int_id=V_79 num_old="(UFO::GC_58*(-2*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projp(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(8)))+spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(8))))+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(8))))*spenso::g(spenso::dind(spenso::cof(3,hedge(3))),spenso::cof(3,hedge(8)))"];
        2[dod=0 int_id=V_71 num_old="UFO::GC_1*spenso::g(spenso::dind(spenso::cof(3,hedge(9))),spenso::cof(3,hedge(10)))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(10)),spenso::mink(4,hedge(0)))"];
        3[dod=0 int_id=V_71 num_old="UFO::GC_1*spenso::g(spenso::dind(spenso::cof(3,hedge(11))),spenso::cof(3,hedge(7)))*spenso::gamma(spenso::bis(4,hedge(11)),spenso::bis(4,hedge(7)),spenso::mink(4,hedge(1)))"];
        2:1-> 3:0 [id=0 dir=none source=2 sink=2  dod=-2 is_cut=0 is_dummy=false lmb_rep="P(0,a___)" name=e0 num_old="1" particle="a"];
        1:8-> 2:9 [id=1 dir=back source=1 sink=0  dod=-1 is_dummy=false lmb_rep="K(0,a___)+K(1,a___)" name=e1 num_old="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(8))),spenso::cof(3,hedge(9)))*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(9)),spenso::mink(4,edge(1,1)))" particle="d~"];
        0:2-> 1:3 [id=2 dir=back source=1 sink=0  dod=-1 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e2 num_old="Q(2,spenso::mink(4,edge(2,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::cof(3,hedge(3)))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3)),spenso::mink(4,edge(2,1)))" particle="d~"];
        0:4-> 1:5 [id=3 dir=none source=2 sink=2  dod=0 is_dummy=false  lmb_rep="K(1,a___)" name=e3 num_old="-spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))+UFO::MZ^-2*Q(3,spenso::mink(4,hedge(4)))*Q(3,spenso::mink(4,hedge(5)))" particle="Z"];
        0:6-> 3:7 [id=4 source=0 sink=1  dod=-1 is_dummy=false lmb_id=1 lmb_rep="-K(0,a___)-K(1,a___)" name=e4 num_old="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(7))),spenso::cof(3,hedge(6)))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,edge(4,1)))" particle="d"];
        2:10-> 3:11 [id=5 dir=back source=1 sink=0  dod=-1 is_dummy=false lmb_rep="-P(0,a___)+K(0,a___)+K(1,a___)" name=e5 num_old="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(10))),spenso::cof(3,hedge(11)))*spenso::gamma(spenso::bis(4,hedge(10)),spenso::bis(4,hedge(11)),spenso::mink(4,edge(5,1)))" particle="d~"];
        }
    ).unwrap();

    let gl_11:Graph = dot!(
        digraph GL11{
            num = "1";
            overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
            // projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
            projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵ(0,spenso::mink(4,hedge(0)))";
            // projector = "spenso::g(spenso::mink(4,hedge(1)),spenso::mink(4,hedge(0)))";
            0[dod=0 int_id=V_79 num_old="(UFO::GC_58*(-2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projp(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6)))+spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6))))+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6))))*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::cof(3,hedge(6)))"];
            1[dod=0 int_id=V_79 num_old="(UFO::GC_58*(-2*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projp(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3)))+spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3))))+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3))))*spenso::g(spenso::dind(spenso::cof(3,hedge(8))),spenso::cof(3,hedge(3)))"];
            2[dod=0 int_id=V_71 num_old="UFO::GC_1*spenso::g(spenso::dind(spenso::cof(3,hedge(10))),spenso::cof(3,hedge(9)))*spenso::gamma(spenso::bis(4,hedge(10)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(0)))"];
            3[dod=0 int_id=V_71 num_old="UFO::GC_1*spenso::g(spenso::dind(spenso::cof(3,hedge(7))),spenso::cof(3,hedge(11)))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(11)),spenso::mink(4,hedge(1)))"];
            2:1-> 3:0 [id=0 dir=none source=2 sink=2  dod=-2 is_cut=0 is_dummy=false lmb_rep="P(0,a___)" name=e0 num_old="1" particle="a"];
            1:8-> 2:9 [id=1 source=0 sink=1  dod=-1 is_dummy=false lmb_rep="K(0,a___)+K(1,a___)" name=e1 num_old="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(9))),spenso::cof(3,hedge(8)))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(1,1)))" particle="d"];
            0:2-> 1:3 [id=2 source=0 sink=1  dod=-1 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e2 num_old="Q(2,spenso::mink(4,edge(2,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(3))),spenso::cof(3,hedge(2)))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(2,1)))" particle="d"];
            0:4-> 1:5 [id=3 dir=none source=2 sink=2  dod=0 is_dummy=false  lmb_rep="K(1,a___)" name=e3 num_old="-spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))+UFO::MZ^-2*Q(3,spenso::mink(4,hedge(4)))*Q(3,spenso::mink(4,hedge(5)))" particle="Z"];
            0:6-> 3:7 [id=4 dir=back source=1 sink=0  dod=-1 is_dummy=false lmb_id=1 lmb_rep="-K(0,a___)-K(1,a___)" name=e4 num_old="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(6))),spenso::cof(3,hedge(7)))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(7)),spenso::mink(4,edge(4,1)))" particle="d~"];
            2:10-> 3:11 [id=5 source=0 sink=1  dod=-1 is_dummy=false lmb_rep="-P(0,a___)+K(0,a___)+K(1,a___)" name=e5 num_old="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(11))),spenso::cof(3,hedge(10)))*spenso::gamma(spenso::bis(4,hedge(11)),spenso::bis(4,hedge(10)),spenso::mink(4,edge(5,1)))" particle="d"];
        }
        ).unwrap();

    let loop_momenta = vec![vec![79, 83, 89, 97], vec![101, 103, 107, 109]];
    let pol_v = vec![(22, vec![43, 47, 53, 59], vec![29, 31, 37, 41])];
    let pol_out = vec![(22, vec![61, 67, 71, 73])];
    let model = load_generic_model("sm");

    let _reps = vec![
        Replacement::new(parse_lit!(UFO::MZ).to_pattern(), Atom::num(13)),
        Replacement::new(parse_lit!(UFO::ee).to_pattern(), Atom::num(17)),
        Replacement::new(parse_lit!(UFO::cw).to_pattern(), Atom::num(19)),
        Replacement::new(parse_lit!(UFO::sw).to_pattern(), Atom::num(23)),
    ];

    let lib = manual_lib(loop_momenta, pol_v, pol_out, &model);
    let mut cpl_reps: Vec<Replacement> = vec![];
    for cpl in model.couplings.values() {
        let [lhs, rhs] = cpl.rep_rule();
        cpl_reps.push(Replacement::new(lhs.to_pattern(), rhs));
    }

    let mut numerator_11 = gl_11.numerator(&gl_11.no_dummy(), &gl_11.empty_subgraph());
    numerator_11.state.expr *= &gl_11.global_prefactor.num * &gl_11.global_prefactor.projector; // * &bare_graph.overall_factor;
    numerator_11.state.expr = numerator_11.state.expr.replace_multiple(&cpl_reps);

    let numerator_color_simplified_11 = numerator_11
        .clone()
        .color_simplify()
        .get_single_atom()
        .unwrap();
    // .replace_multiple(&reps);
    // .replace(parse_lit!(UFO::MZ ^ (-2)))
    // .with(Atom::Zero)
    // // .replace(parse_lit!(UFO::cw))
    // // .with(parse_lit!(UFO::sw));
    // // .replace(parse_lit!(UFO::sw))
    // // .with(Atom::num(1));
    // .replace(parse_lit!(spenso::projp(w__)))
    // .with(parse_lit!(-spenso::projm(w__)));
    // .replace(parse_lit!(spenso::projm(w__)))
    // .with(parse_lit!(spenso::g(w__)));

    let mut numerator_12 = gl_12.numerator(&gl_12.no_dummy(), &gl_12.empty_subgraph());
    numerator_12.state.expr *= &gl_12.global_prefactor.num * &gl_12.global_prefactor.projector; // * &bare_graph.overall_factor;
    numerator_12.state.expr = numerator_12.state.expr.replace_multiple(&cpl_reps);

    let numerator_color_simplified_12 = numerator_12
        .clone()
        .color_simplify()
        .get_single_atom()
        .unwrap();
    // .replace_multiple(&reps);
    // .replace(parse_lit!(UFO::MZ ^ (-2)))
    // .with(Atom::Zero)
    // .replace(parse_lit!(spenso::projp(w__)))
    // .with(parse_lit!(-spenso::projm(w__)));
    // .replace(parse_lit!(UFO::cw))
    // .with(parse_lit!(UFO::sw));
    // .replace(parse_lit!(UFO::sw))
    // .with(Atom::num(1));
    // .replace(parse_lit!(spenso::projp(w__)))
    // .with(parse_lit!(spenso::g(w__)))
    // .replace(parse_lit!(spenso::projm(w__)))
    // .with(parse_lit!(spenso::g(w__)));

    println!("color_simplified_12 {}", numerator_color_simplified_12);
    println!("color_simplified_11 {}", numerator_color_simplified_11);
    println!(
        "color_simplified_12_can {}",
        numerator_color_simplified_12.canonize(Aind::Dummy)
    );
    println!(
        "color_simplified_11_can {}",
        numerator_color_simplified_11.canonize(Aind::Dummy)
    );

    let r = (numerator_color_simplified_11.canonize(Aind::Dummy)
        / &numerator_color_simplified_12.canonize(Aind::Dummy))
        .expand();

    println!("ratio:{r}");

    let samples = [(vec![], lib)];

    let options = NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(
        GraphGroupingOptions {
            test_canonized_numerator: false,
            ..Default::default()
        },
    );

    let pn_12 = ProcessedNumeratorForComparison::from_numerator_symbolic_expression(
        12,
        &gl_12,
        numerator_color_simplified_12,
        &samples,
        &GlobalSettings {
            n_cores: Parallelisation {
                feyngen: 10,
                ..Default::default()
            },
            ..Default::default()
        },
        &options,
    )
    .unwrap();

    let pn_11 = ProcessedNumeratorForComparison::from_numerator_symbolic_expression(
        11,
        &gl_11,
        numerator_color_simplified_11,
        &samples,
        &GlobalSettings {
            n_cores: Parallelisation {
                feyngen: 10,
                ..Default::default()
            },
            ..Default::default()
        },
        &options,
    )
    .unwrap();

    if let Some(a) = pn_11.compare_with_scalar_rescaling(&pn_12) {
        println!("{}", a);
        println!("Expanded: {}", a.expand());
    }
}

#[test]
fn cut_content() {
    let model = load_generic_model("sm");

    let mut coupling = HashMap::new();
    coupling.insert("QED".into(), (6, Some(6)));
    let mut pert = HashMap::new();
    pert.insert("QCD".into(), 1);
    let process_definition = ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![-11, 11],
        final_pdgs_lists: vec![vec![5, -5, 25]],
        loop_count_range: (3, 3),
        symmetrize_initial_states: true,
        symmetrize_final_states: true,
        symmetrize_left_right_states: true,
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 6,
        amplitude_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
        ]),
        cross_section_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
            FeynGenFilter::PerturbativeOrders(pert),
            FeynGenFilter::CouplingOrders(coupling),
            FeynGenFilter::LoopCountRange((3, 3)),
            FeynGenFilter::BlobRange(1..=1),
            FeynGenFilter::SpectatorRange(0..=0),
        ]),
        filter_self_loop: true,
        folder_name: "N/A".into(),
        process_id: 0,
        graph_prefix: "GL".into(),
        loop_momentum_bases: None,
        numerator_grouping: NumeratorAwareGraphGroupingOption::NoGrouping,
        prefactor: GlobalPrefactor::default(),
        selected_graphs: None,
        vetoed_graphs: None,
        filter_zero_flow_edges: true,
    };

    let mut graph = SymbolicaGraph::new();
    #[allow(non_snake_case)]
    let bbH = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_78"),
    };
    let bba = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_73"),
    };
    let bbg = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_76"),
    };
    let epema = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_98"),
    };

    let dummy_external_vertex_rule = ArcVertexRule(Arc::new(VertexRule {
        name: "external".into(),
        couplings: vec![],
        lorentz_structures: vec![],
        particles: vec![],
        color_structures: ColorStructure {
            color_structure: vec![],
        },
    }));

    let v0 = graph.add_node(bbH.clone());
    let v1 = graph.add_node(bbH.clone());
    let v2 = graph.add_node(bbg.clone());
    let v3 = graph.add_node(bba.clone());
    let v4 = graph.add_node(bbg.clone());
    let v5 = graph.add_node(bba.clone());
    let v6 = graph.add_node(epema.clone());
    let v7 = graph.add_node(epema.clone());
    let e1 = NodeColorWithVertexRule {
        external_tag: 1,
        vertex_rule: dummy_external_vertex_rule.clone(),
    };
    let mut e2 = e1.clone();
    e2.external_tag = 2;
    let mut e3 = e1.clone();
    e3.external_tag = 3;
    let mut e4 = e1.clone();
    e4.external_tag = 4;

    let v8 = graph.add_node(e1.clone());
    let v9 = graph.add_node(e2.clone());
    let v10 = graph.add_node(e3.clone());
    let v11 = graph.add_node(e4.clone());

    let h = EdgeColor::from_particle(model.get_particle("H"));

    let b = EdgeColor::from_particle(model.get_particle("b"));
    let bbar = EdgeColor::from_particle(model.get_particle("b").0.get_anti_particle(&model));
    let g = EdgeColor::from_particle(model.get_particle("g"));
    let a = EdgeColor::from_particle(model.get_particle("a"));
    let eplus = EdgeColor::from_particle(model.get_particle("e+"));
    let eminus = EdgeColor::from_particle(model.get_particle("e-"));
    graph.add_edge(v0, v1, true, h).unwrap();
    graph.add_edge(v0, v5, true, b).unwrap();
    graph.add_edge(v1, v3, true, b).unwrap();
    graph.add_edge(v2, v4, true, g).unwrap();
    graph.add_edge(v2, v4, true, b).unwrap();
    graph.add_edge(v3, v2, true, b).unwrap();
    graph.add_edge(v3, v7, true, a).unwrap();
    graph.add_edge(v4, v1, true, b).unwrap();
    graph.add_edge(v5, v0, true, b).unwrap();
    graph.add_edge(v5, v6, true, a).unwrap();
    graph.add_edge(v6, v10, false, eplus).unwrap();
    graph.add_edge(v6, v11, false, eminus).unwrap();
    graph.add_edge(v8, v7, false, eplus).unwrap();
    graph.add_edge(v9, v7, false, eminus).unwrap();

    let (n_unresolved, unresolved_type) = process_definition.unresolved_cut_content(&model);
    assert!(!process_definition.half_edge_filters(
        &model,
        &graph,
        &[],
        n_unresolved,
        &unresolved_type
    ));

    let mut double_double_triangle = SymbolicaGraph::new();
    let v0 = double_double_triangle.add_node(e1.clone());
    let v1 = double_double_triangle.add_node(e2.clone());
    let v2 = double_double_triangle.add_node(e3.clone());
    let v3 = double_double_triangle.add_node(e4.clone());

    let v4 = double_double_triangle.add_node(bba.clone());
    let v5 = double_double_triangle.add_node(bba.clone());
    let v6 = double_double_triangle.add_node(bbH.clone());
    let v7 = double_double_triangle.add_node(bbH.clone());
    let v8 = double_double_triangle.add_node(bbg.clone());
    let v9 = double_double_triangle.add_node(bbg.clone());
    let v10 = double_double_triangle.add_node(bbg.clone());
    let v11 = double_double_triangle.add_node(bbg.clone());
    let v12 = double_double_triangle.add_node(epema.clone());
    let v13 = double_double_triangle.add_node(epema.clone());

    double_double_triangle
        .add_edge(v0, v13, true, eplus)
        .unwrap();
    double_double_triangle
        .add_edge(v1, v13, true, eminus)
        .unwrap();
    double_double_triangle.add_edge(v4, v13, false, a).unwrap();
    double_double_triangle
        .add_edge(v4, v11, true, bbar)
        .unwrap();
    double_double_triangle.add_edge(v4, v10, true, b).unwrap();
    double_double_triangle.add_edge(v10, v11, false, g).unwrap();
    double_double_triangle.add_edge(v6, v11, true, b).unwrap();
    double_double_triangle
        .add_edge(v6, v10, true, bbar)
        .unwrap();
    double_double_triangle.add_edge(v6, v7, true, h).unwrap();
    double_double_triangle.add_edge(v7, v8, true, bbar).unwrap();
    double_double_triangle.add_edge(v7, v9, true, b).unwrap();
    double_double_triangle.add_edge(v5, v8, true, b).unwrap();
    double_double_triangle.add_edge(v5, v9, true, bbar).unwrap();
    double_double_triangle.add_edge(v8, v9, true, g).unwrap();
    double_double_triangle.add_edge(v5, v12, true, a).unwrap();
    double_double_triangle
        .add_edge(v12, v2, true, eminus)
        .unwrap();
    double_double_triangle
        .add_edge(v12, v3, true, eplus)
        .unwrap();

    let mut coupling = HashMap::new();
    coupling.insert("QED".into(), (6, Some(6)));
    let mut pert = HashMap::new();
    pert.insert("QCD".into(), 2);
    let process_definition = ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![-11, 11],
        final_pdgs_lists: vec![vec![5, -5, 25]],
        loop_count_range: (4, 4),
        symmetrize_initial_states: true,
        symmetrize_final_states: true,
        symmetrize_left_right_states: true,
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 6,
        amplitude_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
        ]),
        cross_section_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
            FeynGenFilter::PerturbativeOrders(pert),
            FeynGenFilter::CouplingOrders(coupling),
            FeynGenFilter::LoopCountRange((4, 4)),
            FeynGenFilter::BlobRange(1..=1),
            FeynGenFilter::SpectatorRange(0..=0),
        ]),
        filter_self_loop: true,
        folder_name: "N/A".into(),
        process_id: 0,
        graph_prefix: "GL".into(),
        loop_momentum_bases: None,
        numerator_grouping: NumeratorAwareGraphGroupingOption::NoGrouping,
        prefactor: GlobalPrefactor::default(),
        selected_graphs: None,
        vetoed_graphs: None,
        filter_zero_flow_edges: true,
    };

    let (n_unresolved, unresolved_type) = process_definition.unresolved_cut_content(&model);
    assert!(!process_definition.half_edge_filters(
        &model,
        &double_double_triangle,
        &[],
        n_unresolved,
        &unresolved_type
    ));
}

// fn fs_generate(loop_count: usize) -> Vec<BareGraph> {
//     let model = load_generic_model("sm_dis");

//     let options = dis_cart_prod(&["d", "d~", "g"], loop_count, &model);
//     chain_dis_generate(&options, &model)
// }

pub(crate) fn dis_options_impl(
    init: &[isize],
    final_states: &[Vec<isize>],
    pert: usize,
    loop_count: usize,
    coupling: usize,
) -> ProcessDefinition {
    let mut amp_coupling = HashMap::new();
    amp_coupling.insert("QED".into(), (2, Some(2)));
    amp_coupling.insert("LT".into(), (1, Some(1)));
    let mut xs_coupling = HashMap::new();
    xs_coupling.insert("QCD".into(), (coupling, Some(coupling)));

    let mut xs_pert = HashMap::new();
    xs_pert.insert("QCD".into(), pert);
    ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: init.iter().map(|a| *a as i64).collect(),
        final_pdgs_lists: final_states
            .iter()
            .map(|f| f.iter().map(|a| *a as i64).collect())
            .collect(),
        loop_count_range: (loop_count, loop_count),
        symmetrize_initial_states: true,
        symmetrize_final_states: true,
        symmetrize_left_right_states: false,
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 0,
        amplitude_filters: FeynGenFilters(vec![
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 9000005, -9000005, 12, 14, 16, 2, 4, 6,
                3, 5, 25, 250, 251, 13, 15,
            ]),
            FeynGenFilter::TadpolesFilter(TadpolesFilterOptions {
                veto_tadpoles_attached_to_massive_lines: true,
                veto_tadpoles_attached_to_massless_lines: true,
                veto_only_scaleless_tadpoles: false,
            }),
            FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions {
                veto_snails_attached_to_massive_lines: false,
                veto_snails_attached_to_massless_lines: true,
                veto_only_scaleless_snails: false,
            }),
            FeynGenFilter::CouplingOrders(amp_coupling),
        ]),
        cross_section_filters: FeynGenFilters(vec![
            FeynGenFilter::SewedFilter(SewedFilterOptions {
                filter_tadpoles: true,
            }),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 9000005, -9000005, 12, 14, 16, 2, 4, 6,
                3, 5, 25, 250, 251, 13, 15,
            ]),
            FeynGenFilter::TadpolesFilter(TadpolesFilterOptions {
                veto_tadpoles_attached_to_massive_lines: true,
                veto_tadpoles_attached_to_massless_lines: true,
                veto_only_scaleless_tadpoles: false,
            }),
            FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions {
                veto_snails_attached_to_massive_lines: false,
                veto_snails_attached_to_massless_lines: true,
                veto_only_scaleless_snails: false,
            }),
            FeynGenFilter::PerturbativeOrders(xs_pert),
            FeynGenFilter::CouplingOrders(xs_coupling),
            FeynGenFilter::LoopCountRange((loop_count, loop_count)),
            FeynGenFilter::BlobRange(1..=1),
            FeynGenFilter::SpectatorRange(0..=1),
        ]),
        filter_self_loop: true,
        folder_name: "N/A".into(),
        process_id: 0,
        graph_prefix: "DIS".into(),
        loop_momentum_bases: None,
        numerator_grouping: NumeratorAwareGraphGroupingOption::OnlyDetectZeroes,
        prefactor: GlobalPrefactor::default(),
        selected_graphs: None,
        vetoed_graphs: None,
        filter_zero_flow_edges: true,
    }
}

#[allow(dead_code)]
pub(crate) fn dis_options(
    init: &[&'static str],
    final_states: &[Vec<&'static str>],
    pert: usize,
    model: &Model,
) -> ProcessDefinition {
    let initial_pdgs: Vec<_> = init
        .iter()
        .map(|a| model.get_particle(a).0.pdg_code)
        .collect();

    let final_pdgs_lists: Vec<_> = final_states
        .iter()
        .map(|f| f.iter().map(|a| model.get_particle(a).0.pdg_code).collect())
        .collect();

    let initial_state_mult = init.len();
    let loop_count = pert + 3 - initial_state_mult;
    let coupling = 2 * pert;
    dis_options_impl(&initial_pdgs, &final_pdgs_lists, pert, loop_count, coupling)
}
struct CombinationsWithRepetition<T> {
    elements: Vec<T>,
    k: usize,
    current: Vec<usize>,
    done: bool,
}

impl<T: Clone + Ord> CombinationsWithRepetition<T> {
    fn new(set: HashSet<T>, k: usize) -> Self {
        let elements: Vec<T> = set.into_iter().sorted().collect();
        CombinationsWithRepetition {
            elements,
            k,
            current: vec![0; k],
            done: false,
        }
    }
}

impl<T: Clone> Iterator for CombinationsWithRepetition<T> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let result = self
            .current
            .iter()
            .map(|&i| self.elements[i].clone())
            .collect();

        // Generate next combination (non-decreasing sequence)
        let n = self.elements.len();
        let mut i = self.k - 1;
        while i < self.k && self.current[i] == n - 1 {
            i = i.wrapping_sub(1);
        }
        if i >= self.k {
            self.done = true;
        } else {
            self.current[i] += 1;
            for j in (i + 1)..self.k {
                self.current[j] = self.current[i];
            }
        }

        Some(result)
    }
}

pub(crate) fn dis_cart_prod_impl(
    initial_states: &[isize],
    loop_count: usize,
) -> Vec<ProcessDefinition> {
    let mut options = vec![];

    let initial_states: HashSet<isize> = initial_states.iter().cloned().collect();

    let initial_state_template = vec![11];
    let final_states = initial_states
        .iter()
        .map(|a| {
            let mut temp = initial_state_template.iter().copied().collect_vec();

            temp.extend([*a]);
            temp
        })
        .collect_vec();

    for initial_state_mult in 1..(loop_count + 2) {
        for mut init_states in
            CombinationsWithRepetition::new(initial_states.clone(), initial_state_mult)
        {
            init_states.extend(initial_state_template.clone());

            let init_states = init_states.into_iter().collect_vec();

            // info!("initial states: {:?}\nfinal states: {:?}\ncross_section_orders:{}\nloop_count:{}\nn_unresolved:{}", init_states, final_states,2*loop_count ,loop_count+ 2 - initial_state_mult,loop_count);

            let pert = loop_count;
            let loop_count = pert + 2 - initial_state_mult;
            let coupling = 2 * pert;

            options.push(dis_options_impl(
                &init_states,
                &final_states,
                pert,
                loop_count,
                coupling,
            ));
        }
    }

    options
}

pub(crate) fn dis_cart_prod(
    initial_states: &[&'static str],
    loop_count: usize,
    model: &Model,
) -> Vec<ProcessDefinition> {
    let initial_states: Vec<_> = initial_states
        .iter()
        .map(|a| model.get_particle(a).0.pdg_code)
        .collect();
    dis_cart_prod_impl(&initial_states, loop_count)
}

pub(crate) fn chain_dis_generate(options: &[ProcessDefinition], model: &Model) -> Vec<Graph> {
    options
        .iter()
        .flat_map(|a| {
            a.generate(
                model,
                &GlobalSettings {
                    n_cores: Parallelisation {
                        feyngen: 10,
                        ..Default::default()
                    },
                    ..Default::default()
                },
            )
            .unwrap()
        })
        .collect()
}

#[test]
fn nlo_fs_dis() {
    // setup_logger().unwrap();

    let model = load_generic_model("sm_dis");
    let options = dis_cart_prod(&["d", "d~", "g"], 1, &model);

    for option in options {
        let diagrams = chain_dis_generate(&[option.clone()], &model);

        let process_name = option
            .initial_pdgs
            .iter()
            .map(|a| model.get_particle_from_pdg(*a as isize).0.name.clone())
            .join(",");

        assert_snapshot!(
            format!("number of diagrams for {}", process_name),
            diagrams.len()
        );
        // let dots = diagrams
        //     .iter()
        //     .map(|a| a.dot())
        //     .collect::<Vec<_>>()
        //     .join("\n");
        // assert_snapshot!(format!("dots for {}", process_name), dots);
    }

    // println!("Number of fs: {}", fs_diagrams.len());
}
