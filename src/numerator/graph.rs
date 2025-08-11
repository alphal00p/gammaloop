use linnet::half_edge::{
    involution::{EdgeData, Flow, Orientation},
    subgraph::SubGraph,
};
use spenso::{
    iterators::IteratableTensor,
    network::parsing::ShadowedStructure,
    shadowing::Concretize,
    structure::{PermutedStructure, TensorStructure},
    tensors::parametric::ParamTensor,
};
use symbolica::atom::{Atom, AtomOrView, FunctionBuilder, Symbol};

use crate::{
    new_graph::{Edge, Graph, NumHedgeData},
    utils::GS,
};

use super::aind::Aind;

impl Graph {
    /// Returns the polarizations of the given subgraph. One polarization per half-edge.
    /// If you only want those that are dangling first get the crown of the subgraph.
    pub(crate) fn polarizations_of<S: SubGraph>(&self, subgraph: &S) -> Atom {
        let mut pols = Atom::num(1);

        for h in subgraph.included_iter() {
            let eid = self.underlying[&h];
            pols *= self.underlying.get_edge_data_full(h).polarization(
                &[Atom::num(eid.0 as i32)],
                &self.underlying[h],
                self.underlying.flow(h),
            );
        }

        pols
    }

    pub(crate) fn polarization_parameters_of<S: SubGraph>(&self, subgraph: &S) -> Vec<Atom> {
        let mut pols = Vec::new();

        for h in subgraph.included_iter() {
            let eid = self.underlying[&h];

            let Some(pol) = self
                .underlying
                .get_edge_data_full(h)
                .polarization_structure(
                    &[Atom::num(eid.0 as i32)],
                    &self.underlying[h],
                    self.underlying.flow(h),
                )
            else {
                continue;
            };

            let concrete: ParamTensor<_> = pol
                .structure
                .to_shell()
                .concretize(Some(pol.index_permutation));

            for (_, i) in concrete.iter_flat() {
                pols.push(i.to_owned());
            }
        }

        pols
    }

    pub(crate) fn polarizations(&self) -> Atom {
        self.polarizations_of(&self.underlying.external_filter())
    }
    pub(crate) fn polarization_params(&self) -> Vec<Atom> {
        self.polarization_parameters_of(&self.underlying.external_filter())
    }

    // pub(crate) fn polarizations_values(&self) -> Vec<Atom> {
    //     self.polarizations_of(&self.underlying.external_filter())
    // }
}

pub trait ReversibleEdge {
    fn pdg(&self) -> isize;

    fn pol_symbol(&self, flow: Flow) -> Option<Symbol>;

    fn polarization<'a, T>(&self, add_args: &'a [T], hedge_data: &NumHedgeData, flow: Flow) -> Atom
    where
        &'a T: Into<AtomOrView<'a>>;

    fn polarization_structure<'a, T>(
        &self,
        add_args: &'a [T],
        hedge_data: &NumHedgeData,
        flow: Flow,
    ) -> Option<PermutedStructure<ShadowedStructure<Aind>>>
    where
        &'a T: Into<AtomOrView<'a>>;
}

impl ReversibleEdge for EdgeData<&Edge> {
    fn pdg(&self) -> isize {
        if let Some(p) = self.data.particle() {
            match self.orientation {
                Orientation::Default | Orientation::Undirected => p.pdg_code,
                Orientation::Reversed => -p.pdg_code,
            }
        } else {
            0
        }
    }

    fn pol_symbol(&self, flow: Flow) -> Option<Symbol> {
        if let Some(p) = self.data.particle() {
            match (p.spin, flow) {
                (2, Flow::Sink) => {
                    if self.pdg() > 0 {
                        Some(GS.u)
                    } else {
                        Some(GS.vbar)
                    }
                }
                (2, Flow::Source) => {
                    if self.pdg() > 0 {
                        Some(GS.ubar)
                    } else {
                        Some(GS.v)
                    }
                }
                (3, Flow::Sink) => Some(GS.epsilon),
                (3, Flow::Source) => Some(GS.epsilonbar),
                _ => None,
            }
        } else {
            None
        }
    }

    fn polarization<'a, T>(&self, add_args: &'a [T], hedge_data: &NumHedgeData, flow: Flow) -> Atom
    where
        &'a T: Into<AtomOrView<'a>>,
    {
        if let Some(name) = self.pol_symbol(flow) {
            hedge_data.polarization(FunctionBuilder::new(name).add_args(add_args))
        } else {
            Atom::num(1)
        }
    }

    fn polarization_structure<'a, T>(
        &self,
        add_args: &'a [T],
        hedge_data: &NumHedgeData,
        flow: Flow,
    ) -> Option<PermutedStructure<ShadowedStructure<Aind>>>
    where
        &'a T: Into<AtomOrView<'a>>,
    {
        self.pol_symbol(flow).map(|name| {
            ShadowedStructure::from_iter(
                hedge_data.edge_spin_slots(),
                name,
                Some(add_args.iter().map(|arg| arg.into().into_owned()).collect()),
            )
        })
    }
}

// #[cfg(test)]
// mod test {

//     use env_logger::WriteStyle;
//     use log::{debug, LevelFilter};
//     use spenso::structure::HasStructure;
//     use symbolica::atom::AtomCore;

//     use crate::{
//         dot,
//         new_cs::Amplitude,
//         new_graph::{parse::IntoGraph, FeynmanGraph, Graph},
//         numerator::UnInit,
//         uv::UltravioletGraph,
//         KinematicsSettings, ProcessSettings, Settings,
//     };

//     #[test]
//     fn two_photons() {
//         let _ = env_logger::Builder::new()
//             .filter(None, LevelFilter::Debug)
//             .write_style(WriteStyle::Always)
//             .is_test(true)
//             .try_init();
//         // let _ = env_logger::builder().is_test(true).try_init();

//         let model = crate::tests_from_pytest::load_generic_model("sm");

//         let graph: Graph = dot!(
//         digraph physical_1L_6photons_0 {
//             edge[dod=-1000]
//             node[dod=-1000]
//             ext    [style=invis]
//             ext -> v1 [particle=a id=1];
//             v2 -> ext [particle=a id=0];
//             v1 -> v2 [pdg=1];
//             v2 -> v1 [pdg=1];
//         })
//         .unwrap();
//         println!("{}", graph.dot_serialize());
//         let reps = graph.get_ose_replacements();
//         for r in reps {
//             println!("{r}")
//         }
//         // return;

//         let mut amp: Amplitude = Amplitude::new("name".into());

//         let mut settings = Settings::default();

//         settings.kinematics = KinematicsSettings::random(&graph, 42);
//         amp.add_graph(graph).unwrap();
//         // Amplitude::new(name)

//         let proc_set = ProcessSettings::default();

//         amp.preprocess(&model, &proc_set);
//         let integrand = amp.generate_integrand(settings, &model);

//         // println!("{}", a.factor());

//         //     let mut amp: Amplitude<UnInit> = Amplitude::new("name".into());
//         //     let mut settings = Settings::default();
//         //     for g in graphs {
//         //         settings.kinematics = KinematicsSettings::random(&g, 42);
//         //         amp.add_graph(g).unwrap();
//         //         // Amplitude::new(name)
//         //     }

//         //     let proc_set = ProcessSettings::default();

//         //     amp.preprocess(&model, &proc_set).unwrap();

//         //     let integrand = amp.generate_integrand(settings, &model);
//         // }
//     }

//     #[test]
//     fn six_photons() {
//         let _ = env_logger::Builder::new()
//             .filter(None, LevelFilter::Debug)
//             .write_style(WriteStyle::Always)
//             .is_test(true)
//             .try_init();
//         // let _ = env_logger::builder().is_test(true).try_init();

//         let model = crate::tests_from_pytest::load_generic_model("sm");

//         let graph: Graph = dot!(
//         digraph physical_1L_6photons_0 {
//             ext    [style=invis]
//             ext -> v7 [particle=a];
//             ext -> v8 [particle=a];
//             v9 -> ext [particle=a];
//             v10 -> ext [particle=a];
//             v11 -> ext [particle=a];
//             v12 -> ext [particle=a];
//             v7 -> v8 [pdg=6];
//             v8 -> v9 [pdg=6];
//             v9 -> v10 [pdg=6];
//             v10 -> v11 [pdg=6];
//             v11 -> v12 [pdg=6];
//             v12 -> v7 [pdg=6];
//         })
//         .unwrap();
//         println!("{}", graph.dot_serialize());
//         let reps = graph.get_ose_replacements();
//         for r in reps {
//             println!("{r}")
//         }
//         // return;

//         let mut amp: Amplitude = Amplitude::new("name".into());

//         let mut settings = Settings::default();

//         settings.kinematics = KinematicsSettings::random(&graph, 42);
//         amp.add_graph(graph).unwrap();
//         // Amplitude::new(name)

//         let proc_set = ProcessSettings::default();

//         // amp.graphs[0].generate_cff().unwrap();
//         // let a = amp.graphs[0].build_all_orientations_integrand_atom();
//         // println!("{:>}", a);
//         // debug!("Generating Cff");
//         // amp.graphs[0].generate_cff().unwrap();
//         // debug!("Building Evaluator");
//         // amp.graphs[0].build_evaluator(&model);

//         amp.preprocess(&model, &proc_set);
//         let integrand = amp.generate_integrand(settings, &model);

//         // println!("{}", a.factor());

//         //     let mut amp: Amplitude<UnInit> = Amplitude::new("name".into());
//         //     let mut settings = Settings::default();
//         //     for g in graphs {
//         //         settings.kinematics = KinematicsSettings::random(&g, 42);
//         //         amp.add_graph(g).unwrap();
//         //         // Amplitude::new(name)
//         //     }

//         //     let proc_set = ProcessSettings::default();

//         //     amp.preprocess(&model, &proc_set).unwrap();

//         //     let integrand = amp.generate_integrand(settings, &model);
//         // }
//     }

//     #[test]
//     fn evaluate_pols() {
//         let model = crate::tests_from_pytest::load_generic_model("sm");

//         let graphs: Vec<Graph> = dot!(
//         digraph bxatobx{
//             graph [
//                 // polarizations="1"
//                 color_num="spenso::g(spenso::dind(spenso::cof(3,hedge(0))),spenso::cof(3,hedge(2)))"
//             ]
//             bla    [style=invis]
//             bla -> A:1   [particle=a id=0]
//             bla -> A:2    [particle="b~" id=2]
//             A:0  -> bla  [particle="b~" id=1]
//         })
//         .unwrap();

//         let mut amp: Amplitude = Amplitude::new("name".into());
//         let mut settings = Settings::default();
//         for g in graphs {
//             println!("{}", g.dot_serialize());
//             settings.kinematics = KinematicsSettings::random(&g, 42);
//             amp.add_graph(g).unwrap();
//             // Amplitude::new(name)
//         }

//         let proc_set = ProcessSettings::default();

//         amp.preprocess(&model, &proc_set).unwrap();

//         let integrand = amp.generate_integrand(settings, &model);

//         // integrand.evaluate_sample(sample, wgt, iter, use_f128, max_eval);
//     }

//     #[test]
//     fn pols() {
//         let mut graphs:Vec<Graph> = dot!(
//             digraph bxatobx{
//                 graph [
//                     num="v(0,spenso::bis(4,hedge(0)))*vbar(2,spenso::bis(4,hedge(2)))*ϵ(1,spenso::mink(4,hedge(1)))"
//                 ]
//                 ext    [style=invis]
//                 ext -> A:1   [particle=a id=1]
//                 ext -> A:2    [particle="b~" id=2]
//                 A:0  -> ext  [particle="b~" id=0]
//             }

//             digraph bxatobx2{
//                 graph [
//                     num="v(0,spenso::bis(4,hedge(0)))*vbar(2,spenso::bis(4,hedge(2)))*ϵ(1,spenso::mink(4,hedge(1)))"
//                 ]
//                 ext    [style=invis]
//                 ext -> A:1   [particle=a id=1]
//                 ext -> A:2    [particle="b" id=2 dir=back]
//                 A:0  -> ext  [particle="b" id=0 dir=back]
//             }

//             digraph batob{
//                 graph [
//                     num="ubar(0,spenso::bis(4,hedge(0)))*u(2,spenso::bis(4,hedge(2)))*ϵ(1,spenso::mink(4,hedge(1)))"
//                 ]
//                 ext    [style=invis]
//                 ext -> A:1   [particle=a id=1]
//                 ext -> A:2    [particle="b" id=2]
//                 A:0  -> ext  [particle="b" id=0]
//             }

//             digraph bbato{
//                 graph [
//                     num="vbar(0,spenso::bis(4,hedge(0)))*u(2,spenso::bis(4,hedge(2)))*ϵ(1,spenso::mink(4,hedge(1)))"
//                 ]
//                 ext    [style=invis]
//                 ext -> A:1   [particle=a id=1]
//                 ext -> A:2    [particle="b" id=2]
//                 ext -> A:0  [pdg=-5 id=0]
//             }
//         )
//         .unwrap();

//         for g in &mut graphs {
//             let pols = g.polarizations();
//             // println!("{pols}");
//             let expected = &g.global_prefactor.num;
//             // println!("{expected}");
//             assert_eq!(
//                 expected, &pols,
//                 "\nExpected: {:}\nActual: {:} for graph {:}",
//                 expected, pols, g.name,
//             );
//         }
//     }
// }
