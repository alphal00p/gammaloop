use linnet::half_edge::{
    HedgeGraph,
    involution::{EdgeData, Flow, Orientation},
    subgraph::{SuBitGraph, SubSetLike},
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
    graph::{Edge, Graph, NumHedgeData, edge::ParseEdge},
    model::ArcParticle,
    utils::GS,
};

use super::aind::Aind;

pub trait GeneratePolarizations {
    /// Returns the polarizations of the given subgraph. One polarization per half-edge.
    /// If you only want those that are dangling first get the crown of the subgraph.
    fn generate_polarizations_of<S: SubSetLike>(&self, subgraph: &S) -> Atom;
    fn generate_polarization_parameters_of<S: SubSetLike>(&self, subgraph: &S) -> Vec<Atom>;
    // fn generate_polarizations(&self) -> Atom;
    // fn generate_polarization_params(&self) -> Vec<Atom>;
}

impl Graph {
    pub fn generate_polarizations(&self) -> Atom {
        self.generate_polarizations_of(&self.underlying.external_filter::<SuBitGraph>())
    }
    pub fn generate_polarization_params(&self) -> Vec<Atom> {
        self.generate_polarization_parameters_of(&self.underlying.external_filter::<SuBitGraph>())
    }
}

impl GeneratePolarizations for Graph {
    fn generate_polarizations_of<S: SubSetLike>(&self, subgraph: &S) -> Atom {
        self.underlying.generate_polarizations_of(subgraph)
    }

    fn generate_polarization_parameters_of<S: SubSetLike>(&self, subgraph: &S) -> Vec<Atom> {
        self.underlying
            .generate_polarization_parameters_of(subgraph)
    }

    // fn generate_polarizations(&self) -> Atom {
    //     self.generate_polarizations_of(&self.underlying.external_filter())
    // }
    // fn generate_polarization_params(&self) -> Vec<Atom> {
    //     self.generate_polarization_parameters_of(&self.underlying.external_filter())
    // }

    // pub(crate) fn polarizations_values(&self) -> Vec<Atom> {
    //     self.polarizations_of(&self.underlying.external_filter())
    // }
}

impl<V, E> GeneratePolarizations for HedgeGraph<E, V, NumHedgeData>
where
    for<'a> EdgeData<&'a E>: ReversibleEdge,
{
    fn generate_polarizations_of<S: SubSetLike>(&self, subgraph: &S) -> Atom {
        let mut pols = Atom::num(1);

        for h in subgraph.included_iter() {
            let eid = self[&h];
            pols *= self.get_edge_data_full(h).polarization(
                &[Atom::num(eid.0)],
                &self[h],
                self.flow(h),
            );
        }

        pols
    }

    fn generate_polarization_parameters_of<S: SubSetLike>(&self, subgraph: &S) -> Vec<Atom> {
        let mut pols = Vec::new();

        for h in subgraph.included_iter() {
            let eid = self[&h];

            let Some(pol) = self.get_edge_data_full(h).polarization_structure(
                &[Atom::num(eid.0)],
                &self[h],
                self.flow(h),
            ) else {
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

    // fn generate_polarizations(&self) -> Atom {
    //     self.generate_polarizations_of(&self.external_filter())
    // }
    // fn generate_polarization_params(&self) -> Vec<Atom> {
    //     self.generate_polarization_parameters_of(&self.external_filter())
    // }

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
            // println!("{}{:?}", p.spin, flow);
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

impl ReversibleEdge for EdgeData<&ParseEdge> {
    fn pdg(&self) -> isize {
        if let Some(p) = self.data.particle.particle() {
            match self.orientation {
                Orientation::Default | Orientation::Undirected => p.pdg_code,
                Orientation::Reversed => p.pdg_code,
            }
        } else {
            0
        }
    }

    fn pol_symbol(&self, flow: Flow) -> Option<Symbol> {
        if let Some(p) = self.data.particle.particle() {
            // println!(
            //     "pdg:{},orientation_aware:{},spin:{},flow:{:?}",
            //     p.pdg_code,
            //     self.pdg(),
            //     p.spin,
            //     flow
            // );
            match (p.spin, flow) {
                (2, Flow::Sink) => {
                    if self.pdg() > 0 {
                        // println!("u");
                        Some(GS.u)
                    } else {
                        // println!("vbar");
                        Some(GS.vbar)
                    }
                }
                (2, Flow::Source) => {
                    if self.pdg() > 0 {
                        // println!("ubar");
                        Some(GS.ubar)
                    } else {
                        // println!("v");
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

impl ReversibleEdge for EdgeData<(usize, isize)> {
    fn pdg(&self) -> isize {
        match self.orientation {
            Orientation::Default | Orientation::Undirected => self.data.1,
            Orientation::Reversed => -self.data.1,
        }
    }

    fn pol_symbol(&self, flow: Flow) -> Option<Symbol> {
        match (self.data.0, flow) {
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

impl ReversibleEdge for EdgeData<ArcParticle> {
    fn pdg(&self) -> isize {
        match self.orientation {
            Orientation::Default | Orientation::Undirected => self.data.pdg_code,
            Orientation::Reversed => -self.data.pdg_code,
        }
    }

    fn pol_symbol(&self, flow: Flow) -> Option<Symbol> {
        match (self.data.spin, flow) {
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

#[cfg(test)]
mod test {

    // use env_logger::WriteStyle;
    use idenso::{color::ColorSimplifier, gamma::GammaSimplifier};

    use spenso::network::parsing::{NetworkParse, ParseSettings};
    use symbolica::{
        atom::{Atom, AtomCore},
        parse_lit,
    };

    use crate::{
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::{initialise, test_initialise},
        integrands::process::param_builder::ParamBuilderGraph,
        numerator::aind::Aind,
        processes::{Amplitude, AmplitudeGraph, DotExportSettings},
        settings::{
            GlobalSettings, RuntimeSettings,
            global::GenerationSettings,
            runtime::{LockedRuntimeSettings, kinematic::KinematicsSettings},
        },
        uv::UltravioletGraph,
    };

    #[test]
    fn evaluate_pols() {
        initialise().unwrap();
        let model = crate::utils::test_utils::load_generic_model("sm");

        let graphs: Vec<Graph> = dot!(
        digraph bxatobx{
            graph [
                // polarizations="1"
                color_num="spenso::g(spenso::dind(spenso::cof(3,hedge(0))),spenso::cof(3,hedge(2)))"
            ]
            bla    [style=invis]
            bla -> A:1   [particle=a id=0]
            bla -> A:2    [particle="b~" id=2]
            A:0  -> bla  [particle="b~" id=1]
        })
        .unwrap();

        let mut amp: Amplitude = Amplitude::from_graph_list("name", graphs.clone()).unwrap();
        let mut settings = RuntimeSettings::default();

        for g in graphs {
            println!("{}", g.dot_serialize(&DotExportSettings::default()));
            settings.kinematics = KinematicsSettings::random(&g, 42);

            // Amplitude::new(name)
        }

        let proc_set = GenerationSettings::default();
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();

        let default_runtime_settings = RuntimeSettings::default();
        let locked_runtime_settings = LockedRuntimeSettings::from(&default_runtime_settings);
        amp.preprocess(&model, &proc_set, &locked_runtime_settings, &thread_pool)
            .unwrap();

        amp.build_integrand(
            &model,
            &GlobalSettings::default(),
            (&RuntimeSettings::default()).into(),
            &thread_pool,
        )
        .unwrap();

        // integrand.evaluate_sample(sample, wgt, iter, use_f128, max_eval);
    }

    #[test]
    fn pols() {
        initialise().unwrap();
        let mut graphs: Vec<Graph> = dot!(
            digraph bxatobx{
                graph [
                    num="v(0,spenso::bis(4,hedge(0)))*vbar(2,spenso::bis(4,hedge(2)))*ϵ(1,spenso::mink(4,hedge(1)))"
                ]
                ext    [style=invis]
                ext -> A:1   [particle=a id=1]
                ext -> A:2    [dir=back particle="b~" id=2]
                A:0  -> ext  [dir=back particle="b~" id=0]
            }

            digraph batob{
                graph [
                    num="ubar(0,spenso::bis(4,hedge(0)))*u(2,spenso::bis(4,hedge(2)))*ϵ(1,spenso::mink(4,hedge(1)))"
                ]
                ext    [style=invis]
                ext -> A:1   [particle=a id=1]
                ext -> A:2    [particle="b" id=2]
                A:0  -> ext  [particle="b" id=0]
            }

            digraph bbato{
                graph [
                    num="vbar(0,spenso::bis(4,hedge(0)))*u(2,spenso::bis(4,hedge(2)))*ϵ(1,spenso::mink(4,hedge(1)))"
                ]
                ext    [style=invis]
                ext -> A:1   [particle=a id=1]
                ext -> A:2    [particle="b" id=2]
                ext -> A:0  [dir=back pdg=-5 id=0]
            }

            digraph bbato{
                graph [
                    num="vbar(0,spenso::bis(4,hedge(0)))*u(2,spenso::bis(4,hedge(2)))*ϵbar(1,spenso::mink(4,hedge(1)))"
                ]
                ext    [style=invis]
                A:1 -> ext  [particle=a id=1]
                ext -> A:2    [particle="b" id=2]
                ext -> A:0  [dir=back pdg=-5 id=0]
            }

        )
        .unwrap();

        for g in &mut graphs {
            let pols = g.generate_polarizations();
            println!("{pols}");
            let expected = &g.global_prefactor.num;
            println!("{expected}");
            assert_eq!(
                expected, &pols,
                "\nExpected: {:}\nActual: {:} for graph {:}",
                expected, pols, g.name,
            );
        }
    }

    #[test]
    fn vertex_rule() {
        let mut graphs: Vec<Graph> = dot!(
            digraph dxda{
                ext [style=invis]
                ext->v1:0[particle="d" id=0]
                ext->v1:1[particle="d~" id=1]
                v1:2->ext[particle="a" id=2]
            }

        )
        .unwrap();

        for g in &mut graphs {
            let mut out = String::new();
            g.dot_serialize_fmt(&mut out, &DotExportSettings::default())
                .unwrap();
            println!("{}", out);

            assert!(g.iter_nodes().all(|(_, _, v)| {
                let a = v.vertex_rule.as_ref().unwrap().name.as_str();
                a == "V_71"
            }),);
        }
    }
    #[test]
    fn pslash() {
        let mut graphs: Vec<Graph> = dot!(
            digraph dxda{
                ext [style=invis]
                node[num=1]
                ext->v1:0[particle="d" id=0]
                v1:1->v2:2[particle="d" id=1]
                v2:3->ext[particle="d" id=2]
            }

        )
        .unwrap();

        for g in &mut graphs {
            let mut out = String::new();
            g.dot_serialize_fmt(&mut out, &DotExportSettings::default())
                .unwrap();
            println!("{}", out);
        }
    }

    #[test]
    fn qqx_aaa_tree() {
        test_initialise();

        let mut graph:AmplitudeGraph = dot!(digraph qqx_aaa_tree_1 {
                    num="spenso::g(spenso::dind(spenso::cof(3, hedge(1))), spenso::cof(3, hedge(2)))/3"
                    ext    [style=invis]
                    ext -> v1:1 [particle="d" id=1];
                    ext -> v3:2 [particle="d~" id=2];
                    v1:3 -> ext [particle="a" id=3];
                    v2:4 -> ext [particle="a" id=4];
                    v3:0 -> ext [particle="a" id=0];
                    v1 -> v2 [particle="d" id=5];
                    v2 -> v3 [particle="d" id=6];
        }).unwrap();

        let set = GenerationSettings::default();
        // let vk_settings = set.uv.vakint.true_settings();
        let vk = crate::utils::vakint().unwrap();

        // let model = crate::utils::test_utils::load_generic_model("sm");

        graph.generate_cff().unwrap();
        graph
            .build_integrands(&GenerationSettings::default(), vk)
            .unwrap();

        println!("{}", graph.derived_data.all_mighty_integrand);
    }
    #[test]
    fn tree() {
        initialise().unwrap();
        let mut graphs: Vec<Graph> = dot!(
            digraph qqx_aaa_tree_1 {
                        num="spenso::g(spenso::dind(spenso::cof(3, hedge(1))), spenso::cof(3, hedge(2)))/3"
                        ext    [style=invis]
                        ext -> v1:1 [particle="d" id=1];
                        ext -> v3:2 [dir=back particle="d~" id=2];
                        v1:3 -> ext [particle="a" id=3];
                        v2:4 -> ext [particle="a" id=4];
                        v3:0 -> ext [particle="a" id=0];
                        v1 -> v2 [particle="d" id=5];
                        v2 -> v3 [particle="d" id=6];
            }

            digraph qqx_aaa_tree_1 {
                        ext    [style=invis]
                        ext -> v1:1 [particle="d" id=1];
                        ext -> v3:2 [dir=back particle="d~" id=2];
                        v1:3 -> ext [particle="a" id=3];
                        v2:4 -> ext [particle="a" id=4];
                        v3:0 -> ext [particle="a" id=0];
                        v1 -> v2 [particle="d" id=5];
                        v2 -> v3 [particle="d" id=6];
                        num="spenso::g(spenso::dind(spenso::cof(3, hedge(1))), spenso::cof(3, hedge(2)))/3"
            }


            digraph qqx_aaa_tree_1_glob {
            ext [style=invis];
            ext -> v1:1 [particle="d", id=1];
            ext -> v3:2 [dir=back particle="d~", id=2];
            v1:3 -> ext [particle="a", id=3];
            v2:4 -> ext [particle="a", id=4];
            v3:0 -> ext [particle="a", id=0];
            v1 -> v2 [particle="d", id=5];
            v2 -> v3 [particle="d", id=6];
            num=" UFO::GC_1^3
                *spenso::g(spenso::cof(3,hedge(1)),spenso::dind(spenso::cof(3,hedge(5))))
                *spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(1)),spenso::mink(4,hedge(3)))

                *spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(6))))
                *Q(5,spenso::mink(4,edge(5,1)))
                *spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(5)),spenso::mink(4,edge(5,1)))


                *spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(7))))
                *spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,hedge(4)))

                *spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(8))))
                *spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(7)),spenso::mink(4,edge(6,1)))
                *Q(6,spenso::mink(4,edge(6,1)))

                *spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(2))))
                *spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(8)),spenso::mink(4,hedge(0)))

                   "
             overall_factor="1"
                projector="u(1,spenso::bis(4,hedge(1)))
            *vbar(2,spenso::bis(4,hedge(2)))
            *ϵbar(0,spenso::mink(4,hedge(0)))
            *ϵbar(3,spenso::mink(4,hedge(3)))
            *ϵbar(4,spenso::mink(4,hedge(4)))
            *spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(1))))/3"

            edge [num="1"];
            node [num="1"];
            }


        )
        .unwrap();

        let mut a: Option<Atom> = None;
        for g in &mut graphs {
            let mut out = String::new();

            let new_a = (g
                .numerator(&g.full_filter(), &g.empty_subgraph())
                .state
                .expr
                * &g.global_prefactor.projector
                * &g.global_prefactor.num
                * &g.overall_factor)
                .simplify_color();
            println!("New:{new_a}");
            if let Some(a) = &a {
                println!("Old:{a}");
                assert_eq!(&new_a, a, "{}", (&new_a / a).expand().to_canonical_string());
            } else {
                a = Some(new_a);
            }
            g.dot_serialize_fmt(&mut out, &DotExportSettings::default())
                .unwrap();
            println!("{}", out);
        }

        let mut a = Amplitude::from_graph_list("test", graphs.clone()).unwrap();

        let model = crate::utils::test_utils::load_generic_model("sm");

        let generation_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();

        let default_runtime_settings = RuntimeSettings::default();
        let locked_runtime_settings = LockedRuntimeSettings::from(&default_runtime_settings);
        a.preprocess(
            &model,
            &GenerationSettings::default(),
            &locked_runtime_settings,
            &generation_pool,
        )
        .unwrap();
    }

    #[test]
    fn dod_override() {
        let _gr:Vec<Graph> = dot!(digraph g{
            node [num=1]
            edge [num=1]
            a->b[dod=-100]
            b->c[dod="-100"]
        }
        digraph triangle_ct {
        num="(-1*_gammaloop::P(1,spenso::cind(1))*_gammaloop::P(2,spenso::cind(1))+-1*_gammaloop::P(1,spenso::cind(2))*_gammaloop::P(2,spenso::cind(2))+-1*_gammaloop::P(1,spenso::cind(3))*_gammaloop::P(2,spenso::cind(3))+_gammaloop::P(1,spenso::cind(0))*_gammaloop::P(2,spenso::cind(0)))^-2*(-1*_gammaloop::P(2,spenso::mink(4,python::mu2))+-1*_gammaloop::Q(6,spenso::mink(4,python::mu2))+_gammaloop::P(1,spenso::mink(4,python::mu7)))*(-1*_gammaloop::Q(7,spenso::mink(4,python::mu2))+_gammaloop::P(2,spenso::mink(4,python::mu2)))*-1/4*_gammaloop::G^2*_gammaloop::P(1,spenso::mink(4,python::mu3))*_gammaloop::P(1,spenso::mink(4,python::mu5))*_gammaloop::P(2,spenso::mink(4,python::mu4))*_gammaloop::P(2,spenso::mink(4,python::mu6))*spenso::gamma(spenso::bis(4,_gammaloop::hedge(2)),spenso::bis(4,python::s1),spenso::mink(4,python::mu1))*spenso::gamma(spenso::bis(4,python::s1),spenso::bis(4,python::s2),spenso::mink(4,python::mu2))*spenso::gamma(spenso::bis(4,python::s2),spenso::bis(4,python::s3),spenso::mink(4,python::mu3))*spenso::gamma(spenso::bis(4,python::s3),spenso::bis(4,python::tree_form_factor_spinor_2),spenso::mink(4,python::mu4))*spenso::gamma(spenso::bis(4,python::s4),spenso::bis(4,python::s6),spenso::mink(4,python::mu7))*spenso::gamma(spenso::bis(4,python::s5),spenso::bis(4,python::s4),spenso::mink(4,python::mu6))*spenso::gamma(spenso::bis(4,python::s6),spenso::bis(4,_gammaloop::hedge(1)),spenso::mink(4,python::mu1))*spenso::gamma(spenso::bis(4,python::tree_form_factor_spinor_1),spenso::bis(4,python::s5),spenso::mink(4,python::mu5))";
        overall_factor="1";
        projector="((-1*_gammaloop::P(3,spenso::cind(0))+-1*_gammaloop::P(4,spenso::cind(0))+_gammaloop::P(1,spenso::cind(0)))^2+(-1*_gammaloop::P(3,spenso::cind(1))+-1*_gammaloop::P(4,spenso::cind(1))+_gammaloop::P(1,spenso::cind(1)))^2*-1+(-1*_gammaloop::P(3,spenso::cind(2))+-1*_gammaloop::P(4,spenso::cind(2))+_gammaloop::P(1,spenso::cind(2)))^2*-1+(-1*_gammaloop::P(3,spenso::cind(3))+-1*_gammaloop::P(4,spenso::cind(3))+_gammaloop::P(1,spenso::cind(3)))^2*-1)^-1*((-1*_gammaloop::P(3,spenso::cind(0))+_gammaloop::P(1,spenso::cind(0)))^2+(-1*_gammaloop::P(3,spenso::cind(1))+_gammaloop::P(1,spenso::cind(1)))^2*-1+(-1*_gammaloop::P(3,spenso::cind(2))+_gammaloop::P(1,spenso::cind(2)))^2*-1+(-1*_gammaloop::P(3,spenso::cind(3))+_gammaloop::P(1,spenso::cind(3)))^2*-1)^-1*(-1*_gammaloop::P(3,spenso::mink(4,_gammaloop::edge(5,1)))+_gammaloop::P(1,spenso::mink(4,_gammaloop::edge(5,1))))*(-1*_gammaloop::P(3,spenso::mink(4,_gammaloop::edge(6,1)))+-1*_gammaloop::P(4,spenso::mink(4,_gammaloop::edge(6,1)))+_gammaloop::P(1,spenso::mink(4,_gammaloop::edge(6,1))))*-1/27*_gammaloop::ee^3*_gammaloop::u(1,spenso::bis(4,_gammaloop::hedge(1)))*_gammaloop::vbar(2,spenso::bis(4,_gammaloop::hedge(2)))*_gammaloop::ϵbar(0,spenso::mink(4,_gammaloop::hedge(0)))*_gammaloop::ϵbar(3,spenso::mink(4,_gammaloop::hedge(3)))*_gammaloop::ϵbar(4,spenso::mink(4,_gammaloop::hedge(4)))*spenso::gamma(spenso::bis(4,_gammaloop::hedge(5)),spenso::bis(4,python::tree_form_factor_spinor_1),spenso::mink(4,_gammaloop::hedge(3)))*spenso::gamma(spenso::bis(4,_gammaloop::hedge(6)),spenso::bis(4,_gammaloop::hedge(5)),spenso::mink(4,_gammaloop::edge(5,1)))*spenso::gamma(spenso::bis(4,_gammaloop::hedge(7)),spenso::bis(4,_gammaloop::hedge(6)),spenso::mink(4,_gammaloop::hedge(4)))*spenso::gamma(spenso::bis(4,_gammaloop::hedge(8)),spenso::bis(4,_gammaloop::hedge(7)),spenso::mink(4,_gammaloop::edge(6,1)))*spenso::gamma(spenso::bis(4,python::tree_form_factor_spinor_2),spenso::bis(4,_gammaloop::hedge(8)),spenso::mink(4,_gammaloop::hedge(0)))";
        edge [num="1", dod="-100"];
        node [num="1", dod="-100"];
        ext [style=invis];
        ext -> vl1:1 [particle="d", id=1];
        ext -> vl2:2 [dir=back particle="d~", id=2];
        v1:3 -> ext [id=3];
        v1:4 -> ext [id=4];
        v1:0 -> ext [id=0];
        vl1 -> v1 [particle="d", id=5];
        v1 -> vl2 [particle="d", id=6];
        vl1 -> vl2 [particle="g", id=7, lmb_id=0];
        }
)
        .unwrap();
    }

    #[test]
    fn simplify() {
        test_initialise().unwrap();
        let expr = parse_lit!(
            -GC_1 * GC_11
                ^ 2 * Q(4, mink(dim, gammalooprs::edge(4, 1)))
                    * Q(5, mink(dim, gammalooprs::edge(5, 1)))
                    * g(
                        mink(dim, gammalooprs::hedge(16)),
                        mink(dim, gammalooprs::hedge(17))
                    )
                    * g(
                        coad(8, gammalooprs::hedge(16)),
                        coad(8, gammalooprs::hedge(17))
                    )
                    * g(
                        dind(cof(3, gammalooprs::hedge(4))),
                        cof(3, gammalooprs::hedge(6))
                    )
                    * g(
                        dind(cof(3, gammalooprs::hedge(5))),
                        cof(3, gammalooprs::hedge(4))
                    )
                    * g(
                        dind(cof(3, gammalooprs::hedge(6))),
                        cof(3, gammalooprs::hedge(7))
                    )
                    * gamma(
                        bis(4, gammalooprs::hedge(4)),
                        bis(4, gammalooprs::hedge(6)),
                        mink(dim, gammalooprs::hedge(8))
                    )
                    * gamma(
                        bis(4, gammalooprs::hedge(7)),
                        bis(4, gammalooprs::hedge(13)),
                        mink(dim, gammalooprs::hedge(17))
                    )
                    * gamma(
                        bis(4, gammalooprs::hedge(11)),
                        bis(4, gammalooprs::hedge(5)),
                        mink(dim, gammalooprs::hedge(16))
                    )
                    * gamma(
                        bis(4, gammalooprs::hedge(5)),
                        bis(4, gammalooprs::hedge(4)),
                        mink(dim, gammalooprs::edge(4, 1))
                    )
                    * gamma(
                        bis(4, gammalooprs::hedge(6)),
                        bis(4, gammalooprs::hedge(7)),
                        mink(dim, gammalooprs::edge(5, 1))
                    )
                    * t(
                        coad(8, gammalooprs::hedge(16)),
                        cof(3, gammalooprs::hedge(5)),
                        dind(cof(3, gammalooprs::hedge(11)))
                    )
                    * t(
                        coad(8, gammalooprs::hedge(17)),
                        cof(3, gammalooprs::hedge(13)),
                        dind(cof(3, gammalooprs::hedge(7)))
                    ),
            default_namespace = "spenso"
        );

        let net = expr
            .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
            .unwrap();
        println!("{}", net.dot_pretty());

        let _ = expr.simplify_gamma();
    }
}
