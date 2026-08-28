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
    graph::{Edge, Graph, HedgeData, edge::ParseEdge},
    model::{Model, ParticleId, ParticleIdGammaLoopExt},
    utils::GS,
};

use super::aind::Aind;

pub trait HedgePolarizationData {
    fn polarization(&self, builder: FunctionBuilder) -> Atom;

    fn edge_spin_slots<'a>(
        &'a self,
    ) -> impl Iterator<Item = spenso::structure::representation::LibrarySlot<Aind>> + 'a;
}

impl HedgePolarizationData for HedgeData {
    fn polarization(&self, builder: FunctionBuilder) -> Atom {
        HedgeData::polarization(self, builder)
    }

    fn edge_spin_slots<'a>(
        &'a self,
    ) -> impl Iterator<Item = spenso::structure::representation::LibrarySlot<Aind>> + 'a {
        HedgeData::edge_spin_slots(self)
    }
}

pub trait GeneratePolarizations {
    /// Returns the polarizations of the given subgraph. One polarization per half-edge.
    /// If you only want those that are dangling first get the crown of the subgraph.
    fn generate_polarizations_of<S: SubSetLike>(&self, model: &Model, subgraph: &S) -> Atom;
    fn generate_polarization_parameters_of<S: SubSetLike>(
        &self,
        model: &Model,
        subgraph: &S,
    ) -> Vec<Atom>;
    // fn generate_polarizations(&self) -> Atom;
    // fn generate_polarization_params(&self) -> Vec<Atom>;
}

impl Graph {
    pub fn generate_polarizations(&self, model: &Model) -> Atom {
        self.generate_polarizations_of(model, &self.underlying.external_filter::<SuBitGraph>())
    }
    pub fn generate_polarization_params(&self, model: &Model) -> Vec<Atom> {
        self.generate_polarization_parameters_of(
            model,
            &self.underlying.external_filter::<SuBitGraph>(),
        )
    }
}

impl GeneratePolarizations for Graph {
    fn generate_polarizations_of<S: SubSetLike>(&self, model: &Model, subgraph: &S) -> Atom {
        self.underlying.generate_polarizations_of(model, subgraph)
    }

    fn generate_polarization_parameters_of<S: SubSetLike>(
        &self,
        model: &Model,
        subgraph: &S,
    ) -> Vec<Atom> {
        self.underlying
            .generate_polarization_parameters_of(model, subgraph)
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

impl<V, E> GeneratePolarizations for HedgeGraph<E, V, HedgeData>
where
    for<'a> EdgeData<&'a E>: ReversibleEdge,
{
    fn generate_polarizations_of<S: SubSetLike>(&self, model: &Model, subgraph: &S) -> Atom {
        let mut pols = Atom::num(1);

        for h in subgraph.included_iter() {
            let eid = self[&h];
            pols *= self.get_edge_data_full(h).polarization(
                model,
                &[Atom::num(eid.0)],
                &self[h],
                self.flow(h),
            );
        }

        pols
    }

    fn generate_polarization_parameters_of<S: SubSetLike>(
        &self,
        model: &Model,
        subgraph: &S,
    ) -> Vec<Atom> {
        let mut pols = Vec::new();

        for h in subgraph.included_iter() {
            let eid = self[&h];

            let Some(pol) = self.get_edge_data_full(h).polarization_structure(
                model,
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
    fn pdg(&self, model: &Model) -> isize;

    fn spin(&self, model: &Model) -> Option<i64>;

    fn pol_symbol(&self, model: &Model, flow: Flow) -> Option<Symbol> {
        let pdg = self.pdg(model);
        match (self.spin(model)?, flow) {
            (2, Flow::Sink) => Some(if pdg > 0 { GS.u } else { GS.vbar }),
            (2, Flow::Source) => Some(if pdg > 0 { GS.ubar } else { GS.v }),
            (3, Flow::Sink) => Some(GS.epsilon),
            (3, Flow::Source) => Some(GS.epsilonbar),
            _ => None,
        }
    }

    fn polarization<'a, T, H>(
        &self,
        model: &Model,
        add_args: &'a [T],
        hedge_data: &H,
        flow: Flow,
    ) -> Atom
    where
        &'a T: Into<AtomOrView<'a>>,
        H: HedgePolarizationData,
    {
        if let Some(name) = self.pol_symbol(model, flow) {
            hedge_data.polarization(FunctionBuilder::new(name).add_args(add_args))
        } else {
            Atom::num(1)
        }
    }

    fn polarization_structure<'a, T, H>(
        &self,
        model: &Model,
        add_args: &'a [T],
        hedge_data: &H,
        flow: Flow,
    ) -> Option<PermutedStructure<ShadowedStructure<Aind>>>
    where
        &'a T: Into<AtomOrView<'a>>,
        H: HedgePolarizationData,
    {
        self.pol_symbol(model, flow).map(|name| {
            ShadowedStructure::from_iter(
                hedge_data.edge_spin_slots(),
                name,
                Some(add_args.iter().map(|arg| arg.into().into_owned()).collect()),
            )
        })
    }
}

impl ReversibleEdge for EdgeData<&Edge> {
    fn pdg(&self, model: &Model) -> isize {
        self.data
            .particle()
            .map(|particle| oriented_pdg(particle.resolve(model).pdg_code, self.orientation))
            .unwrap_or(0)
    }

    fn spin(&self, model: &Model) -> Option<i64> {
        self.data
            .particle()
            .map(|particle| particle.resolve(model).spin)
    }
}

impl ReversibleEdge for EdgeData<&ParseEdge> {
    fn pdg(&self, model: &Model) -> isize {
        self.data
            .particle
            .particle()
            .map(|particle| oriented_pdg(particle.resolve(model).pdg_code, self.orientation))
            .unwrap_or(0)
    }

    fn spin(&self, model: &Model) -> Option<i64> {
        self.data
            .particle
            .particle()
            .map(|particle| particle.resolve(model).spin)
    }
}

impl ReversibleEdge for EdgeData<(usize, isize)> {
    fn pdg(&self, _model: &Model) -> isize {
        match self.orientation {
            Orientation::Default | Orientation::Undirected => self.data.1,
            Orientation::Reversed => -self.data.1,
        }
    }

    fn spin(&self, _model: &Model) -> Option<i64> {
        Some(self.data.0 as i64)
    }
}

impl ReversibleEdge for EdgeData<ParticleId> {
    fn pdg(&self, model: &Model) -> isize {
        oriented_pdg(self.data.resolve(model).pdg_code, self.orientation)
    }

    fn spin(&self, model: &Model) -> Option<i64> {
        Some(self.data.resolve(model).spin)
    }
}

fn oriented_pdg(pdg: i64, orientation: Orientation) -> isize {
    let pdg = isize::try_from(pdg).expect("PDG code must fit in an isize");
    match orientation {
        Orientation::Default | Orientation::Undirected => pdg,
        Orientation::Reversed => -pdg,
    }
}

#[cfg(test)]
mod test {

    // use env_logger::WriteStyle;
    use idenso::{dirac::GammaSimplifier, tensor::SymbolicNetParse};

    use spenso::network::parsing::ParseSettings;
    use symbolica::parse_lit;

    use crate::{initialisation::test_initialise, numerator::aind::Aind};

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
