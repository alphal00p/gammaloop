use std::sync::{Arc, Mutex};

use linnet::{
    half_edge::involution::{Flow, Orientation},
    parser::DotEdgeData,
};
use spenso::{
    algebra::complex::Complex,
    structure::{IndexLess, representation::LibraryRep},
};
use symbolica::{domains::float::Complex as SymComplex, prelude::*};

use crate::{
    integrands::process::ParamBuilder,
    model::{
        Model, ModelGammaLoopExt, ParticleGammaLoopExt, ParticleId, ParticleIdGammaLoopExt,
        UFOSymbol,
    },
    momentum::Helicity,
    numerator::aind::Aind,
    utils::{F, FloatLike},
    uv::uv_graph::UVE,
};

use super::parse::ToQuoted;
use crate::graph::Autogen;
use color_eyre::Result;
use eyre::eyre;

/// The model particle of an edge, optionally with a DOT `mass` override.
#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub enum PossibleParticle {
    Particle(ParticleId),
    MassOverriddenParticle { particle: ParticleId, mass: Atom },
}

impl From<ParticleId> for PossibleParticle {
    fn from(particle: ParticleId) -> Self {
        PossibleParticle::Particle(particle)
    }
}

impl PossibleParticle {
    pub fn reverse(&self, model: &Model) -> Self {
        match self {
            PossibleParticle::Particle(particle) => {
                PossibleParticle::Particle(particle.antiparticle(model))
            }
            PossibleParticle::MassOverriddenParticle { particle, mass } => {
                PossibleParticle::MassOverriddenParticle {
                    particle: particle.antiparticle(model),
                    mass: mass.clone(),
                }
            }
        }
    }

    pub fn orientation(&self, model: &Model) -> Orientation {
        let particle = self.particle().resolve(model);
        if !particle.is_fermion() {
            Orientation::Undirected
        } else if particle.pdg_code < 0 {
            Orientation::Reversed
        } else {
            Orientation::Default
        }
    }

    pub fn is_fermion(&self, model: &Model) -> bool {
        self.particle().resolve(model).is_fermion()
    }

    pub fn is_anticommutating(&self, model: &Model) -> bool {
        self.particle().resolve(model).is_anticommutating()
    }

    pub fn is_self_antiparticle(&self, model: &Model) -> bool {
        self.particle().is_self_antiparticle(model)
    }

    pub fn mass_atom(&self, model: &Model) -> Atom {
        match &self {
            PossibleParticle::Particle(particle) => particle.resolve(model).symbolic_mass(model),
            PossibleParticle::MassOverriddenParticle { mass, .. } => mass.clone(),
        }
    }

    pub(crate) fn override_mass(self, mass: Atom) -> Self {
        PossibleParticle::MassOverriddenParticle {
            particle: self.particle(),
            mass,
        }
    }

    pub(crate) fn color_reps(&self, model: &Model, flow: Flow) -> IndexLess {
        self.particle().resolve(model).color_reps(flow)
    }

    pub(crate) fn spin_reps(&self, model: &Model) -> IndexLess<LibraryRep, Aind> {
        self.particle().resolve(model).spin_reps()
    }

    pub(crate) fn particle(&self) -> ParticleId {
        match self {
            PossibleParticle::Particle(particle)
            | PossibleParticle::MassOverriddenParticle { particle, .. } => *particle,
        }
    }

    pub(crate) fn is_massless(&self, model: &Model) -> bool {
        match self {
            PossibleParticle::Particle(particle) => particle.is_massless(model),
            PossibleParticle::MassOverriddenParticle { mass, .. } => mass.is_zero(),
        }
    }

    pub(crate) fn is_massive(&self, model: &Model) -> bool {
        !self.is_massless(model)
    }
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
#[derive(Default)]
pub struct EdgeExtraData {
    /// This controls the power of the momentum in the momtrop sampler. This does *not* affect other aspects of the treatment of the graph by gammaloop.
    pub momtrop_edge_power: Option<Atom>,
    /// This controls the power of the momentum in the vakint evaluation of the graph. This does *not* affect other aspects of the treatment of the graph by gammaloop.
    pub vakint_edge_power: Option<isize>,
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Edge {
    // #[bincode(with_serde)]
    pub name: Autogen<String>,
    pub extra_data: EdgeExtraData,
    // pub edge_type: EdgeType,
    // pub propagator: ArcPropagator,
    pub particle: PossibleParticle,
    pub mass: EdgeMass,
    pub num: Autogen<Atom>,
    // pub spin_num: Atom,
    pub dod: Autogen<i32>,
    pub is_dummy: bool, // #[bincode(with_serde)]
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub enum EdgeMass {
    Zero,
    Value(Complex<F<f64>>),
    ModelVar(Symbol),
    Evaluator(Arc<Mutex<ExpressionEvaluator<Complex<F<f64>>>>>),
}

impl PartialEq for EdgeMass {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (EdgeMass::Zero, EdgeMass::Zero) => true,
            (EdgeMass::Value(a), EdgeMass::Value(b)) => a == b,
            (EdgeMass::ModelVar(a), EdgeMass::ModelVar(b)) => a == b,
            _ => false,
        }
    }
}

impl EdgeMass {
    pub fn from_atom(atom: Atom, model: &Model, paramb: &ParamBuilder) -> Result<Self> {
        if atom.is_zero() {
            return Ok(EdgeMass::Zero);
        } else if let AtomView::Var(v) = atom.as_view() {
            if model.contains_symbol(&UFOSymbol(v.get_symbol())) {
                return Ok(EdgeMass::ModelVar(v.get_symbol()));
            }
        } else if let Ok(a) = SymComplex::<Float>::try_from(&atom) {
            return Ok(EdgeMass::Value(Complex {
                re: F(a.re.to_f64()),
                im: F(a.im.to_f64()),
            }));
        }

        let params: Vec<Atom> = (&paramb.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();

        let a = atom
            .evaluator(&params)
            .function_map(paramb.fn_map.clone())
            .optimization_settings(OptimizationSettings::default())
            .build()
            .map_err(|a| eyre!(a))?;

        Ok(EdgeMass::Evaluator(Arc::new(Mutex::new(a.map_coeff(
            &|r| Complex::new(F::from(&r.re), F::from(&r.im)),
        )))))
    }

    pub fn value<T: FloatLike>(
        &self,
        model: &Model,
        paramb: &ParamBuilder,
    ) -> Option<Complex<F<T>>> {
        match self {
            EdgeMass::Zero => None,
            EdgeMass::Value(v) => Some(*v),
            EdgeMass::ModelVar(s) => model.get_symbol_value(UFOSymbol(*s)),
            EdgeMass::Evaluator(a) => Some(a.lock().unwrap().evaluate_single(&paramb.values[0])),
        }
        .map(|a| a.map_ref(|a| F::from_ff64(*a)))
    }
}

impl UVE for Edge {
    fn mass_atom(&self, model: &Model) -> Atom {
        self.particle
            .mass_atom(model)
            .replace(UFOSymbol::zero().0)
            .with(Atom::Zero)
    }

    fn particle_pdg_code(&self, model: &Model) -> Option<isize> {
        Some(
            self.particle()
                .resolve(model)
                .pdg_code
                .try_into()
                .expect("validated particle PDG codes must fit in isize"),
        )
    }

    fn is_massive(&self, model: &Model) -> bool {
        self.particle.is_massive(model)
    }
}

impl Edge {
    pub fn random_helicity(&self, model: &Model, seed: u64) -> Helicity {
        if let PossibleParticle::Particle(particle) = self.particle {
            particle.resolve(model).random_helicity(seed)
        } else {
            Helicity::ZERO
        }
    }

    pub(crate) fn particle(&self) -> ParticleId {
        self.particle.particle()
    }

    pub(crate) fn mass_value<T: FloatLike>(
        &self,
        model: &Model,
        paramb: &ParamBuilder,
    ) -> Option<Complex<F<T>>> {
        self.mass.value(model, paramb)
    }
}

impl ParseEdge {
    /// Compact DOT attributes, naming the particle so that FeynKit can import them again.
    fn to_dot_data(&self, model: &Model) -> DotEdgeData {
        let mut e = DotEdgeData::empty();
        if let Some(name) = &self.name {
            e.add_statement("name", name.clone());
        }
        match &self.particle {
            PossibleParticle::Particle(particle) => {
                e.add_statement("particle", &particle.resolve(model).name);
            }
            PossibleParticle::MassOverriddenParticle { mass, particle, .. } => {
                e.add_statement("mass", mass.to_quoted());
                e.add_statement("particle", &particle.resolve(model).name);
            }
        }
        if let Some(dod) = &self.dod {
            e.add_statement("dod", *dod);
        }
        if let Some(num) = &self.num {
            e.add_statement("num", num.to_quoted());
        }

        if let Some(mep) = self.momtrop_edge_power.as_ref() {
            e.add_statement("momtrop_edge_power", mep.to_canonical_string());
        }

        if let Some(vak) = self.vakint_edge_power {
            e.add_statement("vakint_edge_power", vak);
        }

        if self.is_dummy {
            e.add_statement("is_dummy", self.is_dummy);
        }
        e
    }
}

impl Edge {
    pub(crate) fn to_parse_edge(&self, include_autogenerated_fields: bool) -> ParseEdge {
        ParseEdge {
            name: self
                .name
                .clone()
                .option_with_generated(include_autogenerated_fields),

            particle: self.particle.clone(),
            dod: self.dod.option_with_generated(include_autogenerated_fields),
            is_dummy: self.is_dummy,
            num: self
                .num
                .clone()
                .option_with_generated(include_autogenerated_fields),
            momtrop_edge_power: self.extra_data.momtrop_edge_power.clone(),
            vakint_edge_power: self.extra_data.vakint_edge_power,
        }
    }

    pub(crate) fn to_dot_data(
        &self,
        model: &Model,
        include_autogenerated_fields: bool,
    ) -> DotEdgeData {
        let parse_edge = self.to_parse_edge(include_autogenerated_fields);
        let mut e = parse_edge.to_dot_data(model);
        if include_autogenerated_fields {
            if self.name.autogenerated {
                e.add_statement("name_autogen", true);
            }
            if self.dod.autogenerated {
                e.add_statement("dod_autogen", true);
            }
            if self.num.autogenerated {
                e.add_statement("num_autogen", true);
            }
        }
        e
    }
}

#[derive(Debug, Clone)]
pub struct ParseEdge {
    pub name: Option<String>,

    pub particle: PossibleParticle,

    pub dod: Option<i32>,
    pub is_dummy: bool,

    /// User provided numerator
    pub num: Option<Atom>,
    pub momtrop_edge_power: Option<Atom>,
    pub vakint_edge_power: Option<isize>,
}

impl UVE for ParseEdge {
    fn mass_atom(&self, model: &Model) -> Atom {
        self.particle.mass_atom(model)
    }

    fn particle_pdg_code(&self, model: &Model) -> Option<isize> {
        Some(
            self.particle
                .particle()
                .resolve(model)
                .pdg_code
                .try_into()
                .expect("validated particle PDG codes must fit in isize"),
        )
    }

    fn is_massive(&self, model: &Model) -> bool {
        self.particle.is_massive(model)
    }
}

impl ParseEdge {
    pub fn new(particle: impl Into<PossibleParticle>) -> Self {
        ParseEdge {
            is_dummy: false,
            name: None,

            particle: particle.into(),
            dod: None,
            num: None,
            momtrop_edge_power: None,
            vakint_edge_power: None,
        }
    }

    pub fn with_label(mut self, label: String) -> Self {
        self.name = Some(label);
        self
    }

    pub fn with_num(mut self, num: Atom) -> Self {
        self.num = Some(num);
        self
    }
}
