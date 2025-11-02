use std::{cell::RefCell, fmt::Display};

use eyre::Context;
use linnet::{
    half_edge::{
        HedgeGraph,
        involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Orientation},
    },
    parser::{DotEdgeData, DotHedgeData, DotVertexData},
};
use spenso::{
    algebra::complex::Complex,
    structure::{IndexLess, ScalarStructure, representation::LibraryRep},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::float::{Complex as SymComplex, Float},
    evaluate::{ExpressionEvaluator, OptimizationSettings},
    parse,
};

use crate::{
    feyngen::diagram_generator::EdgeColor,
    gammaloop_integrand::ParamBuilder,
    model::{ArcParticle, Model, UFOSymbol},
    momentum::Helicity,
    momentum_sample::LoopIndex,
    numerator::aind::{Aind, NewAind},
    utils::{F, FloatLike, GS},
    uv::uv_graph::UVE,
};

use super::parse::{StripParse, ToQuoted};
use color_eyre::Result;
use eyre::eyre;

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub enum PossibleParticle {
    Particle(ArcParticle),
    MassOverriddenParticle { particle: ArcParticle, mass: Atom },
    JustMass { expr: Atom },
}

impl From<ArcParticle> for PossibleParticle {
    fn from(particle: ArcParticle) -> Self {
        PossibleParticle::Particle(particle)
    }
}

impl From<Atom> for PossibleParticle {
    fn from(atom: Atom) -> Self {
        PossibleParticle::JustMass { expr: atom }
    }
}

impl From<()> for PossibleParticle {
    fn from(_: ()) -> Self {
        PossibleParticle::JustMass { expr: Atom::Zero }
    }
}

impl PossibleParticle {
    pub fn reverse(&self, model: &Model) -> Self {
        match self {
            PossibleParticle::Particle(p) => PossibleParticle::Particle(p.get_anti_particle(model)),
            PossibleParticle::MassOverriddenParticle { particle, mass } => {
                PossibleParticle::MassOverriddenParticle {
                    particle: particle.get_anti_particle(model),
                    mass: mass.clone(),
                }
            }
            PossibleParticle::JustMass { expr } => {
                PossibleParticle::JustMass { expr: expr.clone() }
            }
        }
    }

    pub fn orientation(&self) -> Orientation {
        self.particle()
            .map(|a| {
                if a.is_fermion() {
                    if a.pdg_code < 0 {
                        Orientation::Reversed
                    } else {
                        Orientation::Default
                    }
                } else {
                    Orientation::Undirected
                }
            })
            .unwrap_or(Orientation::Undirected)
    }

    pub fn is_fermion(&self) -> bool {
        match self {
            PossibleParticle::Particle(p) => p.is_fermion(),
            PossibleParticle::MassOverriddenParticle { particle, .. } => particle.is_fermion(),
            PossibleParticle::JustMass { .. } => false,
        }
    }

    pub fn is_self_antiparticle(&self) -> bool {
        match self {
            PossibleParticle::Particle(p) => p.is_self_antiparticle(),
            PossibleParticle::MassOverriddenParticle { particle, .. } => {
                particle.is_self_antiparticle()
            }
            PossibleParticle::JustMass { .. } => false,
        }
    }
    pub fn mass_atom(&self) -> Atom {
        match &self {
            PossibleParticle::JustMass { expr, .. } => expr.clone(),
            PossibleParticle::Particle(p) => p.mass.0.into(),
            PossibleParticle::MassOverriddenParticle { mass, .. } => mass.clone(),
        }
    }
    // pub fn just_mass(mass:Atom,)

    pub fn zero() -> Self {
        ().into()
    }

    pub(crate) fn override_mass(self, mass: Option<Atom>) -> Self {
        let Some(mass) = mass else {
            return self;
        };

        match self {
            PossibleParticle::JustMass { .. } => PossibleParticle::JustMass { expr: mass },
            PossibleParticle::MassOverriddenParticle { particle, .. }
            | PossibleParticle::Particle(particle) => {
                PossibleParticle::MassOverriddenParticle { particle, mass }
            }
        }
    }

    pub(crate) fn color_reps(&self, flow: Flow) -> IndexLess {
        match self {
            PossibleParticle::Particle(p) => p.color_reps(flow),
            _ => IndexLess::scalar_structure(),
        }
    }

    pub(crate) fn spin_reps(&self) -> IndexLess<LibraryRep, Aind> {
        match self {
            PossibleParticle::Particle(p) => p.spin_reps(),
            _ => IndexLess::scalar_structure(),
        }
    }

    pub(crate) fn particle(&self) -> Option<ArcParticle> {
        match self {
            PossibleParticle::Particle(particle)
            | PossibleParticle::MassOverriddenParticle { particle, .. } => Some(particle.clone()),
            _ => None,
        }
    }

    pub(crate) fn is_massless(&self) -> bool {
        match self {
            PossibleParticle::JustMass { expr } => expr.is_zero(),
            PossibleParticle::Particle(p) => p.is_massless(),
            PossibleParticle::MassOverriddenParticle { mass, .. } => mass.is_zero(),
        }
    }

    pub(crate) fn is_massive(&self) -> bool {
        !self.is_massless()
    }
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct EdgeExtraData {
    /// This controls the power of the momentum in the momtrop sampler. This does *not* affect other aspects of the treatment of the graph by gammaloop.
    pub momtrop_edge_power: Option<Atom>,
    /// This controls the power of the momentum in the vakint evaluation of the graph. This does *not* affect other aspects of the treatment of the graph by gammaloop.
    pub vakint_edge_power: Option<isize>,
}

impl Default for EdgeExtraData {
    fn default() -> Self {
        EdgeExtraData {
            momtrop_edge_power: None,
            vakint_edge_power: None,
        }
    }
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Edge {
    // #[bincode(with_serde)]
    pub name: String,
    pub extra_data: EdgeExtraData,
    // pub edge_type: EdgeType,
    // pub propagator: ArcPropagator,
    pub particle: PossibleParticle,
    pub mass: EdgeMass,
    pub num: Atom,
    // pub spin_num: Atom,
    pub dod: i32,
    pub is_dummy: bool, // #[bincode(with_serde)]
}

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub enum EdgeMass {
    Zero,
    Value(Complex<F<f64>>),
    ModelVar(Symbol),
    Evaluator(RefCell<ExpressionEvaluator<Complex<F<f64>>>>),
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
                re: F(a.re.into_inner().to_f64()),
                im: F(a.im.into_inner().to_f64()),
            }));
        }

        let params: Vec<Atom> = (&paramb.pairs)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();

        let a = atom
            .evaluator(&paramb.fn_map, &params, OptimizationSettings::default())
            .map_err(|a| eyre!(a))?;

        Ok(EdgeMass::Evaluator(RefCell::new(a.map_coeff(&|r| {
            Complex::new(F::from(&r.re), F::from(&r.im))
        }))))
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
            EdgeMass::Evaluator(a) => Some(a.borrow_mut().evaluate_single(&paramb.values)),
        }
        .map(|a| a.map_ref(|a| F::from_ff64(*a)))
    }
}

impl UVE for Edge {
    fn mass_atom(&self) -> Atom {
        match &self.particle {
            PossibleParticle::JustMass { expr, .. } => expr.clone(),
            PossibleParticle::Particle(p) => p.mass.0.into(),
            PossibleParticle::MassOverriddenParticle { mass, .. } => mass.clone(),
        }
    }
}

impl Edge {
    pub fn random_helicity(&self, seed: u64) -> Helicity {
        if let PossibleParticle::Particle(p) = &self.particle {
            p.random_helicity(seed)
        } else {
            Helicity::Zero
        }
    }

    pub(crate) fn particle(&self) -> Option<ArcParticle> {
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

impl From<&ParseEdge> for DotEdgeData {
    fn from(value: &ParseEdge) -> Self {
        let mut e = DotEdgeData::empty();
        if let Some(name) = &value.name {
            e.add_statement("name", name.clone());
        }
        match &value.particle {
            PossibleParticle::Particle(p) => {
                e.add_statement("particle", format!("\"{}\"", p.name));
            }
            PossibleParticle::JustMass { expr, .. } => {
                e.add_statement("mass", expr.to_quoted());
            }
            PossibleParticle::MassOverriddenParticle { mass, particle, .. } => {
                e.add_statement("mass", mass.to_quoted());
                e.add_statement("particle", format!("\"{}\"", particle.name));
            }
        }
        if let Some(lmb_id) = &value.lmb_id {
            e.add_statement("lmb_id", usize::from(*lmb_id));
        }
        if let Some(cut) = &value.is_cut {
            e.add_statement("is_cut", usize::from(*cut));
        }
        if let Some(dod) = &value.dod {
            e.add_statement("dod", *dod);
        }
        if let Some(num) = &value.num {
            e.add_statement("num", num.to_quoted());
        }

        if let Some(mep) = value.momtrop_edge_power.as_ref() {
            e.add_statement("momtrop_edge_power", mep.to_canonical_string());
        }

        if let Some(vak) = value.vakint_edge_power {
            e.add_statement("vakint_edge_power", vak);
        }

        // e.add_statement("color_num", value.color_num.to_quoted());
        e.add_statement("is_dummy", value.is_dummy);
        e
    }
}

impl From<&Edge> for DotEdgeData {
    fn from(value: &Edge) -> Self {
        let mut e = DotEdgeData::empty();
        e.add_statement("name", value.name.clone());
        match &value.particle {
            PossibleParticle::Particle(p) => {
                e.add_statement("particle", format!("\"{}\"", p.name));
            }
            PossibleParticle::JustMass { expr, .. } => {
                e.add_statement("mass", expr.to_quoted());
            }
            PossibleParticle::MassOverriddenParticle { mass, particle, .. } => {
                e.add_statement("mass", mass.to_quoted());
                e.add_statement("particle", format!("\"{}\"", particle.name));
            }
        }
        e.add_statement("dod", value.dod);
        e.add_statement("num", value.num.to_quoted());
        // e.add_statement("color_num", value.color_num.to_quoted());
        e.add_statement("is_dummy", value.is_dummy);

        value
            .extra_data
            .momtrop_edge_power
            .as_ref()
            .map(|mep| e.add_statement("momtrop_edge_power", mep.to_canonical_string()));
        value
            .extra_data
            .vakint_edge_power
            .map(|vak| e.add_statement("vakint_edge_power", vak));

        e
    }
}

#[derive(Debug, Clone)]
pub struct ParseEdge {
    pub name: Option<String>,
    pub particle: PossibleParticle,
    pub dod: Option<i32>,
    pub is_dummy: bool,
    pub lmb_id: Option<LoopIndex>,
    /// User provided numerator
    pub num: Option<Atom>,
    /// Specifies the incoming initial state hedge
    pub is_cut: Option<Hedge>,
    pub momtrop_edge_power: Option<Atom>,
    pub vakint_edge_power: Option<isize>,
}

impl ParseEdge {
    pub fn from_symbolica_edge(
        model: &Model,
        edge_color: &EdgeColor,
        is_cut: Option<Hedge>,
    ) -> Self {
        let particle = model.get_particle_from_pdg(edge_color.pdg);
        let mut e = ParseEdge::new(particle);
        e.is_cut = is_cut;
        e
    }
}

impl Display for ParseEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        DotEdgeData::from(self).fmt(f)
    }
}

impl UVE for ParseEdge {
    fn mass_atom(&self) -> Atom {
        match &self.particle {
            PossibleParticle::JustMass { expr, .. } => expr.clone(),
            PossibleParticle::Particle(p) => p.mass.0.into(),
            PossibleParticle::MassOverriddenParticle { mass, .. } => mass.clone(),
        }
    }
}

impl ParseEdge {
    pub fn new(particle: impl Into<PossibleParticle>) -> Self {
        ParseEdge {
            is_dummy: false,
            name: None,
            particle: particle.into(),
            dod: None,
            lmb_id: None,
            num: None,
            is_cut: None,
            momtrop_edge_power: None,
            vakint_edge_power: None,
        }
    }

    pub fn is_dummy(mut self) -> Self {
        self.is_dummy = true;
        self
    }

    pub fn with_label(mut self, label: String) -> Self {
        self.name = Some(label);
        self
    }

    pub fn with_dod(mut self, dod: i32) -> Self {
        self.dod = Some(dod);
        self
    }

    pub fn with_lmb_id(mut self, lmb_id: LoopIndex) -> Self {
        self.lmb_id = Some(lmb_id);
        self
    }

    pub fn with_num(mut self, num: Atom) -> Self {
        self.num = Some(num);
        self
    }
}

impl ParseEdge {
    pub(crate) fn localize_ainds(
        atom: impl AtomCore,
        eid: EdgeIndex,
        hedge_pair: HedgePair,
    ) -> Atom {
        let a = atom
            .replace(GS.edgeid)
            .with(Atom::num(usize::from(eid) as i64))
            .replace_map(|term, _, out| {
                if let AtomView::Fun(f) = term
                    && f.get_symbol() == GS.edgeid
                    && f.get_nargs() == 1
                    && let Ok(i) = i64::try_from(f.iter().next().unwrap())
                    && let Ok(u) = u16::try_from(i)
                {
                    **out = eid.aind(u).into();
                }
            });

        match hedge_pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => a
                .replace(GS.sink_id)
                .with(Atom::num(sink.0 as i64))
                .replace_map(|term, _, out| {
                    if let AtomView::Fun(f) = term
                        && f.get_symbol() == GS.sink_id
                        && f.get_nargs() == 1
                        && let Ok(i) = i64::try_from(f.iter().next().unwrap())
                        && let Ok(u) = u16::try_from(i)
                    {
                        **out = sink.aind(u).into();
                    }
                })
                .replace(GS.source_id)
                .with(Atom::num(source.0 as i64))
                .replace_map(|term, _, out| {
                    if let AtomView::Fun(f) = term
                        && f.get_symbol() == GS.source_id
                        && f.get_nargs() == 1
                        && let Ok(i) = i64::try_from(f.iter().next().unwrap())
                        && let Ok(u) = u16::try_from(i)
                    {
                        **out = source.aind(u).into();
                    }
                }),
            HedgePair::Unpaired { hedge, flow } => match flow {
                Flow::Source => a
                    .replace(GS.source_id)
                    .with(Atom::num(hedge.0 as i64))
                    .replace_map(|term, _, out| {
                        if let AtomView::Fun(f) = term
                            && f.get_symbol() == GS.source_id
                            && f.get_nargs() == 1
                            && let Ok(i) = i64::try_from(f.iter().next().unwrap())
                            && let Ok(u) = u16::try_from(i)
                        {
                            **out = hedge.aind(u).into();
                        }
                    }),
                Flow::Sink => a
                    .replace(GS.sink_id)
                    .with(Atom::num(hedge.0 as i64))
                    .replace_map(|term, _, out| {
                        if let AtomView::Fun(f) = term
                            && f.get_symbol() == GS.sink_id
                            && f.get_nargs() == 1
                            && let Ok(i) = i64::try_from(f.iter().next().unwrap())
                            && let Ok(u) = u16::try_from(i)
                        {
                            **out = hedge.aind(u).into();
                        }
                    }),
            },
        }
    }
    pub(crate) fn parse<'a>(
        model: &'a Model,
    ) -> impl FnMut(
        &'a HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData>,
        EdgeIndex,
        HedgePair,
        EdgeData<&'a DotEdgeData>,
    ) -> Result<EdgeData<Self>> {
        |_: &'a HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData>,
         eid: EdgeIndex,
         p: HedgePair,
         e_data: EdgeData<&'a DotEdgeData>| {
            let e = e_data.data;
            let label = e.get::<_, String>("name").transpose()?;

            let lmb_id: Option<LoopIndex> = e
                .get::<_, usize>("lmb_id")
                .transpose()?
                .map(LoopIndex::from);

            let is_cut: Option<Hedge> = e.get::<_, usize>("is_cut").transpose()?.map(Hedge::from);

            let dod = e
                .get::<_, String>("dod")
                .transpose()
                .with_context(|| format!("Error parsing dod"))?
                .map(|a| a.strip_parse())
                .transpose()?;
            let is_dummy = e.get::<_, bool>("is_dummy").transpose()?.unwrap_or(false);

            let num = e
                .get::<_, String>("num")
                .transpose()?
                .map(|a| -> Result<Atom> {
                    Ok(Self::localize_ainds(a.strip_parse::<Atom>()?, eid, p))
                })
                .transpose()?;

            let mass = e
                .get::<_, String>("mass")
                .transpose()?
                .map(|a| a.strip_parse::<Atom>())
                .transpose()?;

            let momtrop_edge_power: Option<Atom> = e
                .get::<_, String>("momtrop_edge_power")
                .transpose()?
                .map(|a| parse!(&a));
            let vakint_edge_power: Option<isize> =
                e.get::<_, isize>("vakint_edge_power").transpose()?;

            let particle: PossibleParticle = if let Some(v) = e.get::<_, isize>("pdg") {
                model.try_get_particle_from_pdg(v?)?.into()
            } else if let Some(v) = e.get::<_, String>("particle") {
                let pname: String = v?.strip_parse()?;
                model.try_get_particle(pname)?.into()
            } else {
                ().into()
            };

            let orientation = if e.local_statements.contains_key("dir") {
                e_data.orientation
            } else {
                particle.orientation()
            };
            Ok(EdgeData::new(
                ParseEdge {
                    is_dummy,
                    dod,
                    lmb_id,
                    particle: particle.override_mass(mass),
                    num,
                    is_cut,
                    name: label,
                    momtrop_edge_power,
                    vakint_edge_power,
                },
                orientation,
            ))
        }
    }
}
