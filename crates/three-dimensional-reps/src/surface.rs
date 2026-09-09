use std::{
    collections::BTreeMap,
    fmt,
    ops::{Add, Mul, Neg, Sub},
};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::involution::EdgeIndex;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    id::{Pattern, Replacement},
    parse,
};
use typed_index_collections::TiVec;

use crate::symbols::{
    cut_energy, external_energy_atom_from_index, numerator_sampling_scale, ose_atom_from_index,
};
use crate::utils::Rational;

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Encode, Decode, Hash, PartialOrd, Ord,
)]
pub struct EsurfaceID(pub usize);

impl From<usize> for EsurfaceID {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<EsurfaceID> for usize {
    fn from(value: EsurfaceID) -> Self {
        value.0
    }
}

impl From<EsurfaceID> for Atom {
    fn from(id: EsurfaceID) -> Self {
        parse!(&format!("η({})", id.0))
    }
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Encode, Decode, Hash, PartialOrd, Ord,
)]
pub struct HsurfaceID(pub usize);

impl From<usize> for HsurfaceID {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<HsurfaceID> for usize {
    fn from(value: HsurfaceID) -> Self {
        value.0
    }
}

impl From<HsurfaceID> for Atom {
    fn from(value: HsurfaceID) -> Self {
        parse!(&format!("H({})", value.0))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, Hash, Default)]
pub struct LinearEnergyExpr {
    pub internal_terms: Vec<(EdgeIndex, Rational)>,
    pub external_terms: Vec<(EdgeIndex, Rational)>,
    pub uniform_scale_coeff: Rational,
    pub constant: Rational,
}

impl LinearEnergyExpr {
    pub fn zero() -> Self {
        Self::default()
    }

    pub fn ose(edge_id: EdgeIndex, coeff: i64) -> Self {
        Self::ose_with_coeff(edge_id, Rational::from(coeff))
    }

    pub fn external(edge_id: EdgeIndex, coeff: i64) -> Self {
        Self::external_with_coeff(edge_id, Rational::from(coeff))
    }

    pub fn uniform_scale(coeff: i64) -> Self {
        Self::uniform_scale_with_coeff(Rational::from(coeff))
    }

    pub fn ose_with_coeff(edge_id: EdgeIndex, coeff: Rational) -> Self {
        if coeff.is_zero() {
            Self::zero()
        } else {
            Self {
                internal_terms: vec![(edge_id, coeff)],
                ..Self::default()
            }
        }
    }

    pub fn external_with_coeff(edge_id: EdgeIndex, coeff: Rational) -> Self {
        if coeff.is_zero() {
            Self::zero()
        } else {
            Self {
                external_terms: vec![(edge_id, coeff)],
                ..Self::default()
            }
        }
    }

    pub fn uniform_scale_with_coeff(coeff: Rational) -> Self {
        if coeff.is_zero() {
            Self::zero()
        } else {
            Self {
                uniform_scale_coeff: coeff,
                ..Self::default()
            }
        }
    }

    pub fn canonical(mut self) -> Self {
        self.internal_terms = collect_linear_terms(self.internal_terms);
        self.external_terms = collect_linear_terms(self.external_terms);
        self
    }

    pub fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        let internal = self
            .internal_terms
            .iter()
            .fold(Atom::new(), |acc, (edge_id, coeff)| {
                let energy = if cut_edges.contains(edge_id) {
                    cut_energy(*edge_id)
                } else {
                    ose_atom_from_index(*edge_id)
                };
                acc + Atom::num(coeff.clone()) * energy
            });

        let external = self
            .external_terms
            .iter()
            .fold(Atom::new(), |acc, (edge_id, coeff)| {
                acc + Atom::num(coeff.clone()) * external_energy_atom_from_index(*edge_id)
            });

        let scale = if self.uniform_scale_coeff.is_zero() {
            Atom::new()
        } else {
            Atom::num(self.uniform_scale_coeff.clone()) * numerator_sampling_scale()
        };

        internal + external + scale + Atom::num(self.constant.clone())
    }

    pub fn uses_uniform_scale(&self) -> bool {
        !self.uniform_scale_coeff.is_zero()
    }

    pub fn is_zero(&self) -> bool {
        self.internal_terms.is_empty()
            && self.external_terms.is_empty()
            && self.uniform_scale_coeff.is_zero()
            && self.constant.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.internal_terms.is_empty()
            && self.external_terms.is_empty()
            && self.uniform_scale_coeff.is_zero()
            && self.constant.is_one()
    }

    pub fn remap_internal_edges(self, edge_map: &std::collections::BTreeMap<usize, usize>) -> Self {
        self.remap_energy_edges(edge_map, &std::collections::BTreeMap::new())
    }

    pub fn remap_energy_edges(
        self,
        internal_edge_map: &std::collections::BTreeMap<usize, usize>,
        external_edge_map: &std::collections::BTreeMap<usize, usize>,
    ) -> Self {
        Self {
            internal_terms: self
                .internal_terms
                .into_iter()
                .map(|(edge_id, coeff)| {
                    (
                        EdgeIndex(
                            internal_edge_map
                                .get(&edge_id.0)
                                .copied()
                                .unwrap_or(edge_id.0),
                        ),
                        coeff,
                    )
                })
                .collect(),
            external_terms: self
                .external_terms
                .into_iter()
                .map(|(edge_id, coeff)| {
                    (
                        EdgeIndex(
                            external_edge_map
                                .get(&edge_id.0)
                                .copied()
                                .unwrap_or(edge_id.0),
                        ),
                        coeff,
                    )
                })
                .collect(),
            uniform_scale_coeff: self.uniform_scale_coeff,
            constant: self.constant,
        }
        .canonical()
    }

    pub fn scale(self, coeff: i64) -> Self {
        self.scale_rational(Rational::from(coeff))
    }

    pub fn scale_rational(mut self, coeff: Rational) -> Self {
        if coeff.is_zero() {
            return Self::zero();
        }
        for (_, item_coeff) in self
            .internal_terms
            .iter_mut()
            .chain(&mut self.external_terms)
        {
            *item_coeff *= &coeff;
        }
        self.uniform_scale_coeff *= &coeff;
        self.constant *= coeff;
        self.canonical()
    }

    pub fn substitute_external_energy(
        mut self,
        edge_id: EdgeIndex,
        replacement: &LinearEnergyExpr,
    ) -> Self {
        let mut replacement_terms = LinearEnergyExpr::zero();
        let mut retained_external_terms = Vec::with_capacity(self.external_terms.len());

        for (term_edge_id, coeff) in self.external_terms {
            if term_edge_id == edge_id {
                replacement_terms = replacement_terms + replacement.clone().scale_rational(coeff);
            } else {
                retained_external_terms.push((term_edge_id, coeff));
            }
        }

        self.external_terms = retained_external_terms;
        (self + replacement_terms).canonical()
    }
}

impl Add for LinearEnergyExpr {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.internal_terms.extend(rhs.internal_terms);
        self.external_terms.extend(rhs.external_terms);
        self.uniform_scale_coeff += rhs.uniform_scale_coeff;
        self.constant += rhs.constant;
        self.canonical()
    }
}

impl Sub for LinearEnergyExpr {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Neg for LinearEnergyExpr {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.scale(-1)
    }
}

impl Mul<i64> for LinearEnergyExpr {
    type Output = Self;

    fn mul(self, rhs: i64) -> Self::Output {
        self.scale(rhs)
    }
}

impl fmt::Display for LinearEnergyExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms = Vec::<(String, bool)>::new();
        if !self.constant.is_zero() {
            let magnitude = self.constant.abs();
            terms.push((magnitude.to_string(), self.constant.is_negative()));
        }

        let mut push_term = |label: String, coeff: &Rational| {
            let magnitude = coeff.abs();
            let term = if magnitude.is_one() {
                label
            } else {
                format!("{magnitude}*{label}")
            };
            terms.push((term, coeff.is_negative()));
        };

        for (edge_id, coeff) in &self.internal_terms {
            push_term(format!("OSE[{}]", edge_id.0), coeff);
        }
        for (edge_id, coeff) in &self.external_terms {
            push_term(format!("E[{}]", edge_id.0), coeff);
        }
        if !self.uniform_scale_coeff.is_zero() {
            push_term("M".to_string(), &self.uniform_scale_coeff);
        }

        let Some((first, negative)) = terms.first() else {
            return f.write_str("0");
        };

        if *negative {
            write!(f, "- {first}")?;
        } else {
            f.write_str(first)?;
        }
        for (term, negative) in terms.iter().skip(1) {
            if *negative {
                write!(f, " - {term}")?;
            } else {
                write!(f, " + {term}")?;
            }
        }
        Ok(())
    }
}

fn collect_linear_terms(terms: Vec<(EdgeIndex, Rational)>) -> Vec<(EdgeIndex, Rational)> {
    let mut collected = BTreeMap::<EdgeIndex, Rational>::new();
    for (edge_id, coeff) in terms {
        *collected.entry(edge_id).or_default() += coeff;
    }
    collected
        .into_iter()
        .filter(|(_, coeff)| !coeff.is_zero())
        .collect()
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, Hash)]
pub enum LinearSurfaceKind {
    Esurface,
    Hsurface,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, Hash)]
pub enum SurfaceOrigin {
    Physical,
    Helper,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, Hash)]
pub struct LinearSurface {
    pub kind: LinearSurfaceKind,
    pub expression: LinearEnergyExpr,
    pub origin: SurfaceOrigin,
    pub numerator_only: bool,
}

impl SurfaceAtom for LinearSurface {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        self.expression.to_atom(cut_edges)
    }
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Encode, Decode, Hash, PartialOrd, Ord,
)]
pub struct LinearSurfaceID(pub usize);

impl From<usize> for LinearSurfaceID {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<LinearSurfaceID> for usize {
    fn from(value: LinearSurfaceID) -> Self {
        value.0
    }
}

pub type LinearSurfaceCollection = TiVec<LinearSurfaceID, LinearSurface>;

impl From<LinearSurfaceID> for Atom {
    fn from(value: LinearSurfaceID) -> Self {
        parse!(&format!("L({})", value.0))
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct UnitSurface {}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct InfiniteSurface {}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HybridSurface<E = (), H = ()> {
    Esurface(E),
    Hsurface(H),
    Linear(LinearSurface),
    Unit(UnitSurface),
    Infinite(InfiniteSurface),
}

#[derive(Debug, Clone)]
pub enum HybridSurfaceRef<'a, E = (), H = ()> {
    Esurface(&'a E),
    Hsurface(&'a H),
    Linear(&'a LinearSurface),
    Unit(UnitSurface),
    Infinite(InfiniteSurface),
}

impl<E: SurfaceAtom, H: SurfaceAtom> HybridSurfaceRef<'_, E, H> {
    pub fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        match self {
            HybridSurfaceRef::Esurface(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Hsurface(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Linear(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Unit(_) => Atom::num(1),
            HybridSurfaceRef::Infinite(_) => parse!("η_inf"),
        }
    }
}

#[derive(
    Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Encode, Decode, Hash, PartialOrd, Ord,
)]
pub enum HybridSurfaceID {
    Esurface(EsurfaceID),
    Hsurface(HsurfaceID),
    Linear(LinearSurfaceID),
    Unit,
    Infinite,
}

pub type EsurfaceCollection<E> = TiVec<EsurfaceID, E>;
pub type HsurfaceCollection<H> = TiVec<HsurfaceID, H>;
pub type HsurfaceCache<T> = TiVec<HsurfaceID, T>;
pub type HybridSurfaceCollection<E, H> = TiVec<HybridSurfaceID, HybridSurface<E, H>>;
pub type HybridSurfaceCache<T> = TiVec<HybridSurfaceID, T>;

impl From<HybridSurfaceID> for Atom {
    fn from(id: HybridSurfaceID) -> Atom {
        match id {
            HybridSurfaceID::Esurface(id) => Atom::from(id),
            HybridSurfaceID::Hsurface(id) => Atom::from(id),
            HybridSurfaceID::Linear(id) => Atom::from(id),
            HybridSurfaceID::Unit => Atom::num(1),
            HybridSurfaceID::Infinite => parse!("η_inf"),
        }
    }
}

pub trait SurfaceAtom {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom;
}

impl SurfaceAtom for () {
    fn to_atom(&self, _cut_edges: &[EdgeIndex]) -> Atom {
        Atom::new()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct SurfaceCache<E = (), H = ()> {
    pub esurface_cache: EsurfaceCollection<E>,
    pub hsurface_cache: HsurfaceCollection<H>,
    pub linear_surface_cache: LinearSurfaceCollection,
}

impl<E, H> SurfaceCache<E, H> {
    pub fn new() -> Self {
        Self {
            esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
            hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
            linear_surface_cache: LinearSurfaceCollection::from_iter(std::iter::empty()),
        }
    }
}

impl<E, H> Default for SurfaceCache<E, H> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E: SurfaceAtom, H: SurfaceAtom> SurfaceCache<E, H> {
    pub fn substitute_energies(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom {
        let replacement_rules = self.get_all_replacements(cut_edges);
        atom.replace_multiple(&replacement_rules)
    }

    pub fn iter_all_surfaces(
        &'_ self,
    ) -> impl Iterator<Item = (HybridSurfaceID, HybridSurfaceRef<'_, E, H>)> + '_ {
        let esurface_id_iter = self.esurface_cache.iter_enumerated().map(|(id, esurface)| {
            (
                HybridSurfaceID::Esurface(id),
                HybridSurfaceRef::Esurface(esurface),
            )
        });

        let hsurface_id_iter = self.hsurface_cache.iter_enumerated().map(|(id, hsurface)| {
            (
                HybridSurfaceID::Hsurface(id),
                HybridSurfaceRef::Hsurface(hsurface),
            )
        });

        let linear_surface_id_iter =
            self.linear_surface_cache
                .iter_enumerated()
                .map(|(id, surface)| {
                    (
                        HybridSurfaceID::Linear(id),
                        HybridSurfaceRef::Linear(surface),
                    )
                });

        esurface_id_iter
            .chain(hsurface_id_iter)
            .chain(linear_surface_id_iter)
    }

    pub fn get_all_replacements(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(surface.to_atom(cut_edges));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }

    pub fn get_surface(&self, surface_id: HybridSurfaceID) -> HybridSurfaceRef<'_, E, H> {
        match surface_id {
            HybridSurfaceID::Esurface(id) => HybridSurfaceRef::Esurface(&self.esurface_cache[id]),
            HybridSurfaceID::Hsurface(id) => HybridSurfaceRef::Hsurface(&self.hsurface_cache[id]),
            HybridSurfaceID::Linear(id) => HybridSurfaceRef::Linear(&self.linear_surface_cache[id]),
            HybridSurfaceID::Unit => HybridSurfaceRef::Unit(UnitSurface {}),
            HybridSurfaceID::Infinite => HybridSurfaceRef::Infinite(InfiniteSurface {}),
        }
    }
}
