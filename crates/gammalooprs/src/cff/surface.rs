use super::{
    esurface::{Esurface, EsurfaceID},
    hsurface::{Hsurface, HsurfaceID},
};
use crate::utils::{GS, cut_energy, external_energy_atom_from_index, ose_atom_from_index};
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use linnet::half_edge::involution::EdgeIndex;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, fmt};
use symbolica::{atom::Atom, parse};
use typed_index_collections::TiVec;
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, Hash)]
pub struct RationalCoefficient {
    pub numerator: i64,
    pub denominator: i64,
}

impl RationalCoefficient {
    pub const fn zero() -> Self {
        Self {
            numerator: 0,
            denominator: 1,
        }
    }

    pub const fn one() -> Self {
        Self {
            numerator: 1,
            denominator: 1,
        }
    }

    pub fn new(numerator: i64, denominator: i64) -> Self {
        assert_ne!(denominator, 0, "rational denominator cannot be zero");
        let mut numerator = numerator;
        let mut denominator = denominator;
        if denominator < 0 {
            numerator = -numerator;
            denominator = -denominator;
        }
        let divisor = gcd_i64(numerator, denominator);
        Self {
            numerator: numerator / divisor,
            denominator: denominator / divisor,
        }
    }

    pub fn to_atom(self) -> Atom {
        Atom::num(self.numerator) / Atom::num(self.denominator)
    }

    pub fn is_zero(self) -> bool {
        self.numerator == 0
    }

    pub fn is_one(self) -> bool {
        self.numerator == self.denominator
    }
}

impl Default for RationalCoefficient {
    fn default() -> Self {
        Self::zero()
    }
}

impl fmt::Display for RationalCoefficient {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}

fn gcd_i64(lhs: i64, rhs: i64) -> i64 {
    let mut a = lhs.unsigned_abs();
    let mut b = rhs.unsigned_abs();
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    i64::try_from(a.max(1)).unwrap()
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, Hash, Default)]
pub struct LinearEnergyExpr {
    pub internal_terms: Vec<(EdgeIndex, i64)>,
    pub external_terms: Vec<(EdgeIndex, i64)>,
    pub uniform_scale_coeff: i64,
    pub constant: RationalCoefficient,
}

impl LinearEnergyExpr {
    pub fn zero() -> Self {
        Self::default()
    }

    pub fn ose(edge_id: EdgeIndex, coeff: i64) -> Self {
        if coeff == 0 {
            Self::zero()
        } else {
            Self {
                internal_terms: vec![(edge_id, coeff)],
                ..Self::default()
            }
        }
    }

    pub fn external(edge_id: EdgeIndex, coeff: i64) -> Self {
        if coeff == 0 {
            Self::zero()
        } else {
            Self {
                external_terms: vec![(edge_id, coeff)],
                ..Self::default()
            }
        }
    }

    pub fn uniform_scale(coeff: i64) -> Self {
        if coeff == 0 {
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
        self.constant =
            RationalCoefficient::new(self.constant.numerator, self.constant.denominator);
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
                acc + Atom::num(*coeff) * energy
            });

        let external = self
            .external_terms
            .iter()
            .fold(Atom::new(), |acc, (edge_id, coeff)| {
                acc + Atom::num(*coeff) * external_energy_atom_from_index(*edge_id)
            });

        let scale = if self.uniform_scale_coeff == 0 {
            Atom::new()
        } else {
            Atom::num(self.uniform_scale_coeff) * Atom::var(GS.numerator_sampling_scale)
        };

        internal + external + scale + self.constant.to_atom()
    }

    pub fn uses_uniform_scale(&self) -> bool {
        self.uniform_scale_coeff != 0
    }
}

impl fmt::Display for LinearEnergyExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms = Vec::<(String, bool)>::new();
        if !self.constant.is_zero() {
            let magnitude =
                RationalCoefficient::new(self.constant.numerator.abs(), self.constant.denominator);
            terms.push((magnitude.to_string(), self.constant.numerator < 0));
        }

        let mut push_term = |label: String, coeff: i64| {
            let magnitude = coeff.abs();
            let term = if magnitude == 1 {
                label
            } else {
                format!("{magnitude}*{label}")
            };
            terms.push((term, coeff < 0));
        };

        for (edge_id, coeff) in &self.internal_terms {
            push_term(format!("OSE[{}]", edge_id.0), *coeff);
        }
        for (edge_id, coeff) in &self.external_terms {
            push_term(format!("E[{}]", edge_id.0), *coeff);
        }
        if self.uniform_scale_coeff != 0 {
            push_term("M".to_string(), self.uniform_scale_coeff);
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

fn collect_linear_terms(terms: Vec<(EdgeIndex, i64)>) -> Vec<(EdgeIndex, i64)> {
    let mut collected = BTreeMap::<EdgeIndex, i64>::new();
    for (edge_id, coeff) in terms {
        *collected.entry(edge_id).or_default() += coeff;
    }
    collected
        .into_iter()
        .filter(|(_, coeff)| *coeff != 0)
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

impl LinearSurface {
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        self.expression.to_atom(cut_edges)
    }
}

#[derive(
    From, Into, Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Encode, Decode, Hash,
)]
pub struct LinearSurfaceID(pub usize);
pub type LinearSurfaceCollection = TiVec<LinearSurfaceID, LinearSurface>;

impl From<LinearSurfaceID> for Atom {
    fn from(value: LinearSurfaceID) -> Self {
        parse!(&format!("L({})", value.0))
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]

/// A esurface that is equal to 1, useful for represnting a single vertex
pub struct UnitSurface {}

/// Esurface whose inverse is equal to 0, useful for setting surfaces to zero
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct InfiniteSurface {}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HybridSurface {
    Esurface(Esurface),
    Hsurface(Hsurface),
    Linear(LinearSurface),
    Unit(UnitSurface),
    Infinite(InfiniteSurface),
}

impl HybridSurface {
    #[allow(dead_code)]
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        match self {
            HybridSurface::Esurface(surface) => surface.to_atom(cut_edges),
            HybridSurface::Hsurface(surface) => surface.to_atom(cut_edges),
            HybridSurface::Linear(surface) => surface.to_atom(cut_edges),
            HybridSurface::Unit(_) => Atom::num(1),
            HybridSurface::Infinite(_) => parse!("η_inf"),
        }
    }
}

#[derive(Debug, Clone)]
pub enum HybridSurfaceRef<'a> {
    Esurface(&'a Esurface),
    Hsurface(&'a Hsurface),
    Linear(&'a LinearSurface),
    Unit(UnitSurface),
    Infinite(InfiniteSurface),
}

impl HybridSurfaceRef<'_> {
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        match self {
            HybridSurfaceRef::Esurface(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Hsurface(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Linear(surface) => surface.to_atom(cut_edges),
            HybridSurfaceRef::Unit(_) => Atom::num(1),
            HybridSurfaceRef::Infinite(_) => parse!("η_inf"),
        }
    }
}

#[derive(From, Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Encode, Decode, Hash)]
pub enum HybridSurfaceID {
    Esurface(EsurfaceID),
    Hsurface(HsurfaceID),
    Linear(LinearSurfaceID),
    Unit,
    Infinite,
}

pub type HybridSurfaceCollection = TiVec<HybridSurfaceID, HybridSurface>;
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
