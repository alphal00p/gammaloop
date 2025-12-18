#![allow(dead_code)]

use std::borrow::Borrow;

use crate::{
    cff::esurface::EsurfaceID,
    settings::global::OrientationPattern,
    utils::{W_, ose_atom_from_index},
};
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use itertools::{EitherOrBoth, Itertools};
use linnet::half_edge::{
    GVEdgeAttrs, HedgeGraph,
    involution::{EdgeVec, Orientation},
    nodestore::NodeStorageOps,
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol},
    function,
    id::{Pattern, Replacement},
    symbol,
};
use typed_index_collections::TiVec;

use super::{
    cut_expression::SuperGraphOrientationID, generation::SurfaceCache, surface::HybridSurfaceID,
    tree::Tree,
};

use crate::utils::GS;

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    From,
    Into,
    Hash,
    PartialEq,
    Eq,
    Copy,
    Encode,
    Decode,
    PartialOrd,
    Ord,
)]
pub struct AmplitudeOrientationID(pub usize);

#[derive(
    Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq, Copy, Encode, Decode,
)]
pub struct SubgraphOrientationID(pub usize);

impl GraphOrientation for EdgeVec<Orientation> {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        self
    }
}

pub trait GraphOrientation {
    fn orientation(&self) -> &EdgeVec<Orientation>;

    fn orientation_thetas(&self) -> Atom {
        let mut thetas = Atom::num(1);

        for (e, h) in self.orientation() {
            match h {
                Orientation::Default => {
                    thetas *= GS.sign_theta(GS.sign(e));
                }
                Orientation::Reversed => {
                    thetas *= GS.sign_theta(-GS.sign(e));
                }
                _ => {}
            }
        }
        thetas
    }

    fn orientation_delta(&self) -> Atom
    where
        Self: Sized,
    {
        GS.orientation_delta(self)
    }

    fn select<'a>(&self, atom: impl Into<AtomOrView<'a>>) -> Atom {
        let theta_reps = vec![
            Replacement::new(GS.sign_theta(1).to_pattern(), Atom::num(1)),
            Replacement::new(GS.sign_theta(-1).to_pattern(), Atom::Zero),
        ];

        let mut reps = Vec::new();

        for (e, h) in self.orientation() {
            match h {
                Orientation::Default => {
                    reps.push(Replacement::new(GS.sign(e).to_pattern(), Atom::num(1)));
                }
                Orientation::Reversed => {
                    reps.push(Replacement::new(GS.sign(e).to_pattern(), Atom::num(-1)));
                }
                _ => {}
            }
        }

        let orientation = self.orientation();
        atom.into()
            .replace_multiple(&reps)
            .replace_multiple(&theta_reps)
            .replace_map(|term, _ctx, out| {
                if let AtomView::Fun(f) = term {
                    if f.get_symbol() == GS.orientation_delta {
                        if f.iter().zip_longest(orientation.into_iter()).all(
                            |either| match either {
                                EitherOrBoth::Both(a, (_, o)) => {
                                    if let Ok(a) = i64::try_from(a) {
                                        match o {
                                            Orientation::Default => a >= 0,
                                            Orientation::Reversed => a <= 0,
                                            Orientation::Undirected => true,
                                        }
                                    } else {
                                        false
                                    }
                                }
                                EitherOrBoth::Left(_) | EitherOrBoth::Right(_) => false,
                            },
                        ) {
                            **out = Atom::num(1);
                        } else {
                            **out = Atom::Zero;
                        }
                    }
                }
            })
    }
}

pub trait OrientationID: From<usize> + Into<usize> {
    fn symbol() -> Symbol;

    fn atom(self) -> Atom {
        let id: usize = self.into();
        function!(Self::symbol(), id as i64)
    }

    fn select<'a>(self, atom: impl Into<AtomOrView<'a>>) -> Atom {
        atom.into()
            .as_view()
            .replace(self.atom())
            .with(Atom::num(1))
            .replace(function!(Self::symbol(), W_.x_))
            .with(Atom::Zero)
    }
}

impl OrientationID for AmplitudeOrientationID {
    fn symbol() -> Symbol {
        symbol!("amp_sigma")
    }
}

impl OrientationID for SubgraphOrientationID {
    fn symbol() -> Symbol {
        symbol!("subg_sigma")
    }
}

impl OrientationID for SuperGraphOrientationID {
    fn symbol() -> Symbol {
        symbol!("superg_sigma")
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationData {
    pub orientation: EdgeVec<Orientation>,
}

impl GraphOrientation for OrientationData {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.orientation
    }
}

impl OrientationData {
    pub(crate) fn dot<E, V, N: NodeStorageOps<NodeData = V>>(&self, graph: &HedgeGraph<E, V, N>) {
        let mut writer = String::new();
        writer.push_str("digraph {{");

        writer.push_str(&format!(
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";",
        ));

        for (hedge_pair, id, _) in graph.iter_edges() {
            let attr = GVEdgeAttrs {
                color: None,
                label: None,
                other: None,
            };
            writer.push_str("  ");

            let attr = hedge_pair.fill_color(attr);
            hedge_pair
                .add_data(graph)
                .dot_fmt(
                    &mut writer,
                    graph,
                    id,
                    |_| None,
                    |a| a.to_string(),
                    self.orientation[id],
                    attr,
                )
                .unwrap();
        }
        writer.push_str("}}");
    }

    pub(crate) fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.orientation
            .borrow()
            .into_iter()
            .filter_map(|(edge_index, orientation)| {
                if matches!(orientation, Orientation::Reversed) {
                    let energy_atom = ose_atom_from_index(edge_index);

                    let neg_energy_atom = -&energy_atom;

                    Some(Replacement::new(
                        Pattern::from(energy_atom),
                        Pattern::from(neg_energy_atom),
                    ))
                } else {
                    None
                }
            })
            .collect()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationExpression {
    pub data: OrientationData,
    pub expression: Tree<HybridSurfaceID>,
}
#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct CFFExpression<O: OrientationID> {
    pub orientations: TiVec<O, OrientationExpression>,
    pub surfaces: SurfaceCache,
}

impl GraphOrientation for OrientationExpression {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.data.orientation
    }
}

impl<O: OrientationID> CFFExpression<O>
where
    usize: From<O>,
{
    pub(crate) fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
            surfaces: SurfaceCache {
                esurface_cache: TiVec::new(),
                hsurface_cache: TiVec::new(),
            },
        }
    }

    pub fn to_atom(&self, pattern: OrientationPattern) -> Atom {
        self.orientations
            .iter()
            .filter_map(|orientation| {
                if pattern.alt_filter(orientation.orientation()) {
                    Some(orientation.expression.to_atom_inv())
                } else {
                    None
                }
            })
            .reduce(|a, b| a + b)
            .unwrap_or_default()
    }

    pub(crate) fn get_orientation_atoms(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<AmplitudeOrientationID, Atom> {
        self.orientations
            .iter()
            .map(|orientation| {
                if pattern.alt_filter(orientation.orientation()) {
                    orientation.expression.to_atom_inv()
                } else {
                    Atom::new()
                }
            })
            .collect()
    }

    pub fn get_orientation_atoms_with_data(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<AmplitudeOrientationID, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = if pattern.alt_filter(orientation.orientation()) {
                    orientation.expression.to_atom_inv()
                } else {
                    Atom::new()
                };
                let data = orientation.data.clone();
                (atom, data)
            })
            .collect()
    }

    pub(crate) fn get_orientation_atom(&self, orientation_id: O) -> Atom {
        self.orientations[orientation_id].expression.to_atom_inv()
    }

    pub(crate) fn num_unfolded_terms(&self) -> usize {
        self.orientations
            .iter()
            .map(|o| o.expression.get_bottom_layer().len())
            .sum()
    }

    pub(crate) fn get_orientations_with_esurface(&self, esurface_id: EsurfaceID) -> Vec<O> {
        self.orientations
            .iter_enumerated()
            .filter_map(|(id, orientation)| {
                if orientation
                    .expression
                    .iter_nodes()
                    .any(|node| node.data == HybridSurfaceID::Esurface(esurface_id))
                {
                    Some(id)
                } else {
                    None
                }
            })
            .collect()
    }
}

impl From<CFFExpression<AmplitudeOrientationID>> for CFFExpression<SubgraphOrientationID> {
    fn from(value: CFFExpression<AmplitudeOrientationID>) -> Self {
        Self {
            orientations: value.orientations.raw.into(),
            surfaces: value.surfaces,
        }
    }
}

impl From<CFFExpression<SubgraphOrientationID>> for CFFExpression<AmplitudeOrientationID> {
    fn from(value: CFFExpression<SubgraphOrientationID>) -> Self {
        Self {
            orientations: value.orientations.raw.into(),
            surfaces: value.surfaces,
        }
    }
}
