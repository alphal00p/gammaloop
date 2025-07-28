use std::borrow::Borrow;

use crate::{
    cff::esurface::EsurfaceID,
    utils::{ose_atom_from_index, W_},
};
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use itertools::{EitherOrBoth, Itertools};
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    nodestore::NodeStorageOps,
    GVEdgeAttrs, HedgeGraph,
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
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
    Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq, Copy, Encode, Decode,
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

    fn orientation_delta(&self) -> Atom {
        let mut fnbld = FunctionBuilder::new(GS.sign_delta);

        for (_, h) in self.orientation() {
            match h {
                Orientation::Default => fnbld = fnbld.add_arg(1),
                Orientation::Reversed => fnbld = fnbld.add_arg(-1),
                Orientation::Undirected => fnbld = fnbld.add_arg(0),
            }
        }
        fnbld.finish()
    }
    fn select<'a>(&self, atom: impl Into<AtomOrView<'a>>) -> Atom {
        let orientation = self.orientation();
        atom.into().as_view().replace_map(|term, _ctx, out| {
            if let AtomView::Fun(f) = term {
                if f.get_symbol() == GS.sign_delta {
                    if f.iter()
                        .zip_longest(orientation.into_iter())
                        .all(|either| match either {
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
                        })
                    {
                        *out = Atom::num(1);
                    } else {
                        *out = Atom::Zero;
                    }
                    return true;
                }

                if f.get_symbol() == GS.sign {
                    if f.get_nargs() == 1 {
                        let arg = f.iter().next().unwrap();
                        if let Ok(a) = i64::try_from(arg) {
                            if let Ok(a) = usize::try_from(a) {
                                match orientation[EdgeIndex::from(a)] {
                                    Orientation::Default => {
                                        *out = Atom::num(1);
                                        return true;
                                    }
                                    Orientation::Reversed => {
                                        *out = Atom::num(-1);
                                        return true;
                                    }
                                    Orientation::Undirected => {}
                                }
                            }
                        }
                    }
                }
            }
            false
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
    pub fn dot<E, V, N: NodeStorageOps<NodeData = V>>(&self, graph: &HedgeGraph<E, V, N>) {
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
                .dot_fmt(&mut writer, graph, id, |_| None, self.orientation[id], attr)
                .unwrap();
        }
        writer.push_str("}}");
    }

    pub fn get_ose_replacements(&self) -> Vec<Replacement> {
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
    pub fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
            surfaces: SurfaceCache {
                esurface_cache: TiVec::new(),
                hsurface_cache: TiVec::new(),
            },
        }
    }

    pub fn to_atom(&self) -> Atom {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.to_atom_inv())
            .reduce(|a, b| a + b)
            .unwrap_or_default()
    }

    pub fn get_orientation_atoms(&self) -> TiVec<AmplitudeOrientationID, Atom> {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.to_atom_inv())
            .collect()
    }

    pub fn get_orientation_atoms_with_data(
        &self,
    ) -> TiVec<AmplitudeOrientationID, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = orientation.expression.to_atom_inv();
                let data = orientation.data.clone();
                (atom, data)
            })
            .collect()
    }

    pub fn get_orientation_atom(&self, orientation_id: O) -> Atom {
        self.orientations[orientation_id].expression.to_atom_inv()
    }

    pub fn num_unfolded_terms(&self) -> usize {
        self.orientations
            .iter()
            .map(|o| o.expression.get_bottom_layer().len())
            .sum()
    }

    pub fn get_orientations_with_esurface(&self, esurface_id: EsurfaceID) -> Vec<O> {
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
