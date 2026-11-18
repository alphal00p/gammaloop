//! Fixtures shared by the uv unit tests: the thermal rings of
//! `tests/resources/graphs` and the generation settings of the direct
//! local-3D route.

use std::collections::BTreeSet;

use color_eyre::Result;
use linnet::half_edge::involution::EdgeIndex;

use crate::graph::Graph;
use crate::graph::parse::IntoGraph;
use crate::settings::global::{MediumMode, MediumSettings};
use crate::utils::load_generic_model;
use crate::uv::{ApproximationType, RenormalizationPrescriptionSettings, UVgenerationSettings};

/// The thermal rings of `tests/resources/graphs`: two polarization bubbles
/// A = {e2, e3} and B = {e4, e5} (T = {e4} for the tadpole) joined by the
/// gluons e0 and e1, which carry q, with loop momenta [q, p, k] = [e0, e2, e4].
/// They are named after their second loop: the quark ring, the mixed ring of
/// a quark and a gluon bubble, the ghost ring and the tadpole ring.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(in crate::uv) enum Ring {
    Quark,
    Mixed,
    Ghost,
    Tadpole,
}

impl Ring {
    pub(in crate::uv) const ALL: [Self; 4] = [Self::Quark, Self::Mixed, Self::Ghost, Self::Tadpole];

    /// The ring with the Feynman-rule numerators of the generic `sm` model.
    pub(in crate::uv) fn physical_graph(self) -> Result<Graph> {
        let dot = match self {
            Self::Quark => {
                include_str!("../../../../tests/resources/graphs/thermal_quark_ring.dot")
            }
            Self::Mixed => {
                include_str!("../../../../tests/resources/graphs/thermal_gluon_quark_ring.dot")
            }
            Self::Ghost => {
                include_str!("../../../../tests/resources/graphs/thermal_ghost_quark_ring.dot")
            }
            Self::Tadpole => {
                include_str!("../../../../tests/resources/graphs/thermal_tadpole_quark_ring.dot")
            }
        };
        let graph: Graph = dot.into_graph(&load_generic_model("sm"))?;
        assert_eq!(
            graph.loop_momentum_basis.loop_edges.raw,
            [0, 2, 4].map(EdgeIndex)
        );
        Ok(graph)
    }

    /// The ring topology with unit numerators.
    pub(in crate::uv) fn topology(self) -> Result<Graph> {
        let second = match self {
            Self::Quark => "b0 -> b1 [id=4 particle=d lmb_id=2]; b1 -> b0 [id=5 particle=d];",
            Self::Mixed => {
                "b0 -> b1 [id=4 dir=none particle=g lmb_id=2]; b1 -> b0 [id=5 dir=none particle=g];"
            }
            Self::Ghost => "b0 -> b1 [id=4 particle=ghG lmb_id=2]; b1 -> b0 [id=5 particle=ghG];",
            Self::Tadpole => "b0 -> b0 [id=4 dir=none particle=g lmb_id=2];",
        };
        let close = if matches!(self, Self::Tadpole) {
            "b0"
        } else {
            "b1"
        };
        format!(
            "digraph zero_mode_ring {{
                node [num=1]; edge [num=1];
                a1 -> b0 [id=0 dir=none particle=g lmb_id=0];
                {close} -> a0 [id=1 dir=none particle=g];
                a0 -> a1 [id=2 particle=d lmb_id=1];
                a1 -> a0 [id=3 particle=d];
                {second}
            }}"
        )
        .into_graph(&load_generic_model("sm"))
    }

    /// The edges of the second polarization loop.
    pub(in crate::uv) fn second_loop_edges(self) -> &'static [usize] {
        match self {
            Self::Tadpole => &[4],
            Self::Quark | Self::Mixed | Self::Ghost => &[4, 5],
        }
    }

    /// The class representatives whose energies the physical numerator
    /// names: the three-gluon vertices carry q and k, the ghost vertices only
    /// the ghost momenta k and k - q, and the quark-gluon and four-gluon
    /// vertices none.
    pub(in crate::uv) fn numerator_class_energies(self) -> &'static [usize] {
        match self {
            Self::Mixed => &[0, 4],
            Self::Ghost => &[4],
            Self::Quark | Self::Tadpole => &[],
        }
    }
}

/// The edges with the ids `edges`.
pub(in crate::uv) fn edge_set(edges: &[usize]) -> BTreeSet<EdgeIndex> {
    edges.iter().copied().map(EdgeIndex).collect()
}

/// The direct local-3D generation settings without integrated counterterms.
pub(in crate::uv) fn local_3d_settings() -> UVgenerationSettings {
    UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    }
}

/// The settings of the ring run cards with one prescription for every
/// degree: the direct local-3D hedge-poset route without integrated
/// counterterms.
pub(in crate::uv) fn ring_uv_settings(prescription: ApproximationType) -> UVgenerationSettings {
    UVgenerationSettings {
        renormalization_prescription: RenormalizationPrescriptionSettings {
            log_divergent: prescription,
            massive_power_divergent: prescription,
            massless_power_divergent: prescription,
            ..Default::default()
        },
        ..local_3d_settings()
    }
}

pub(in crate::uv) fn medium(mode: MediumMode, vacuum_subtraction: bool) -> MediumSettings {
    MediumSettings {
        mode,
        vacuum_subtraction,
    }
}

pub(in crate::uv) fn thermal_medium(vacuum_subtraction: bool) -> MediumSettings {
    medium(MediumMode::ThermodynamicEquilibrium, vacuum_subtraction)
}
