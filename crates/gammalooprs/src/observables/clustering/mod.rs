mod genkt;
mod types;

use super::events::GenericEvent;
use crate::utils::FloatLike;
use bincode_trait_derive::{Decode, Encode};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

pub use types::{ClusteringResult, Jet};

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum JetAlgorithm {
    Kt,
    CambridgeAachen,
    #[default]
    AntiKt,
}

#[derive(Debug, Clone, Default)]
pub struct JetClustering {
    algorithm: JetAlgorithm,
    clustered_pdgs: Vec<isize>,
    d_r: f64,
    min_jpt: f64,
}

impl JetClustering {
    pub fn new(
        algorithm: JetAlgorithm,
        d_r: f64,
        min_jpt: f64,
        clustered_pdgs: Vec<isize>,
    ) -> Self {
        Self {
            algorithm,
            clustered_pdgs,
            d_r,
            min_jpt,
        }
    }

    pub fn process_event<T: FloatLike>(&self, event: &GenericEvent<T>) -> ClusteringResult<T> {
        let candidates = event
            .cut_info
            .particle_pdgs
            .1
            .iter()
            .copied()
            .zip(event.kinematic_configuration.1.iter().cloned())
            .enumerate()
            .filter(|(_, (pdg, _))| is_qcd_like_jet_pdg(&self.clustered_pdgs, *pdg))
            .map(|(index, (_pdg, momentum))| types::PseudoJet::from_momentum(index, momentum))
            .collect::<SmallVec<[_; 8]>>();
        self.cluster_candidates(candidates)
    }

    pub fn cluster_momenta<T: FloatLike>(
        &self,
        momenta: &[crate::momentum::FourMomentum<crate::utils::F<T>>],
    ) -> ClusteringResult<T> {
        let candidates = momenta
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, momentum)| types::PseudoJet::from_momentum(index, momentum))
            .collect::<SmallVec<[_; 8]>>();
        self.cluster_candidates(candidates)
    }

    fn cluster_candidates<T: FloatLike>(
        &self,
        candidates: SmallVec<[types::PseudoJet<T>; 8]>,
    ) -> ClusteringResult<T> {
        genkt::cluster_candidates(self.algorithm, self.d_r, self.min_jpt, candidates)
    }
}

fn is_qcd_like_jet_pdg(clustered_pdgs: &[isize], pdg: isize) -> bool {
    clustered_pdgs.binary_search(&pdg).is_ok()
}
