use super::JetAlgorithm;
use super::types::{ClusteringResult, PseudoJet};
use crate::utils::{F, FloatLike};
use feynkit_kinematics::JetDefinition;
use smallvec::SmallVec;

pub(crate) fn cluster_candidates<T: FloatLike>(
    algorithm: JetAlgorithm,
    d_r: f64,
    min_jpt: f64,
    active: SmallVec<[PseudoJet<T>; 8]>,
) -> ClusteringResult<T> {
    if active.is_empty() {
        return ClusteringResult::default();
    }

    // Runtime settings currently store clustering parameters as `f64`, so this is
    // the intentional configuration precision boundary for these two values.
    let min_pt = F::from_f64(min_jpt);

    if d_r <= 0.0 || !d_r.is_finite() {
        return unrecombined_candidates(active, &min_pt);
    }

    let definition =
        JetDefinition::new(algorithm.into(), F::from_f64(d_r)).with_minimum_pt(min_pt.clone());
    match definition.cluster_indexed(active.iter().cloned().map(PseudoJet::into_feynkit)) {
        Ok(result) => ClusteringResult::from_feynkit(result),
        // GammaLoop event processing predates fallible clustering. Reuse its
        // existing unrecombined-candidate policy for malformed event data
        // instead of panicking or reclassifying the whole event as zero jets.
        Err(_) => unrecombined_candidates(active, &min_pt),
    }
}

fn unrecombined_candidates<T: FloatLike>(
    mut active: SmallVec<[PseudoJet<T>; 8]>,
    min_pt: &F<T>,
) -> ClusteringResult<T> {
    let mut result = ClusteringResult {
        jets: active
            .drain(..)
            .filter(|jet| jet.pt() >= min_pt.clone())
            .map(PseudoJet::jet)
            .collect(),
    };
    result.sort_by_pt();
    result
}
