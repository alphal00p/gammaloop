use super::JetAlgorithm;
use super::types::{ClusteringResult, Jet, PseudoJet};
use crate::utils::{F, FloatLike};
use smallvec::SmallVec;
use std::cmp::Ordering;
use symbolica::domains::float::Real;

#[derive(Clone, Copy)]
enum BestAction {
    Beam(usize),
    Pair(usize, usize),
}

pub(crate) fn cluster_candidates<T: FloatLike>(
    algorithm: JetAlgorithm,
    d_r: f64,
    min_jpt: f64,
    mut active: SmallVec<[PseudoJet<T>; 8]>,
) -> ClusteringResult<T> {
    let mut jets = SmallVec::<[Jet<T>; 8]>::new();
    if active.is_empty() {
        return ClusteringResult { jets };
    }

    // Runtime settings currently store clustering parameters as `f64`, so this is
    // the intentional configuration precision boundary for these two values.
    let min_pt = F::from_f64(min_jpt);

    if d_r <= 0.0 || !d_r.is_finite() {
        jets.extend(
            active
                .drain(..)
                .filter(|jet| jet.pt() >= min_pt)
                .map(PseudoJet::jet),
        );
        sort_jets(&mut jets);
        return ClusteringResult { jets };
    }

    let r = F::from_f64(d_r);
    let inv_r2 = (r.clone() * r).inv();

    while !active.is_empty() {
        let mut best_action = BestAction::Beam(0);
        let mut best_distance: Option<F<T>> = None;

        for i in 0..active.len() {
            update_best(
                &mut best_action,
                &mut best_distance,
                BestAction::Beam(i),
                beam_distance(&active[i], algorithm),
            );

            for j in (i + 1)..active.len() {
                update_best(
                    &mut best_action,
                    &mut best_distance,
                    BestAction::Pair(i, j),
                    pair_distance(&active[i], &active[j], &inv_r2, algorithm),
                );
            }
        }

        match best_action {
            BestAction::Beam(index) => {
                let jet = active.remove(index);
                if jet.pt() >= min_pt {
                    jets.push(jet.jet());
                }
            }
            BestAction::Pair(i, j) => {
                let merged = PseudoJet::merged(&active[i], &active[j]);
                active[i] = merged;
                active.remove(j);
            }
        }
    }

    sort_jets(&mut jets);
    ClusteringResult { jets }
}

fn update_best<T: FloatLike>(
    best_action: &mut BestAction,
    best_distance: &mut Option<F<T>>,
    action: BestAction,
    candidate_distance: F<T>,
) {
    let should_update = best_distance
        .as_ref()
        .is_none_or(|current| candidate_distance < current.clone());
    if should_update {
        *best_action = action;
        *best_distance = Some(candidate_distance);
    }
}

fn beam_distance<T: FloatLike>(jet: &PseudoJet<T>, algorithm: JetAlgorithm) -> F<T> {
    scale_for_algorithm(jet.pt2(), algorithm)
}

fn pair_distance<T: FloatLike>(
    left: &PseudoJet<T>,
    right: &PseudoJet<T>,
    inv_r2: &F<T>,
    algorithm: JetAlgorithm,
) -> F<T> {
    let left_scale = scale_for_algorithm(left.pt2(), algorithm);
    let right_scale = scale_for_algorithm(right.pt2(), algorithm);
    let distance_scale = if left_scale <= right_scale {
        left_scale
    } else {
        right_scale
    };
    distance_scale * delta_r2(left, right) * inv_r2.clone()
}

fn scale_for_algorithm<T: FloatLike>(pt2: F<T>, algorithm: JetAlgorithm) -> F<T> {
    match algorithm {
        JetAlgorithm::Kt => pt2,
        JetAlgorithm::CambridgeAachen => pt2.one(),
        JetAlgorithm::AntiKt => {
            let tiny = pt2.epsilon();
            if pt2 < tiny {
                pt2.max_value()
            } else {
                pt2.inv()
            }
        }
    }
}

fn delta_r2<T: FloatLike>(left: &PseudoJet<T>, right: &PseudoJet<T>) -> F<T> {
    let delta_y = left.rapidity() - right.rapidity();
    let delta_phi = wrapped_delta_phi(left.phi(), right.phi());
    delta_y.square() + delta_phi.square()
}

fn wrapped_delta_phi<T: FloatLike>(lhs: F<T>, rhs: F<T>) -> F<T> {
    let pi = lhs.PI();
    let two_pi = lhs.TAU();
    let delta = (lhs - rhs).norm();
    if delta > pi { two_pi - delta } else { delta }
}

fn sort_jets<T: FloatLike>(jets: &mut SmallVec<[Jet<T>; 8]>) {
    jets.sort_unstable_by(compare_jets_descending_pt);
}

fn compare_jets_descending_pt<T: FloatLike>(lhs: &Jet<T>, rhs: &Jet<T>) -> Ordering {
    rhs.pt()
        .partial_cmp(&lhs.pt())
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            lhs.rapidity()
                .partial_cmp(&rhs.rapidity())
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| lhs.phi().partial_cmp(&rhs.phi()).unwrap_or(Ordering::Equal))
        .then_with(|| lhs.constituents.cmp(&rhs.constituents))
}
