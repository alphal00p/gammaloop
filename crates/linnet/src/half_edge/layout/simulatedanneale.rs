use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::parser::GlobalData;

// ---------- Traits ----------
pub trait Energy<S> {
    /// Return the absolute energy of `next`.
    /// When a previous state/energy is provided, implementations may reuse or update it
    /// instead of recomputing from scratch.
    fn energy(&self, prev: Option<(&S, f64)>, next: &S) -> f64;

    /// Called once a proposal is accepted so implementations can clear or update caches.
    fn on_accept(&self, _state: &mut S) {}
}

pub trait Neighbor<S> {
    /// Propose a new state from `s` given step scale and temperature.
    fn propose(&self, s: &S, rng: &mut impl Rng, step: f64, temp: f64) -> S;
}

pub trait Schedule {
    /// Called at the start of each epoch. May mutate step & temp. Return `false` to stop.
    fn begin_epoch(
        &mut self,
        epoch: usize,
        last_accept_ratio: f64,
        step: &mut f64,
        temp: &mut f64,
    ) -> bool;
    /// Number of Metropolis steps per epoch.
    fn steps_per_epoch(&self) -> usize;
    /// Total max epochs budget (upper bound; you can early stop in begin_epoch).
    fn max_epochs(&self) -> usize;
}

pub struct SAConfig {
    pub seed: u64,
    pub temp: f64,
    pub step: f64,
}

pub struct SAStats {
    pub best_energy: f64,
    pub total_steps: usize,
    pub epochs_run: usize,
    pub last_accept_ratio: f64,
    pub last_temp: f64,
    pub last_step: f64,
}

pub fn anneal<S: Clone, N: Neighbor<S>, E: Energy<S>, Sch: Schedule, R: SeedableRng + Rng>(
    init: S,
    mut cfg: SAConfig,
    neigh: &N,
    energy: &E,
    sched: &mut Sch,
) -> (S, SAStats) {
    let mut rng = R::seed_from_u64(cfg.seed);
    let mut cur = init.clone();
    let mut current_energy = energy.energy(None, &cur);

    let mut best = cur.clone();
    let mut best_energy = current_energy;

    let mut total_steps = 0usize;
    let mut last_accept_ratio = 0.0;

    for ep in 0..sched.max_epochs() {
        if !sched.begin_epoch(ep, last_accept_ratio, &mut cfg.step, &mut cfg.temp) {
            return (
                best,
                SAStats {
                    best_energy,
                    total_steps,
                    epochs_run: ep,
                    last_accept_ratio,
                    last_temp: cfg.temp,
                    last_step: cfg.step,
                },
            );
        }

        let steps = sched.steps_per_epoch();
        let mut accepted = 0usize;

        for _ in 0..steps {
            total_steps += 1;
            let cand = neigh.propose(&cur, &mut rng, cfg.step, cfg.temp);
            let candidate_energy = energy.energy(Some((&cur, current_energy)), &cand);
            let delta = candidate_energy - current_energy;

            let accept = delta <= 0.0 || rng.gen::<f64>() < (-delta / cfg.temp).exp();
            if accept {
                let mut accepted_state = cand;
                energy.on_accept(&mut accepted_state);
                cur = accepted_state;
                current_energy = candidate_energy;
                accepted += 1;
                if candidate_energy < best_energy {
                    best_energy = candidate_energy;
                    best = cur.clone();
                }
            }
        }

        last_accept_ratio = accepted as f64 / steps as f64;
    }

    (
        best,
        SAStats {
            best_energy,
            total_steps,
            epochs_run: sched.max_epochs(),
            last_accept_ratio,
            last_temp: cfg.temp,
            last_step: cfg.step,
        },
    )
}

#[derive(Debug, Serialize, Deserialize, Clone, Copy)]
pub struct GeoSchedule {
    pub steps: usize,      // e.g. 200
    pub epochs: usize,     // e.g. 8
    pub cool: f64,         // e.g. 0.85
    pub accept_floor: f64, // e.g. 0.15
    pub step_shrink: f64,  // e.g. 0.7
    pub early_tol: f64,    // e.g. 1e-6 (relative best improvement)
    _prev_best: Option<f64>,
}

impl GeoSchedule {
    pub fn add_to_global(&self, global_data: &mut GlobalData) {
        global_data
            .statements
            .insert("steps".to_string(), self.steps.to_string());
        global_data
            .statements
            .insert("epochs".to_string(), self.epochs.to_string());
        global_data
            .statements
            .insert("cool".to_string(), self.cool.to_string());
        global_data
            .statements
            .insert("accept_floor".to_string(), self.accept_floor.to_string());
        global_data
            .statements
            .insert("step_shrink".to_string(), self.step_shrink.to_string());
        global_data
            .statements
            .insert("early_tol".to_string(), self.early_tol.to_string());
    }

    pub fn parse(global_data: &GlobalData) -> Self {
        let steps = global_data
            .statements
            .get("steps")
            .and_then(|s| s.parse().ok())
            .unwrap_or(200);
        let epochs = global_data
            .statements
            .get("epochs")
            .and_then(|s| s.parse().ok())
            .unwrap_or(8);
        let cool = global_data
            .statements
            .get("cool")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.85);
        let accept_floor = global_data
            .statements
            .get("accept_floor")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.15);
        let step_shrink = global_data
            .statements
            .get("step_shrink")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.7);
        let early_tol = global_data
            .statements
            .get("early_tol")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1e-6);
        Self::new(steps, epochs, cool, accept_floor, step_shrink, early_tol)
    }

    pub fn new(
        steps: usize,
        epochs: usize,
        cool: f64,
        accept_floor: f64,
        step_shrink: f64,
        early_tol: f64,
    ) -> Self {
        Self {
            steps,
            epochs,
            cool,
            accept_floor,
            step_shrink,
            early_tol,
            _prev_best: None,
        }
    }
}

impl Default for GeoSchedule {
    fn default() -> Self {
        Self::new(200, 8, 0.85, 0.15, 0.7, 1e-6)
    }
}

impl Schedule for GeoSchedule {
    fn begin_epoch(
        &mut self,
        epoch: usize,
        last_accept: f64,
        step: &mut f64,
        temp: &mut f64,
    ) -> bool {
        if epoch > 0 {
            *temp *= self.cool;
            if last_accept < self.accept_floor {
                *step *= self.step_shrink;
            }
            // Early-stop hook: caller can feed best energy back by resetting prev_best, or ignore.
        }
        epoch < self.epochs
    }
    fn steps_per_epoch(&self) -> usize {
        self.steps
    }
    fn max_epochs(&self) -> usize {
        self.epochs
    }
}
