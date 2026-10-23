//! Frozen common-family workload shared by preparation and measurement.
use std::collections::BTreeSet;

use eyre::{Result, ensure};
use serde::Deserialize;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Root {
    pub(super) id: String,
    pub(super) powers: Vec<i64>,
    pub(super) rationale: String,
    #[serde(default)]
    pub(super) public_input: Option<String>,
    /// Complete owner-mapped input union, authenticated by public evaluation.
    /// Empty is permitted while constructing an explicitly provisional pilot.
    #[serde(default)]
    pub(super) protected_roots: Vec<Vec<i64>>,
    /// Complete common-family input union before native owner routing.
    #[serde(default)]
    pub(super) lowered_roots: Vec<Vec<i64>>,
}

/// Selection is input, never inferred from held-out timing results. All powers
/// use the authenticated application's common family and its ordered axes.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Workload {
    pub(super) schema: String,
    pub(super) data_manifest_sha256: String,
    pub(super) source_census_sha256: String,
    pub(super) training: Vec<Root>,
    pub(super) tuning: Vec<Root>,
    pub(super) held_out: Vec<Root>,
}

impl Workload {
    pub(super) fn validate_final_selection(&self) -> Result<()> {
        self.validate()?;
        ensure!(
            self.tuning
                .iter()
                .chain(&self.held_out)
                .all(|root| { !root.lowered_roots.is_empty() && !root.protected_roots.is_empty() }),
            "final selection requires complete frozen public-lowering and protected owner-root sets"
        );
        Ok(())
    }

    pub(super) fn validate(&self) -> Result<()> {
        ensure!(
            self.schema == "vakint.five-loop-cache-workload.v1",
            "unknown workload schema"
        );
        for hash in [&self.data_manifest_sha256, &self.source_census_sha256] {
            ensure!(
                hash.len() == 64 && hash.bytes().all(|b| b.is_ascii_hexdigit()),
                "invalid workload SHA256"
            );
        }
        ensure!(!self.training.is_empty(), "empty training workload");
        ensure!(!self.tuning.is_empty(), "empty tuning suite");
        ensure!(
            self.held_out.len() >= 12,
            "freeze at least twelve held-out inputs before preparation"
        );
        let mut identities = BTreeSet::new();
        let mut training = BTreeSet::new();
        let mut tuning = BTreeSet::new();
        let mut held_out = BTreeSet::new();
        let mut protected = BTreeSet::new();
        for (kind, roots) in [
            ("training", &self.training),
            ("tuning", &self.tuning),
            ("held_out", &self.held_out),
        ] {
            let mut keys = BTreeSet::new();
            for root in roots {
                ensure!(
                    !root.id.is_empty() && identities.insert(&root.id),
                    "duplicate or empty root identifier {:?}",
                    root.id
                );
                ensure!(
                    !root.rationale.trim().is_empty(),
                    "missing selection rationale for {}",
                    root.id
                );
                ensure!(
                    kind == "training"
                        || root
                            .public_input
                            .as_ref()
                            .is_some_and(|input| !input.trim().is_empty()),
                    "missing public input for {kind} root {}",
                    root.id
                );
                ensure!(
                    root.powers.len() == 15,
                    "{} needs fifteen ordered family powers",
                    root.id
                );
                ensure!(
                    keys.insert(&root.powers),
                    "duplicate {kind} root {}",
                    root.id
                );
                let mut root_protected = BTreeSet::new();
                for powers in &root.protected_roots {
                    ensure!(
                        powers.len() == 15 && root_protected.insert(powers),
                        "invalid or duplicate protected root for {}",
                        root.id
                    );
                    if kind != "training" {
                        protected.insert(powers);
                    }
                }
                let mut lowered = BTreeSet::new();
                ensure!(
                    root.lowered_roots
                        .iter()
                        .all(|powers| powers.len() == 15 && lowered.insert(powers)),
                    "invalid or duplicate lowered root for {}",
                    root.id
                );
                match kind {
                    "training" => {
                        training.insert(&root.powers);
                    }
                    "held_out" => {
                        held_out.insert(&root.powers);
                    }
                    _ => {
                        tuning.insert(&root.powers);
                    }
                }
            }
        }
        ensure!(
            training.is_disjoint(&held_out),
            "a held-out input was selected for training"
        );
        ensure!(
            tuning.is_disjoint(&held_out),
            "a held-out input occurs in the tuning suite"
        );
        ensure!(
            training.is_disjoint(&tuning),
            "a tuning input was selected for training"
        );
        ensure!(
            training.is_disjoint(&protected),
            "a protected tuning or held-out owner root was selected for training"
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn selected() -> Workload {
        let root = |id: String, power: i64| Root {
            id,
            powers: [vec![power], vec![1; 14]].concat(),
            rationale: "independently selected closure variation".into(),
            public_input: Some(format!("vakint::I5L({power})")),
            protected_roots: vec![[vec![power], vec![1; 14]].concat()],
            lowered_roots: vec![[vec![power], vec![1; 14]].concat()],
        };
        Workload {
            schema: "vakint.five-loop-cache-workload.v1".into(),
            data_manifest_sha256: "a".repeat(64),
            source_census_sha256: "b".repeat(64),
            training: vec![root("training".into(), 1)],
            tuning: vec![root("tuning".into(), 2)],
            held_out: (0..12)
                .map(|i| root(format!("held-out-{i}"), i + 3))
                .collect(),
        }
    }

    #[test]
    fn final_selection_requires_frozen_public_and_owner_roots() {
        let mut workload = selected();
        workload.validate_final_selection().unwrap();
        workload.held_out[0].lowered_roots.clear();
        workload.validate().unwrap();
        assert!(workload.validate_final_selection().is_err());
        workload = selected();
        workload.tuning[0].protected_roots.clear();
        assert!(workload.validate_final_selection().is_err());
    }

    #[test]
    fn frozen_held_out_inputs_cannot_be_training_or_tuning_roots() {
        let mut workload = selected();
        workload.validate().unwrap();
        workload.held_out[0].powers = workload.training[0].powers.clone();
        assert!(
            workload
                .validate()
                .unwrap_err()
                .to_string()
                .contains("selected for training")
        );
        workload = selected();
        workload.held_out[0].powers = workload.tuning[0].powers.clone();
        assert!(
            workload
                .validate()
                .unwrap_err()
                .to_string()
                .contains("tuning suite")
        );
    }

    #[test]
    fn selection_requires_twelve_distinct_held_out_inputs_and_provenance() {
        let mut workload = selected();
        workload.held_out.pop();
        assert!(workload.validate().is_err());
        workload = selected();
        workload.source_census_sha256 = "unrecorded".into();
        assert!(workload.validate().is_err());
        workload = selected();
        workload.training[0].rationale.clear();
        assert!(workload.validate().is_err());
    }

    #[test]
    fn tuning_is_nonempty_and_disjoint_from_training() {
        let mut workload = selected();
        workload.tuning.clear();
        assert!(
            workload
                .validate()
                .unwrap_err()
                .to_string()
                .contains("empty tuning")
        );
        workload = selected();
        workload.tuning[0].powers = workload.training[0].powers.clone();
        assert!(
            workload
                .validate()
                .unwrap_err()
                .to_string()
                .contains("tuning input was selected")
        );
    }

    #[test]
    fn measurement_suites_require_public_inputs_but_training_does_not() {
        let mut workload = selected();
        workload.training[0].public_input = None;
        workload.validate().unwrap();
        workload.tuning[0].public_input = None;
        assert!(
            workload
                .validate()
                .unwrap_err()
                .to_string()
                .contains("missing public input for tuning")
        );
        workload = selected();
        workload.held_out[0].public_input = Some(" \n ".into());
        assert!(
            workload
                .validate()
                .unwrap_err()
                .to_string()
                .contains("missing public input for held_out")
        );
    }

    #[test]
    fn protected_public_owner_roots_cannot_be_selected_for_training() {
        for held_out in [false, true] {
            let mut workload = selected();
            let training = workload.training[0].powers.clone();
            let root = if held_out {
                &mut workload.held_out[0]
            } else {
                &mut workload.tuning[0]
            };
            root.protected_roots.push(training);
            assert!(
                workload
                    .validate()
                    .unwrap_err()
                    .to_string()
                    .contains("protected")
            );
        }
    }
}
