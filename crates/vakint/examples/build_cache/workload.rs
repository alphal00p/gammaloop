//! Frozen common-family workload shared by preparation and measurement.
use std::{collections::BTreeSet, path::PathBuf};

use eyre::{Result, ensure};
use serde::{Deserialize, Serialize};

/// A completed measurement and a lowering-only qualification have different
/// authorities. Historical references without a tag retain their strict rules.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(super) enum OwnerProofAuthority {
    #[default]
    CompletedPublicEvaluation,
    QualifiedLowering,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(super) struct OwnerProofReference {
    #[serde(default)]
    pub(super) authority: OwnerProofAuthority,
    pub(super) path: PathBuf,
    pub(super) sha256: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Root {
    pub(super) id: String,
    pub(super) powers: Vec<i64>,
    pub(super) rationale: String,
    #[serde(default)]
    pub(super) public_input: Option<String>,
    /// Complete owner-mapped input union, authenticated by the pinned authority.
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
    #[serde(default)]
    pub(super) owner_proof: Option<OwnerProofReference>,
}

impl Workload {
    pub(super) fn validate_final_selection(&self) -> Result<()> {
        self.validate()?;
        ensure!(
            self.owner_proof.is_some(),
            "final selection requires a pinned owner proof"
        );
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
        if let Some(proof) = &self.owner_proof {
            ensure!(
                !proof.path.as_os_str().is_empty()
                    && proof.sha256.len() == 64
                    && proof.sha256.bytes().all(|byte| byte.is_ascii_hexdigit()),
                "invalid owner-proof reference"
            );
        }
        ensure!(!self.training.is_empty(), "empty training workload");
        ensure!(!self.tuning.is_empty(), "empty tuning suite");
        ensure!(
            self.held_out.len() >= 12,
            "freeze at least twelve held-out inputs before preparation"
        );
        validate_root_groups(&[
            ("training", &self.training),
            ("tuning", &self.tuning),
            ("held_out", &self.held_out),
        ])
    }
}

/// The input suite is independently frozen before a LIGHT training selection.
/// It deliberately is not a Workload: no training roots have been selected.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct FrozenPublicSuites {
    pub(super) schema: String,
    pub(super) frozen_utc: String,
    /// Historical topology-selection provenance; distinct from the authenticated
    /// lowering/data/program evidence in qualification_origin.
    pub(super) topology_source_sha256: String,
    pub(super) ordered_case_ids: Vec<String>,
    pub(super) qualification_origin: serde_json::Value,
    pub(super) owner_proof: OwnerProofReference,
    pub(super) tuning: Vec<Root>,
    pub(super) held_out: Vec<Root>,
    pub(super) scope: String,
}

impl FrozenPublicSuites {
    pub(super) fn validate(&self) -> Result<()> {
        ensure!(
            self.schema == "vakint.five-loop-frozen-public-suites.v1"
                && !self.frozen_utc.trim().is_empty()
                && !self.scope.trim().is_empty()
                && self.topology_source_sha256.len() == 64
                && self
                    .topology_source_sha256
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit())
                && self.qualification_origin.is_object()
                && self.owner_proof.authority == OwnerProofAuthority::QualifiedLowering,
            "invalid frozen public suite provenance"
        );
        ensure!(
            self.tuning.len() == 1 && self.held_out.len() == 12,
            "qualified suite requires exactly one tuning and twelve held-out inputs"
        );
        validate_root_groups(&[("tuning", &self.tuning), ("held_out", &self.held_out)])?;
        ensure!(
            self.ordered_case_ids
                == self
                    .tuning
                    .iter()
                    .chain(&self.held_out)
                    .map(|root| root.id.clone())
                    .collect::<Vec<_>>()
                && self
                    .tuning
                    .iter()
                    .chain(&self.held_out)
                    .all(|root| !root.lowered_roots.is_empty() && !root.protected_roots.is_empty()),
            "frozen suite order or complete root sets differ"
        );
        Ok(())
    }
}

/// Shared input invariants; keeping training optional here does not relax the
/// Workload requirement that an actual final training selection be nonempty.
fn validate_root_groups(groups: &[(&str, &Vec<Root>)]) -> Result<()> {
    let mut identities = BTreeSet::new();
    let mut training = BTreeSet::new();
    let mut tuning = BTreeSet::new();
    let mut held_out = BTreeSet::new();
    let mut protected = BTreeSet::new();
    for &(kind, roots) in groups {
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

#[cfg(test)]
pub(super) mod tests {
    use super::*;

    pub(crate) fn selected() -> Workload {
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
            owner_proof: Some(OwnerProofReference {
                authority: OwnerProofAuthority::CompletedPublicEvaluation,
                path: "owner-proof.json".into(),
                sha256: "c".repeat(64),
            }),
        }
    }

    #[test]
    fn proof_authorities_are_explicit_with_strict_legacy_default() {
        let reference = |authority: Option<&str>| {
            let mut value = serde_json::json!({"path":"proof.json","sha256":"a".repeat(64)});
            if let Some(authority) = authority {
                value["authority"] = authority.into();
            }
            serde_json::from_value::<OwnerProofReference>(value)
        };
        assert_eq!(
            reference(None).unwrap().authority,
            OwnerProofAuthority::CompletedPublicEvaluation
        );
        assert_eq!(
            reference(Some("qualified_lowering")).unwrap().authority,
            OwnerProofAuthority::QualifiedLowering
        );
        assert!(reference(Some("measurement_or_lowering_whatever_is_available")).is_err());
    }

    #[test]
    fn final_selection_requires_frozen_public_and_owner_roots() {
        let mut workload = selected();
        workload.validate_final_selection().unwrap();
        workload.held_out[0].lowered_roots.clear();
        workload.validate().unwrap();
        assert!(workload.validate_final_selection().is_err());
        workload = selected();
        workload.owner_proof = None;
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
