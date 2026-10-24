//! Pinned native owner maps from complete successful public evaluations only.
use super::{
    digest,
    workload::{Root, Workload},
};
use eyre::{Result, ensure};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

pub(super) struct OwnerProof {
    document: Value,
    identity: Value,
    correlations: Value,
}

fn string<'a>(value: &'a Value, field: &str) -> Result<&'a str> {
    value[field]
        .as_str()
        .ok_or_else(|| eyre::eyre!("missing owner-proof {field}"))
}
fn sha256<'a>(value: &'a Value, field: &str) -> Result<&'a str> {
    let hash = string(value, field)?;
    ensure!(
        hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit()),
        "invalid or missing owner-proof SHA256 {field}"
    );
    Ok(hash)
}
fn available_executable(document: &Value, base: &Path, path: &str, hash: &str) -> Result<()> {
    let expected = sha256(document, hash)?;
    if let Some(path) = document[path]
        .as_str()
        .filter(|path| !path.trim().is_empty())
    {
        let path = resolved(base, path);
        if path.exists() {
            ensure!(
                digest(&path)? == expected,
                "available executable differs from provenance"
            );
        }
    }
    Ok(())
}
fn compiled_bindings() -> [(&'static str, Option<&'static str>); 5] {
    [
        (
            "main_revision",
            option_env!("VAKINT_CACHE_BUILD_GIT_REVISION"),
        ),
        (
            "native_revision",
            option_env!("VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
        ),
        (
            "kernel_revision",
            option_env!("VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
        ),
        (
            "main_source_snapshot_sha256",
            option_env!("VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"),
        ),
        (
            "native_source_snapshot_sha256",
            option_env!("VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256"),
        ),
    ]
}
fn cohort(document: &Value, binding: &Value, expected: &[(&str, Option<&str>); 5]) -> Result<()> {
    let compiled = &document["compiled_cohort"];
    ensure!(
        compiled.is_object(),
        "missing final compiled-cohort receipt"
    );
    for field in ["main_revision", "native_revision", "kernel_revision"] {
        let revision = string(compiled, field)?;
        ensure!(
            revision.len() == 40 && revision.bytes().all(|byte| byte.is_ascii_hexdigit()),
            "missing actual compiled revision {field}"
        );
    }
    for field in [
        "main_source_snapshot_sha256",
        "native_source_snapshot_sha256",
    ] {
        sha256(compiled, field)?;
    }
    for field in [
        "mapper_source_sha256",
        "mapper_executable_sha256",
        "measurement_executable_sha256",
        "data_manifest_sha256",
        "loop_manifest_sha256",
        "application_sha256",
    ] {
        ensure!(
            sha256(compiled, field)? == sha256(document, field)?,
            "compiled-cohort binding differs: {field}"
        );
    }
    sha256(document, "cohort_receipt_sha256")?;
    ensure!(
        binding.is_object()
            && binding["artifact_digest"]
                .as_array()
                .is_some_and(|bytes| bytes.len() == 32
                    && bytes
                        .iter()
                        .all(|byte| byte.as_u64().is_some_and(|byte| byte <= 255)))
            && !string(binding, "family_fingerprint")?.trim().is_empty()
            && compiled["semantic_binding"] == *binding,
        "missing or mismatched native program/family binding"
    );
    for &(field, actual) in expected {
        ensure!(
            actual.is_some_and(
                |actual| !actual.trim().is_empty() && compiled[field].as_str() == Some(actual)
            ),
            "owner proof differs from actual compiled consumer: {field}"
        );
    }
    Ok(())
}
pub(super) fn validate_completed_outputs(
    receipt: &Value,
    directory: &Path,
    root: &Root,
) -> Result<()> {
    if receipt["success"] != true {
        // Failed attempts remain excluded without retry; they are never adopted.
        return Ok(());
    }
    ensure!(
        receipt["worker_success"] == true
            && receipt["fresh_process"] == true
            && receipt["supervision"]["success"] == true
            && receipt["supervision"]["exit_code"] == 0
            && receipt["supervision"]["deadline_exhausted"] == false,
        "incomplete successful worker receipt"
    );
    for (name, field) in [
        ("value.symbolica", "value_sha256"),
        ("guards.json", "guards_sha256"),
    ] {
        let expected = receipt[field]
            .as_str()
            .filter(|hash| hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit()))
            .ok_or_else(|| eyre::eyre!("missing completed-output SHA256 {field}"))?;
        ensure!(
            digest(&directory.join(name))? == expected,
            "completed output changed: {name}"
        );
    }
    sha256(receipt, "input_sha256")?;
    if receipt["identity"]["owner_proof"].is_object() {
        let pinned = &receipt["identity"]["owner_proof"];
        ensure!(
            receipt["native_keys_authenticated"] == true
                && receipt["actual_lowered_roots"] == json!(root.lowered_roots)
                && receipt["owner_proof"] == *pinned
                && receipt["cache_statistics"]["binding"]
                    == pinned["compiled_cohort"]["semantic_binding"],
            "resumed public key/program evidence differs from the pinned case"
        );
    }
    Ok(())
}

fn array<'a>(value: &'a Value, field: &str) -> Result<&'a Vec<Value>> {
    value[field]
        .as_array()
        .ok_or_else(|| eyre::eyre!("missing owner-proof {field}"))
}
fn roots(value: &Value) -> Result<BTreeSet<Vec<i64>>> {
    let keys: Vec<Vec<i64>> = serde_json::from_value(value.clone())?;
    let set = keys.iter().cloned().collect::<BTreeSet<_>>();
    ensure!(
        keys.iter().all(|key| key.len() == 15) && keys == set.iter().cloned().collect::<Vec<_>>(),
        "owner roots must be sorted unique fifteen-axis keys"
    );
    Ok(set)
}
fn resolved(base: &Path, path: &str) -> PathBuf {
    let path = Path::new(path);
    if path.is_absolute() {
        path.into()
    } else {
        base.join(path)
    }
}
fn evidence(document: &Value, base: &Path, path: &str, hash: &str) -> Result<PathBuf> {
    let path = resolved(base, string(document, path)?);
    ensure!(
        digest(&path)? == sha256(document, hash)?,
        "owner-proof evidence changed: {}",
        path.display()
    );
    Ok(path)
}

impl OwnerProof {
    fn read(workload: &Workload, workload_path: &Path) -> Result<Option<(Value, PathBuf, Value)>> {
        let Some(reference) = &workload.owner_proof else {
            return Ok(None);
        };
        let path = resolved(
            workload_path.parent().unwrap_or(Path::new(".")),
            &reference.path.to_string_lossy(),
        );
        let bytes = fs::read(&path)?;
        ensure!(
            format!("{:x}", Sha256::digest(&bytes)) == reference.sha256,
            "pinned owner-proof hash differs"
        );
        let document: Value = serde_json::from_slice(&bytes)?;
        ensure!(
            document["schema"] == "vakint.public-lowering-owner-map-proof.v1",
            "unknown completed owner-proof schema"
        );
        let identity = json!({"path":path,"sha256":reference.sha256,"mapper_executable_sha256":document["mapper_executable_sha256"],"mapper_source_sha256":document["mapper_source_sha256"],"cohort_receipt_sha256":document["cohort_receipt_sha256"],"compiled_cohort":document["compiled_cohort"]});
        Ok(Some((document, path, identity)))
    }

    pub(super) fn pinned_identity(
        workload: &Workload,
        workload_path: &Path,
        artifacts: &Path,
    ) -> Result<Value> {
        Self::pinned_identity_with_bindings(
            workload,
            workload_path,
            artifacts,
            &compiled_bindings(),
        )
    }

    fn pinned_identity_with_bindings(
        workload: &Workload,
        workload_path: &Path,
        artifacts: &Path,
        expected: &[(&str, Option<&str>); 5],
    ) -> Result<Value> {
        let Some((document, path, identity)) = Self::read(workload, workload_path)? else {
            return Ok(Value::Null);
        };
        Self::validate_data(&document, artifacts)?;
        let binding = document["semantic_binding"].clone();
        let proof = Self::validate(document, workload, &binding, identity, expected)?;
        proof.validate_evidence(path.parent().unwrap_or(Path::new(".")), workload)?;
        Ok(proof.identity)
    }

    pub(super) fn load(
        workload: &Workload,
        workload_path: &Path,
        artifacts: &Path,
        binding: &Value,
    ) -> Result<Option<Self>> {
        let Some((document, path, identity)) = Self::read(workload, workload_path)? else {
            return Ok(None);
        };
        let base = path.parent().unwrap_or(Path::new("."));
        Self::validate_data(&document, artifacts)?;
        let proof = Self::validate(document, workload, binding, identity, &compiled_bindings())?;
        proof.validate_evidence(base, workload)?;
        ensure!(
            digest(&path)? == workload.owner_proof.as_ref().unwrap().sha256,
            "owner proof changed during validation"
        );
        Ok(Some(proof))
    }

    fn validate_data(document: &Value, artifacts: &Path) -> Result<()> {
        for (relative, field) in [
            ("manifest.json", "data_manifest_sha256"),
            ("5L/manifest.json", "loop_manifest_sha256"),
            ("5L/ibps/application.json", "application_sha256"),
        ] {
            ensure!(
                digest(&artifacts.join(relative))? == sha256(document, field)?,
                "owner-proof data or application differs"
            );
        }
        Ok(())
    }

    fn validate_evidence(&self, base: &Path, workload: &Workload) -> Result<()> {
        let mut frozen = Vec::new();
        for (path, hash) in [
            ("frozen_suite_path", "frozen_suite_sha256"),
            ("reserve_path", "reserve_sha256"),
            ("protocol_path", "protocol_sha256"),
        ] {
            let path = evidence(&self.document, base, path, hash)?;
            frozen.push(serde_json::from_slice::<Value>(&fs::read(path)?)?);
        }
        let cohort_path = evidence(
            &self.document,
            base,
            "cohort_receipt_path",
            "cohort_receipt_sha256",
        )?;
        let receipt: Value = serde_json::from_slice(&fs::read(cohort_path)?)?;
        ensure!(
            receipt == self.document["compiled_cohort"],
            "compiled-cohort receipt differs"
        );
        // Executable hashes are mandatory provenance, but archived ELFs are not.
        // If a binary is available, still authenticate its bytes.
        available_executable(
            &self.document,
            base,
            "mapper_executable_path",
            "mapper_executable_sha256",
        )?;
        available_executable(
            &self.document,
            base,
            "measurement_executable_path",
            "measurement_executable_sha256",
        )?;
        for root in workload.tuning.iter().chain(&workload.held_out) {
            let row = self.row(&root.id)?;
            let workload_path = evidence(row, base, "workload_path", "workload_sha256")?;
            evidence(row, base, "source_census_path", "source_census_sha256")?;
            let source_workload: Value = serde_json::from_slice(&fs::read(workload_path)?)?;
            let role = string(row, "logical_cohort")?;
            let source_cases = array(
                &source_workload,
                if role == "tuning" {
                    "tuning"
                } else {
                    "held_out"
                },
            )?;
            let declarations = source_cases
                .iter()
                .filter(|case| case["id"] == root.id)
                .collect::<Vec<_>>();
            let fixed_cases = if role == "tuning" {
                array(&frozen[0], "tuning")?.iter().collect::<Vec<_>>()
            } else {
                array(&frozen[0], "held_out")?
                    .iter()
                    .chain(array(&frozen[1], "held_out_reserves")?)
                    .collect()
            };
            let fixed = fixed_cases
                .into_iter()
                .filter(|case| case["id"] == root.id)
                .collect::<Vec<_>>();
            ensure!(
                declarations.len() == 1
                    && fixed.len() == 1
                    && declarations[0]["public_input"] == json!(root.public_input)
                    && declarations[0]["powers"] == json!(root.powers)
                    && fixed[0]["public_input"] == declarations[0]["public_input"]
                    && fixed[0]["powers"] == declarations[0]["powers"],
                "historical source identity or role differs from the pinned case"
            );
            let receipt_path = evidence(row, base, "receipt_path", "receipt_sha256")?;
            let receipt: Value = serde_json::from_slice(&fs::read(&receipt_path)?)?;
            ensure!(
                receipt["success"] == true
                    && receipt["worker_success"] == true
                    && receipt["fresh_process"] == true
                    && receipt["variant"] == "uncached"
                    && receipt["sample"] == 0
                    && receipt["supervision"]["success"] == true
                    && receipt["supervision"]["exit_code"] == 0
                    && receipt["supervision"]["deadline_exhausted"] == false,
                "owner proof requires a completed successful public report"
            );
            ensure!(
                receipt["case_id"] == root.id
                    && sha256(&receipt, "input_sha256")? == sha256(row, "input_sha256")?
                    && receipt["actual_lowered_roots"] == row["actual_lowered_roots"]
                    && receipt["cache_statistics"]["binding"] == self.document["semantic_binding"]
                    && receipt["identity"]["data_manifest_sha256"] == workload.data_manifest_sha256
                    && sha256(&receipt["identity"], "workload_sha256")?
                        == sha256(row, "workload_sha256")?
                    && sha256(&receipt["identity"], "source_census_sha256")?
                        == sha256(row, "source_census_sha256")?
                    && sha256(&receipt["identity"], "executable_sha256")?
                        == sha256(&self.document, "measurement_executable_sha256")?
                    && receipt["identity"]["cache_sha256"].is_null(),
                "owner-proof public report differs"
            );
            ensure!(
                digest(&receipt_path.with_file_name("value.symbolica"))?
                    == string(&receipt, "value_sha256")?
                    && digest(&receipt_path.with_file_name("guards.json"))?
                        == string(&receipt, "guards_sha256")?,
                "completed public result or guards changed"
            );
        }
        Ok(())
    }

    fn validate(
        document: Value,
        workload: &Workload,
        binding: &Value,
        identity: Value,
        expected: &[(&str, Option<&str>); 5],
    ) -> Result<Self> {
        workload.validate_final_selection()?;
        cohort(&document, binding, expected)?;
        ensure!(
            document["semantic_binding"] == *binding
                && binding.is_object()
                && binding["physical_arity"] == 15
                && binding["layer"] == "saved-application-route-apply-v1"
                && document["data_manifest_sha256"] == workload.data_manifest_sha256,
            "owner proof differs from the loaded application"
        );
        let mut indexed = BTreeMap::new();
        for row in array(&document, "rows")? {
            ensure!(
                indexed.insert(string(row, "id")?, row).is_none(),
                "duplicate owner-proof case"
            );
        }
        let mut signatures = BTreeSet::new();
        let mut selected = Vec::new();
        for (cohort, cases) in [
            ("tuning", &workload.tuning),
            ("held-out", &workload.held_out),
        ] {
            for root in cases {
                let row = indexed
                    .get(root.id.as_str())
                    .ok_or_else(|| eyre::eyre!("owner proof omits {}", root.id))?;
                let source = format!(
                    "{:x}",
                    Sha256::digest(
                        root.public_input
                            .as_ref()
                            .ok_or_else(|| eyre::eyre!("missing public input"))?
                            .as_bytes()
                    )
                );
                ensure!(
                    row["logical_cohort"] == cohort
                        && row["public_input_source_sha256"] == source
                        && row["reduction_succeeded"] == true
                        && row["complete_owner_transport"] == true
                        && row["mapping_error"].is_null()
                        && row["application_binding"] == *binding,
                    "incomplete or mismatched owner case {}",
                    root.id
                );
                let raw = roots(&row["actual_lowered_roots"])?;
                let protected = roots(&row["protected_owner_roots"])?;
                ensure!(
                    !raw.is_empty()
                        && !protected.is_empty()
                        && raw == roots(&json!(root.lowered_roots))?
                        && protected == roots(&json!(root.protected_roots))?,
                    "frozen owner sets differ"
                );
                let mut covered = BTreeSet::new();
                let mut endpoints = BTreeSet::new();
                for transport in array(row, "native_transports")? {
                    let key: Vec<i64> = serde_json::from_value(transport["powers"].clone())?;
                    ensure!(covered.insert(key), "duplicate native source transport");
                    let mut local = BTreeSet::new();
                    for endpoint in array(transport, "owner_endpoints")? {
                        let key: Vec<i64> = serde_json::from_value(endpoint["powers"].clone())?;
                        ensure!(
                            key.len() == 15
                                && local.insert(key.clone())
                                && !string(endpoint, "coefficient")?.trim().is_empty()
                                && string(endpoint, "coefficient")?.trim() != "0"
                                && endpoint["recorded_coverage_admitted"].is_boolean(),
                            "invalid native endpoint"
                        );
                        endpoints.insert(key);
                    }
                    ensure!(
                        array(transport, "nonzero_conditions")?
                            .iter()
                            .all(|guard| guard
                                .as_str()
                                .is_some_and(
                                    |guard| !guard.trim().is_empty() && guard.trim() != "0"
                                )),
                        "invalid transport guard"
                    );
                }
                ensure!(
                    covered == raw && endpoints == protected,
                    "native transports do not cover the exact frozen sets"
                );
                let signature = format!(
                    "{:x}",
                    Sha256::digest(serde_json::to_vec(&row["protected_owner_roots"])?)
                );
                ensure!(
                    row["protected_root_set_sha256"] == signature,
                    "owner-set signature differs"
                );
                if cohort == "held-out" {
                    signatures.insert(signature);
                }
                ensure!(
                    workload
                        .training
                        .iter()
                        .all(|root| !protected.contains(&root.powers)),
                    "training contains a protected owner root"
                );
                selected.push((root.id.clone(), cohort, protected));
            }
        }
        ensure!(
            signatures.len() >= 12,
            "require twelve distinct nonempty held-out owner-root sets"
        );
        let mut correlations = Vec::new();
        for (index, (a, cohort_a, left)) in selected.iter().enumerate() {
            for (b, cohort_b, right) in selected.iter().skip(index + 1) {
                let shared = left.intersection(right).collect::<Vec<_>>();
                if !shared.is_empty() {
                    correlations.push(json!({"case_a":a,"cohort_a":cohort_a,"case_b":b,"cohort_b":cohort_b,"shared_owner_roots":shared,"identical_owner_sets":left==right}));
                }
            }
        }
        Ok(Self {
            document,
            identity,
            correlations: json!(correlations),
        })
    }
    fn row(&self, id: &str) -> Result<&Value> {
        array(&self.document, "rows")?
            .iter()
            .find(|row| row["id"] == id)
            .ok_or_else(|| eyre::eyre!("missing owner case"))
    }
    pub(super) fn observe(
        &self,
        id: &str,
        raw: &BTreeSet<Vec<i64>>,
        protected: &[Vec<i64>],
    ) -> Result<()> {
        let row = self.row(id)?;
        ensure!(
            raw == &roots(&row["actual_lowered_roots"])?
                && roots(&json!(protected))? == roots(&row["protected_owner_roots"])?,
            "observed public roots differ from the pinned owner proof"
        );
        Ok(())
    }
    pub(super) fn identity(&self) -> &Value {
        &self.identity
    }
    pub(super) fn correlations(&self) -> &Value {
        &self.correlations
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    fn expected_bindings() -> [(&'static str, Option<&'static str>); 5] {
        [
            (
                "main_revision",
                Some("1111111111111111111111111111111111111111"),
            ),
            (
                "native_revision",
                Some("2222222222222222222222222222222222222222"),
            ),
            (
                "kernel_revision",
                Some("3333333333333333333333333333333333333333"),
            ),
            (
                "main_source_snapshot_sha256",
                Some("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
            ),
            (
                "native_source_snapshot_sha256",
                Some("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
            ),
        ]
    }
    fn cohort(document: &Value, binding: &Value) -> Result<()> {
        super::cohort(document, binding, &expected_bindings())
    }
    fn validate(
        document: Value,
        workload: &Workload,
        binding: &Value,
        identity: Value,
    ) -> Result<OwnerProof> {
        OwnerProof::validate(document, workload, binding, identity, &expected_bindings())
    }
    fn pinned_identity(workload: &Workload, path: &Path, artifacts: &Path) -> Result<Value> {
        OwnerProof::pinned_identity_with_bindings(workload, path, artifacts, &expected_bindings())
    }
    #[test]
    fn absent_or_empty_current_bindings_cannot_accept_a_valid_proof() {
        let workload = super::super::workload::tests::selected();
        let (document, binding) = document(&workload);
        for index in 0..5 {
            for absent in [None, Some(""), Some(" \n\t")] {
                let mut expected = expected_bindings();
                expected[index].1 = absent;
                assert!(super::cohort(&document, &binding, &expected).is_err());
            }
        }
    }

    fn document(workload: &Workload) -> (Value, Value) {
        let binding = json!({"artifact_digest":vec![1;32],"family_fingerprint":"fixture","physical_arity":15,"layer":"saved-application-route-apply-v1"});
        let rows=[("tuning",&workload.tuning),("held-out",&workload.held_out)].into_iter().flat_map(|(cohort,cases)|cases.iter().map(move |root|(cohort,root))).map(|(cohort,root)|{
            json!({"id":root.id,"logical_cohort":cohort,"public_input_source_sha256":format!("{:x}",Sha256::digest(root.public_input.as_ref().unwrap().as_bytes())),
                "reduction_succeeded":true,"complete_owner_transport":true,"mapping_error":null,"application_binding":binding,
                "actual_lowered_roots":root.lowered_roots,"protected_owner_roots":root.protected_roots,
                "protected_root_set_sha256":format!("{:x}",Sha256::digest(serde_json::to_vec(&root.protected_roots).unwrap())),
                "native_transports":root.lowered_roots.iter().map(|key|json!({"powers":key,"owner_endpoints":root.protected_roots.iter().map(|key|json!({"powers":key,"coefficient":"1","recorded_coverage_admitted":true})).collect::<Vec<_>>(),"nonzero_conditions":[]})).collect::<Vec<_>>()})
        }).collect::<Vec<_>>();
        (
            {
                let hash = "a".repeat(64);
                let compiled = json!({"main_revision":"1".repeat(40),"native_revision":"2".repeat(40),"kernel_revision":"3".repeat(40),
                    "main_source_snapshot_sha256":hash,"native_source_snapshot_sha256":hash,
                    "mapper_source_sha256":hash,"mapper_executable_sha256":hash,"measurement_executable_sha256":hash,
                    "data_manifest_sha256":workload.data_manifest_sha256,"loop_manifest_sha256":hash,"application_sha256":hash,"semantic_binding":binding});
                json!({"schema":"vakint.public-lowering-owner-map-proof.v1","semantic_binding":binding,
                    "data_manifest_sha256":workload.data_manifest_sha256,"loop_manifest_sha256":hash,"application_sha256":hash,
                    "mapper_source_sha256":hash,"mapper_executable_sha256":hash,"measurement_executable_sha256":hash,
                    "cohort_receipt_path":"cohort.json","cohort_receipt_sha256":hash,"compiled_cohort":compiled,"rows":rows})
            },
            binding,
        )
    }
    pub(crate) fn fixture() -> OwnerProof {
        let workload = super::super::workload::tests::selected();
        let (document, binding) = document(&workload);
        validate(
            document,
            &workload,
            &binding,
            json!({"sha256":"synthetic fixture"}),
        )
        .unwrap()
    }
    #[test]
    fn missing_null_or_mismatched_cohort_bindings_never_authenticate() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        cohort(&original, &binding).unwrap();
        for field in [
            "mapper_source_sha256",
            "mapper_executable_sha256",
            "measurement_executable_sha256",
            "data_manifest_sha256",
            "loop_manifest_sha256",
            "application_sha256",
        ] {
            for value in [Value::Null, json!(""), json!("wrong")] {
                let mut changed = original.clone();
                changed[field] = value.clone();
                changed["compiled_cohort"][field] = value;
                assert!(cohort(&changed, &binding).is_err(), "{field}");
            }
        }
        for field in [
            "main_revision",
            "native_revision",
            "kernel_revision",
            "main_source_snapshot_sha256",
            "native_source_snapshot_sha256",
        ] {
            let mut changed = original.clone();
            changed["compiled_cohort"][field] = Value::Null;
            assert!(cohort(&changed, &binding).is_err(), "{field}");
        }
        let mut changed = original;
        changed["compiled_cohort"]["semantic_binding"] = Value::Null;
        assert!(cohort(&changed, &binding).is_err());
    }

    #[test]
    fn self_consistent_other_compiled_cohort_is_rejected() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        for field in ["main_revision", "native_revision", "kernel_revision"] {
            let mut changed = original.clone();
            changed["compiled_cohort"][field] = json!("f".repeat(40));
            assert!(cohort(&changed, &binding).is_err());
        }
        for field in [
            "main_source_snapshot_sha256",
            "native_source_snapshot_sha256",
        ] {
            let mut changed = original.clone();
            changed["compiled_cohort"][field] = json!("f".repeat(64));
            assert!(cohort(&changed, &binding).is_err());
        }
    }

    #[test]
    fn source_only_provenance_does_not_require_archived_elf_files() {
        let document = json!({"binary_path":"certainly-absent-owner-proof-test-binary",
            "binary_sha256":"a".repeat(64)});
        available_executable(&document, Path::new("."), "binary_path", "binary_sha256").unwrap();
        let mut missing = document;
        missing["binary_sha256"] = Value::Null;
        assert!(
            available_executable(&missing, Path::new("."), "binary_path", "binary_sha256").is_err()
        );
    }

    #[test]
    fn resumed_success_cannot_omit_or_change_output_evidence() {
        let workload = super::super::workload::tests::selected();
        let root = &workload.held_out[0];
        let missing = json!({"success":true,"worker_success":true,"fresh_process":true,
            "supervision":{"success":true,"exit_code":0,"deadline_exhausted":false}});
        assert!(validate_completed_outputs(&missing, Path::new("."), root).is_err());
        let failed = json!({"success":false,"worker_success":false});
        validate_completed_outputs(&failed, Path::new("."), root).unwrap();
    }

    #[test]
    fn completed_output_sidecars_are_rechecked_before_resume() {
        let workload = super::super::workload::tests::selected();
        let root = &workload.held_out[0];
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("owner-proof-output-{}-{stamp}", std::process::id()));
        fs::create_dir(&path).unwrap();
        fs::write(path.join("value.symbolica"), b"complete").unwrap();
        fs::write(path.join("guards.json"), b"[]").unwrap();
        let receipt = json!({"success":true,"worker_success":true,"fresh_process":true,
            "input_sha256":"a".repeat(64),"supervision":{"success":true,"exit_code":0,"deadline_exhausted":false},
            "value_sha256":digest(&path.join("value.symbolica")).unwrap(),"guards_sha256":digest(&path.join("guards.json")).unwrap()});
        validate_completed_outputs(&receipt, &path, root).unwrap();
        fs::write(path.join("value.symbolica"), b"changed").unwrap();
        assert!(validate_completed_outputs(&receipt, &path, root).is_err());
        fs::remove_dir_all(path).unwrap();
    }

    #[test]
    fn pinned_resume_rechecks_receipt_chain_without_archived_binaries() {
        let mut workload = super::super::workload::tests::selected();
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let base =
            std::env::temp_dir().join(format!("owner-proof-chain-{}-{stamp}", std::process::id()));
        let artifacts = base.join("data");
        fs::create_dir_all(artifacts.join("5L/ibps")).unwrap();
        for path in [
            "manifest.json",
            "5L/manifest.json",
            "5L/ibps/application.json",
        ] {
            fs::write(artifacts.join(path), path.as_bytes()).unwrap();
        }
        workload.data_manifest_sha256 = digest(&artifacts.join("manifest.json")).unwrap();
        let (mut document, binding) = document(&workload);
        for (path, field) in [
            ("5L/manifest.json", "loop_manifest_sha256"),
            ("5L/ibps/application.json", "application_sha256"),
        ] {
            let hash = digest(&artifacts.join(path)).unwrap();
            document[field] = json!(hash);
            document["compiled_cohort"][field] = json!(hash);
        }
        fs::write(
            base.join("cohort.json"),
            serde_json::to_vec(&document["compiled_cohort"]).unwrap(),
        )
        .unwrap();
        document["cohort_receipt_sha256"] = json!(digest(&base.join("cohort.json")).unwrap());
        let declared = |root: &Root| json!({"id":root.id,"powers":root.powers,"public_input":root.public_input});
        let suite = json!({"tuning":workload.tuning.iter().map(declared).collect::<Vec<_>>(),
            "held_out":workload.held_out.iter().map(declared).collect::<Vec<_>>()});
        for (name, path_field, hash_field, value) in [
            (
                "suite.json",
                "frozen_suite_path",
                "frozen_suite_sha256",
                suite,
            ),
            (
                "reserves.json",
                "reserve_path",
                "reserve_sha256",
                json!({"held_out_reserves":[]}),
            ),
            (
                "protocol.json",
                "protocol_path",
                "protocol_sha256",
                json!({"fixture":true}),
            ),
        ] {
            fs::write(base.join(name), serde_json::to_vec(&value).unwrap()).unwrap();
            document[path_field] = json!(name);
            document[hash_field] = json!(digest(&base.join(name)).unwrap());
        }
        fs::write(base.join("source-census.json"), b"source-only fixture").unwrap();
        document["mapper_executable_path"] = json!("absent-mapper-ELF");
        document["measurement_executable_path"] = json!("absent-measurement-ELF");
        for (index, row) in document["rows"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .enumerate()
        {
            let directory = base.join(format!("case-{index}"));
            fs::create_dir(&directory).unwrap();
            fs::write(directory.join("value.symbolica"), b"complete").unwrap();
            fs::write(directory.join("guards.json"), b"[]").unwrap();
            let root = workload
                .tuning
                .iter()
                .chain(&workload.held_out)
                .find(|root| row["id"] == root.id)
                .unwrap();
            let source_workload = json!({"tuning":if row["logical_cohort"]=="tuning" {vec![declared(root)]}else{vec![]},
                "held_out":if row["logical_cohort"]=="held-out" {vec![declared(root)]}else{vec![]}});
            fs::write(
                directory.join("workload.json"),
                serde_json::to_vec(&source_workload).unwrap(),
            )
            .unwrap();
            row["workload_path"] = json!(format!("case-{index}/workload.json"));
            row["workload_sha256"] = json!(digest(&directory.join("workload.json")).unwrap());
            row["source_census_path"] = json!("source-census.json");
            row["source_census_sha256"] = json!(digest(&base.join("source-census.json")).unwrap());
            let receipt = json!({"case_id":row["id"],"variant":"uncached","sample":0,
                "success":true,"worker_success":true,"fresh_process":true,"input_sha256":"b".repeat(64),
                "actual_lowered_roots":row["actual_lowered_roots"],"cache_statistics":{"binding":binding},
                "identity":{"data_manifest_sha256":workload.data_manifest_sha256,"executable_sha256":"a".repeat(64),"cache_sha256":null,
                    "workload_sha256":row["workload_sha256"],"source_census_sha256":row["source_census_sha256"]},
                "supervision":{"success":true,"exit_code":0,"deadline_exhausted":false},
                "value_sha256":digest(&directory.join("value.symbolica")).unwrap(),"guards_sha256":digest(&directory.join("guards.json")).unwrap()});
            fs::write(
                directory.join("measurement.json"),
                serde_json::to_vec(&receipt).unwrap(),
            )
            .unwrap();
            row["input_sha256"] = receipt["input_sha256"].clone();
            row["receipt_path"] = json!(format!("case-{index}/measurement.json"));
            row["receipt_sha256"] = json!(digest(&directory.join("measurement.json")).unwrap());
        }
        let proof_path = base.join("proof.json");
        fs::write(&proof_path, serde_json::to_vec(&document).unwrap()).unwrap();
        let reference = workload.owner_proof.as_mut().unwrap();
        reference.path = proof_path.clone();
        reference.sha256 = digest(&proof_path).unwrap();
        pinned_identity(&workload, &base.join("workload.json"), &artifacts).unwrap();
        fs::write(
            base.join("case-0/value.symbolica"),
            b"changed after successful receipt",
        )
        .unwrap();
        assert!(pinned_identity(&workload, &base.join("workload.json"), &artifacts).is_err());
        fs::write(base.join("case-0/value.symbolica"), b"complete").unwrap();
        fs::write(base.join("source-census.json"), b"changed source evidence").unwrap();
        assert!(pinned_identity(&workload, &base.join("workload.json"), &artifacts).is_err());
        fs::write(base.join("source-census.json"), b"source-only fixture").unwrap();
        let source_path = base.join("case-0/workload.json");
        let mut wrong_role: Value =
            serde_json::from_slice(&fs::read(&source_path).unwrap()).unwrap();
        wrong_role["held_out"] = wrong_role["tuning"].clone();
        wrong_role["tuning"] = json!([]);
        fs::write(&source_path, serde_json::to_vec(&wrong_role).unwrap()).unwrap();
        let source_hash = digest(&source_path).unwrap();
        document["rows"][0]["workload_sha256"] = json!(source_hash);
        let receipt_path = base.join("case-0/measurement.json");
        let mut wrong_receipt: Value =
            serde_json::from_slice(&fs::read(&receipt_path).unwrap()).unwrap();
        wrong_receipt["identity"]["workload_sha256"] = json!(source_hash);
        fs::write(&receipt_path, serde_json::to_vec(&wrong_receipt).unwrap()).unwrap();
        document["rows"][0]["receipt_sha256"] = json!(digest(&receipt_path).unwrap());
        fs::write(&proof_path, serde_json::to_vec(&document).unwrap()).unwrap();
        workload.owner_proof.as_mut().unwrap().sha256 = digest(&proof_path).unwrap();
        assert!(pinned_identity(&workload, &base.join("workload.json"), &artifacts).is_err());
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn requires_complete_matching_sources_and_loaded_binding() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        assert!(validate(original.clone(), &workload, &binding, Value::Null).is_ok());
        for (field, value) in [
            ("reduction_succeeded", json!(false)),
            ("complete_owner_transport", json!(false)),
            ("public_input_source_sha256", json!("wrong")),
        ] {
            let mut changed = original.clone();
            changed["rows"][0][field] = value;
            assert!(validate(changed, &workload, &binding, Value::Null).is_err());
        }
        let mut wrong = binding.clone();
        wrong["artifact_digest"][0] = json!(2);
        assert!(validate(original, &workload, &wrong, Value::Null).is_err());
    }
    #[test]
    fn rejects_missing_duplicate_or_inconsistent_native_transports() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        let mut missing = original.clone();
        missing["rows"][0]["native_transports"] = json!([]);
        assert!(validate(missing, &workload, &binding, Value::Null).is_err());
        let mut duplicate = original.clone();
        let row = duplicate["rows"][0].clone();
        duplicate["rows"].as_array_mut().unwrap().push(row);
        assert!(validate(duplicate, &workload, &binding, Value::Null).is_err());
        let mut changed = original;
        changed["rows"][0]["native_transports"][0]["owner_endpoints"][0]["powers"][0] = json!(100);
        assert!(validate(changed, &workload, &binding, Value::Null).is_err());
    }
    #[test]
    fn aliases_do_not_count_as_twelve_and_correlations_are_reported() {
        let mut workload = super::super::workload::tests::selected();
        workload.held_out[11].protected_roots = workload.held_out[0].protected_roots.clone();
        let (aliased, binding) = document(&workload);
        assert!(validate(aliased, &workload, &binding, Value::Null).is_err());
        let mut workload = super::super::workload::tests::selected();
        workload.tuning[0].protected_roots = workload.held_out[0].protected_roots.clone();
        let (correlated, binding) = document(&workload);
        let proof = validate(correlated, &workload, &binding, Value::Null).unwrap();
        assert_eq!(proof.correlations().as_array().unwrap().len(), 1);
    }
}
