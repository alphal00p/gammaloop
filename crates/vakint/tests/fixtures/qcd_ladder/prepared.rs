//! Verified scalar-input persistence; the saved preparation is not current-run work.

use super::{
    LadderFixture, PREPARED_INPUT_MANIFEST, PreparationTimings, Projection, Result, Vakint,
};
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, io::Read, path::Path};
use symbolica::atom::{Atom, AtomCore};

const INPUT_FILE: &str = "vakint-input.symbolica";
const CONDITIONS_FILE: &str = "projector-conditions.symbolica";
const SCHEMA: &str = "vakint.qcd-ladder.prepared-input.v2";
// Bump this definition when changing physical conventions. The source binding
// also covers native generation/assembly in mod.rs, excluding persistence IO.
const DEFINITION: &str = "massless-qcd-rung-family-native-uv-v2";

fn five_loops() -> usize {
    5
}

fn is_five_loops(value: &usize) -> bool {
    *value == 5
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreparedFile {
    pub bytes: u64,
    pub sha256: String,
}

impl PreparedFile {
    fn from_bytes(bytes: &[u8]) -> Self {
        Self {
            bytes: bytes.len() as u64,
            sha256: format!("{:x}", Sha256::digest(bytes)),
        }
    }

    fn from_file(path: &Path) -> std::io::Result<Self> {
        let mut file = std::fs::File::open(path)?;
        let mut digest = Sha256::new();
        let mut buffer = [0_u8; 65_536];
        let mut bytes = 0_u64;
        loop {
            let count = file.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            digest.update(&buffer[..count]);
            bytes += count as u64;
        }
        Ok(Self {
            bytes,
            sha256: format!("{:x}", digest.finalize()),
        })
    }

    fn is_present_regular(directory: &Path, filename: &str) -> Result<bool> {
        match std::fs::symlink_metadata(directory.join(filename)) {
            Ok(metadata) if metadata.is_file() => Ok(true),
            Ok(_) => Err(format!(
                "prepared input {filename} must be a regular file, not a symlink or directory"
            )
            .into()),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(false),
            Err(error) => Err(error.into()),
        }
    }

    fn read(&self, directory: &Path, filename: &str) -> Result<String> {
        if !Self::is_present_regular(directory, filename)? {
            return Err(format!("prepared input {filename} is missing").into());
        }
        let path = directory.join(filename);
        if std::fs::metadata(&path)?.len() != self.bytes {
            return Err(format!("prepared input {filename} has an unexpected byte length").into());
        }
        let bytes = std::fs::read(path)?;
        if Self::from_bytes(&bytes).sha256 != self.sha256 {
            return Err(format!("prepared input {filename} failed SHA-256 validation").into());
        }
        Ok(String::from_utf8(bytes)?)
    }
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreparedInputProducer {
    /// All producer sources, including its IO/assembly module, kept as provenance.
    pub fixture_sources_sha256: BTreeMap<String, String>,
    pub source_commit: Option<String>,
    pub executable_sha256: Option<String>,
}

/// Self-contained identity and provenance of one complete scalar fixture.
///
/// Scientific sources match this loop family, or the exact delivered four-rung
/// v1 identity. `mod.rs` also records producer provenance; persistence/reporting
/// changes do not alter the physical definition.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreparedInputManifest {
    pub schema: String,
    pub definition: String,
    #[serde(default = "five_loops", skip_serializing_if = "is_five_loops")]
    pub loop_count: usize,
    pub projection: Projection,
    pub scientific_sources_sha256: BTreeMap<String, String>,
    pub scientific_assembly_sha256: String,
    pub producer: PreparedInputProducer,
    pub input: PreparedFile,
    pub nonzero_conditions_file: PreparedFile,
    pub nonzero_conditions: Vec<String>,
    pub contributions: usize,
    pub scalar_terms: usize,
    pub preparation_seconds: f64,
    pub timings: PreparationTimings,
}

impl PreparedInputManifest {
    /// Exact identity of the delivered four-rung definition, before its purely
    /// loop-count parameterization. Never accept arbitrary stale source hashes.
    fn legacy_scientific_sources() -> BTreeMap<String, String> {
        [
            (
                "physical.rs",
                "754b9a5b2b0b9701fc666c29ce6728229d78ace9cb10948033c17d9e42275436",
            ),
            (
                "projector.rs",
                "bbd9a9a8ce08e0c29a632c74586f7b51dc04cc2d1ed8802a6ef6a43bdbf5a53d",
            ),
            (
                "vacuum.rs",
                "05cdaa106c6f27fa16ea7c473f4a83ac0c4e7100eba9cc58ee4358e95aedb905",
            ),
            (
                "ladder.dot",
                "990490ae32441fcbc279390b97aedf433d226a7b4986faea1eaf87e15ca79ece",
            ),
            (
                "qcd.json",
                "229179d56fec310801f9db3d5697d6f3947b6a3c112d9b856ba0fe6fc1fe6a3b",
            ),
            (
                "phi3.json",
                "5e7b76223b63752eac1200744ee945be2bc9b125890b1430271faea9917b1dd3",
            ),
        ]
        .into_iter()
        .map(|(name, hash)| (name.to_owned(), hash.to_owned()))
        .collect()
    }

    fn scientific_sources() -> BTreeMap<String, String> {
        [
            ("physical.rs", include_bytes!("physical.rs").as_slice()),
            ("projector.rs", include_bytes!("projector.rs").as_slice()),
            ("vacuum.rs", include_bytes!("vacuum.rs").as_slice()),
            ("ladder.dot", include_bytes!("ladder.dot").as_slice()),
            ("qcd.json", include_bytes!("../qcd.json").as_slice()),
            ("phi3.json", include_bytes!("../phi3.json").as_slice()),
        ]
        .into_iter()
        .map(|(name, bytes)| (name.to_owned(), PreparedFile::from_bytes(bytes).sha256))
        .collect()
    }

    fn scientific_assembly() -> String {
        let source = include_str!("mod.rs");
        for marker in [
            "use std::{",
            "use symbolica::{",
            "impl Projection {",
            "type Shape =",
            "    std::fs::write(evidence.join(\"physical-power-counting.txt\")",
        ] {
            assert_eq!(
                source.matches(marker).count(),
                1,
                "fixture scientific source marker must be unique"
            );
        }
        let imports = source
            .split_once("use std::{")
            .expect("fixture import source marker")
            .1;
        let (imports_before_symbolica, symbolica) = imports
            .split_once("use symbolica::{")
            .expect("fixture symbolic import marker");
        let symbolica = symbolica
            .split_once("\n};")
            .expect("fixture symbolic import ending")
            .0;
        let projection = source
            .split_once("impl Projection {")
            .expect("fixture projection source marker")
            .1
            .split_once("\n}\n")
            .expect("fixture projection source ending")
            .0;
        let assembly = source
            .split_once("type Shape =")
            .expect("fixture assembly source marker")
            .1
            .split_once("    std::fs::write(evidence.join(\"physical-power-counting.txt\")")
            .expect("fixture assembly source ending")
            .0;
        PreparedFile::from_bytes(
            format!("imports\n{imports_before_symbolica}use symbolica::{{{symbolica}\nprojection\n{projection}\nassembly\n{assembly}").as_bytes(),
        )
        .sha256
    }

    fn validate(&self, expected_projection: Projection, expected_loop_count: usize) -> Result<()> {
        let legacy = self.schema == "vakint.qcd-ladder.prepared-input.v1"
            && self.definition == "massless-qcd-four-rung-native-uv-v1"
            && self.loop_count == 5;
        if !legacy && (self.schema != SCHEMA || self.definition != DEFINITION) {
            return Err("unrecognized QCD ladder prepared-input schema or definition".into());
        }
        if self.loop_count == 0 || self.loop_count != expected_loop_count {
            return Err("prepared input does not use the requested QCD ladder loop count".into());
        }
        if self.projection != expected_projection {
            return Err("prepared input does not use the requested QCD ladder projection".into());
        }
        let sources = if legacy {
            Self::legacy_scientific_sources()
        } else {
            Self::scientific_sources()
        };
        let assembly = if legacy {
            "16be3edc83d1cfad237e7762269bc02d2335e21173b01ee881b5c5834b4aa7e9".to_owned()
        } else {
            Self::scientific_assembly()
        };
        if self.scientific_sources_sha256 != sources || self.scientific_assembly_sha256 != assembly
        {
            return Err(
                "prepared input scientific source definition does not match this fixture".into(),
            );
        }
        if self.producer.fixture_sources_sha256.len() != 7
            || self
                .scientific_sources_sha256
                .iter()
                .any(|(name, hash)| self.producer.fixture_sources_sha256.get(name) != Some(hash))
            || !self.producer.fixture_sources_sha256.contains_key("mod.rs")
        {
            return Err("prepared input producer source inventory is inconsistent".into());
        }
        if self.contributions == 0
            || (self.projection == Projection::ZeroMomentumMetric && self.contributions != 1)
            || (self.loop_count == 5
                && self.projection == Projection::QuadraticTransverse
                && self.contributions != 21)
            || self.scalar_terms < self.contributions
            || self.input.bytes == 0
        {
            return Err("prepared input is missing complete QCD ladder contributions".into());
        }
        if [
            self.preparation_seconds,
            self.timings.graph_generation_seconds,
            self.timings.scalar_contraction_seconds,
            self.timings.taylor_seconds,
            self.timings.angular_seconds,
            self.timings.ingress_seconds,
            self.timings.total_seconds,
        ]
        .into_iter()
        .any(|value| !value.is_finite() || value < 0.0)
        {
            return Err("prepared input has invalid historical preparation timings".into());
        }
        let mut conditions = self.nonzero_conditions.clone();
        conditions.sort();
        let mut expected = match self.projection {
            Projection::QuadraticTransverse => vec![
                "1".to_owned(),
                "rustred::{}::d".to_owned(),
                "-1+rustred::{}::d".to_owned(),
                "2+rustred::{}::d".to_owned(),
                "-2*rustred::{}::d+rustred::{}::d^2+rustred::{}::d^3".to_owned(),
            ],
            Projection::ZeroMomentumMetric => vec![],
        };
        expected.sort();
        if conditions != expected {
            return Err(
                "prepared input does not retain the complete projector nonzero conditions".into(),
            );
        }
        Ok(())
    }
}

impl LadderFixture {
    /// An unfinished bundle must not trigger silent expensive regeneration.
    pub fn has_prepared_input(directory: &Path) -> Result<bool> {
        if PreparedFile::is_present_regular(directory, PREPARED_INPUT_MANIFEST)? {
            return Ok(true);
        }
        for filename in [INPUT_FILE, CONDITIONS_FILE] {
            if PreparedFile::is_present_regular(directory, filename)? {
                return Err("prepared scalar payloads exist without their manifest; choose a new evidence directory or restore the manifest".into());
            }
        }
        Ok(false)
    }

    /// Reuse a complete, checked scalar input without graph, trace or Taylor work.
    /// Corrupt or stale data is an error; callers must never silently regenerate it.
    pub fn load_prepared(directory: &Path, expected_projection: Projection) -> Result<Self> {
        Self::load_prepared_at_loop_count(directory, expected_projection, 5)
    }

    pub fn load_prepared_at_loop_count(
        directory: &Path,
        expected_projection: Projection,
        expected_loop_count: usize,
    ) -> Result<Self> {
        if !PreparedFile::is_present_regular(directory, PREPARED_INPUT_MANIFEST)? {
            return Err("QCD ladder prepared-input manifest is missing".into());
        }
        let manifest: PreparedInputManifest =
            serde_json::from_slice(&std::fs::read(directory.join(PREPARED_INPUT_MANIFEST))?)?;
        manifest.validate(expected_projection, expected_loop_count)?;
        let input_text = manifest.input.read(directory, INPUT_FILE)?;
        let conditions_text = manifest
            .nonzero_conditions_file
            .read(directory, CONDITIONS_FILE)?;
        if conditions_text != manifest.nonzero_conditions.join("\n") {
            return Err("prepared input condition payload disagrees with its manifest".into());
        }
        Vakint::initialize_vakint_symbols();
        let input = Atom::parse(&input_text, "qcd_ladder", Default::default())?;
        if input.to_canonical_string() != input_text {
            return Err(
                "prepared scalar input did not preserve its canonical symbolic identity".into(),
            );
        }
        let nonzero_conditions = manifest
            .nonzero_conditions
            .iter()
            .map(|condition| Atom::parse(condition, "qcd_ladder", Default::default()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        if nonzero_conditions
            .iter()
            .zip(&manifest.nonzero_conditions)
            .any(|(condition, text)| condition.to_canonical_string() != *text)
        {
            return Err("prepared input condition lost its canonical symbolic identity".into());
        }
        Ok(Self {
            loop_count: manifest.loop_count,
            input,
            preparation_seconds: manifest.preparation_seconds,
            scalar_terms: manifest.scalar_terms,
            contributions: manifest.contributions,
            nonzero_conditions,
            projection: manifest.projection,
            timings: manifest.timings.clone(),
            prepared_input: Some(manifest),
        })
    }

    /// Write a reusable scalar input, committing its manifest after its payloads.
    /// A different existing bundle is never overwritten implicitly.
    pub fn save_prepared(&self, directory: &Path) -> Result<PreparedInputManifest> {
        let input = self.input.to_canonical_string();
        let conditions = self
            .nonzero_conditions
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>();
        let conditions_text = conditions.join("\n");
        let manifest = if let Some(manifest) = &self.prepared_input {
            manifest.clone()
        } else {
            let mut sources = PreparedInputManifest::scientific_sources();
            sources.insert(
                "mod.rs".to_owned(),
                PreparedFile::from_bytes(include_bytes!("mod.rs")).sha256,
            );
            PreparedInputManifest {
                schema: SCHEMA.to_owned(),
                definition: DEFINITION.to_owned(),
                loop_count: self.loop_count,
                projection: self.projection,
                scientific_sources_sha256: PreparedInputManifest::scientific_sources(),
                scientific_assembly_sha256: PreparedInputManifest::scientific_assembly(),
                producer: PreparedInputProducer {
                    fixture_sources_sha256: sources,
                    source_commit: None,
                    // A debug executable can contain gigabytes of DWARF.
                    // Its supervising test/build receipt owns that identity;
                    // scientific source hashes still bind every fresh bundle.
                    executable_sha256: if cfg!(debug_assertions) {
                        None
                    } else {
                        Some(PreparedFile::from_file(&std::env::current_exe()?)?.sha256)
                    },
                },
                input: PreparedFile::from_bytes(input.as_bytes()),
                nonzero_conditions_file: PreparedFile::from_bytes(conditions_text.as_bytes()),
                nonzero_conditions: conditions,
                contributions: self.contributions,
                scalar_terms: self.scalar_terms,
                preparation_seconds: self.preparation_seconds,
                timings: self.timings.clone(),
            }
        };
        manifest.validate(self.projection, self.loop_count)?;
        if manifest.input.sha256 != PreparedFile::from_bytes(input.as_bytes()).sha256
            || manifest.nonzero_conditions_file.sha256
                != PreparedFile::from_bytes(conditions_text.as_bytes()).sha256
        {
            return Err("prepared fixture was modified after loading".into());
        }
        let manifest_path = directory.join(PREPARED_INPUT_MANIFEST);
        if PreparedFile::is_present_regular(directory, PREPARED_INPUT_MANIFEST)? {
            let previous =
                Self::load_prepared_at_loop_count(directory, self.projection, self.loop_count)?;
            if previous.input != self.input
                || previous.nonzero_conditions != self.nonzero_conditions
            {
                return Err(
                    "a different prepared input already exists; choose a new evidence directory"
                        .into(),
                );
            }
            return Ok(previous
                .prepared_input
                .expect("loaded fixture has provenance"));
        }
        for (filename, contents) in [(INPUT_FILE, &input), (CONDITIONS_FILE, &conditions_text)] {
            let path = directory.join(filename);
            if PreparedFile::is_present_regular(directory, filename)?
                && std::fs::read(&path)? != contents.as_bytes()
            {
                return Err(format!("different prepared payload {filename} already exists; choose a new evidence directory").into());
            }
        }
        std::fs::create_dir_all(directory)?;
        for (filename, contents) in [(INPUT_FILE, &input), (CONDITIONS_FILE, &conditions_text)] {
            if !PreparedFile::is_present_regular(directory, filename)? {
                let mut file = std::fs::OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(directory.join(filename))?;
                std::io::Write::write_all(&mut file, contents.as_bytes())?;
                file.sync_all()?;
            }
        }
        let mut staged = tempfile::NamedTempFile::new_in(directory)?;
        std::io::Write::write_all(&mut staged, &serde_json::to_vec_pretty(&manifest)?)?;
        staged.as_file().sync_all()?;
        staged.persist_noclobber(manifest_path)?;
        Ok(manifest)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> LadderFixture {
        Vakint::initialize_vakint_symbols();
        LadderFixture {
            loop_count: 5,
            input: Atom::parse(
                "1+qcd_ladder_control::{}::x",
                "qcd_ladder",
                Default::default(),
            )
            .unwrap(),
            preparation_seconds: 1.0,
            scalar_terms: 21,
            contributions: 21,
            nonzero_conditions: [
                "1",
                "rustred::{}::d",
                "-1+rustred::{}::d",
                "2+rustred::{}::d",
                "-2*rustred::{}::d+rustred::{}::d^2+rustred::{}::d^3",
            ]
            .into_iter()
            .map(|text| Atom::parse(text, "qcd_ladder", Default::default()).unwrap())
            .collect(),
            projection: Projection::QuadraticTransverse,
            timings: PreparationTimings::default(),
            prepared_input: None,
        }
    }

    fn edit_manifest(directory: &Path, edit: impl FnOnce(&mut serde_json::Value)) {
        let path = directory.join(PREPARED_INPUT_MANIFEST);
        let mut value = serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        edit(&mut value);
        std::fs::write(path, serde_json::to_vec_pretty(&value).unwrap()).unwrap();
    }

    #[test]
    fn prepared_scalar_partial_payloads_do_not_trigger_regeneration_or_overwrite() {
        let directory = tempfile::tempdir().unwrap();
        assert!(!LadderFixture::has_prepared_input(directory.path()).unwrap());
        std::fs::write(directory.path().join(INPUT_FILE), "partial").unwrap();
        assert!(LadderFixture::has_prepared_input(directory.path()).is_err());
        assert!(fixture().save_prepared(directory.path()).is_err());
        assert_eq!(
            std::fs::read_to_string(directory.path().join(INPUT_FILE)).unwrap(),
            "partial"
        );
        assert!(!directory.path().join(PREPARED_INPUT_MANIFEST).exists());
    }

    #[cfg(unix)]
    #[test]
    fn prepared_scalar_rejects_symlink_payload_without_touching_target() {
        let directory = tempfile::tempdir().unwrap();
        let target = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(target.path(), "private evidence").unwrap();
        std::os::unix::fs::symlink(target.path(), directory.path().join(INPUT_FILE)).unwrap();
        assert!(fixture().save_prepared(directory.path()).is_err());
        assert!(LadderFixture::has_prepared_input(directory.path()).is_err());
        assert_eq!(
            std::fs::read_to_string(target.path()).unwrap(),
            "private evidence"
        );
    }

    #[test]
    fn prepared_scalar_roundtrip_retains_guards_and_historical_provenance() {
        let directory = tempfile::tempdir().unwrap();
        let original = fixture();
        original.save_prepared(directory.path()).unwrap();
        let loaded = LadderFixture::load_prepared(directory.path(), original.projection).unwrap();
        assert_eq!(loaded.input, original.input);
        assert_eq!(loaded.nonzero_conditions, original.nonzero_conditions);
        assert_eq!(loaded.contributions, 21);
        assert_eq!(loaded.preparation_seconds, 1.0);
        assert!(loaded.prepared_input.is_some());
        let second = tempfile::tempdir().unwrap();
        let saved = loaded.save_prepared(second.path()).unwrap();
        assert_eq!(
            saved.input.sha256,
            original
                .save_prepared(directory.path())
                .unwrap()
                .input
                .sha256
        );
    }

    #[test]
    fn prepared_scalar_loop_count_is_part_of_its_identity() {
        let directory = tempfile::tempdir().unwrap();
        let mut original = fixture();
        original.loop_count = 1;
        original.save_prepared(directory.path()).unwrap();
        assert!(LadderFixture::load_prepared(directory.path(), original.projection).is_err());
        let loaded =
            LadderFixture::load_prepared_at_loop_count(directory.path(), original.projection, 1)
                .unwrap();
        assert_eq!(loaded.loop_count, 1);
        assert_eq!(loaded.input, original.input);
    }

    #[test]
    fn delivered_five_loop_prepared_identity_remains_readable() {
        let directory = tempfile::tempdir().unwrap();
        let original = fixture();
        original.save_prepared(directory.path()).unwrap();
        edit_manifest(directory.path(), |manifest| {
            manifest["schema"] = "vakint.qcd-ladder.prepared-input.v1".into();
            manifest["definition"] = "massless-qcd-four-rung-native-uv-v1".into();
            manifest["scientific_sources_sha256"] =
                serde_json::to_value(PreparedInputManifest::legacy_scientific_sources()).unwrap();
            manifest["scientific_assembly_sha256"] =
                "16be3edc83d1cfad237e7762269bc02d2335e21173b01ee881b5c5834b4aa7e9".into();
            let producer = manifest["producer"]["fixture_sources_sha256"]["mod.rs"].clone();
            manifest["producer"]["fixture_sources_sha256"] =
                manifest["scientific_sources_sha256"].clone();
            manifest["producer"]["fixture_sources_sha256"]["mod.rs"] = producer;
        });
        let loaded = LadderFixture::load_prepared(directory.path(), original.projection).unwrap();
        assert_eq!(loaded.input, original.input);
        assert_eq!(loaded.loop_count, 5);
        assert!(
            LadderFixture::load_prepared_at_loop_count(directory.path(), original.projection, 4)
                .is_err()
        );
    }

    #[test]
    fn prepared_scalar_rejects_corruption_and_never_regenerates() {
        let directory = tempfile::tempdir().unwrap();
        fixture().save_prepared(directory.path()).unwrap();
        std::fs::write(directory.path().join(INPUT_FILE), "0").unwrap();
        let error = LadderFixture::load_prepared(directory.path(), Projection::QuadraticTransverse)
            .err()
            .unwrap()
            .to_string();
        assert!(error.contains("byte length"));
        assert_eq!(
            std::fs::read_to_string(directory.path().join(INPUT_FILE)).unwrap(),
            "0"
        );
    }

    #[test]
    fn prepared_scalar_rejects_changed_projection_and_scientific_source() {
        let directory = tempfile::tempdir().unwrap();
        fixture().save_prepared(directory.path()).unwrap();
        assert!(
            LadderFixture::load_prepared(directory.path(), Projection::ZeroMomentumMetric).is_err()
        );
        edit_manifest(directory.path(), |value| {
            value["scientific_sources_sha256"]["physical.rs"] = "stale".into()
        });
        let error = LadderFixture::load_prepared(directory.path(), Projection::QuadraticTransverse)
            .err()
            .unwrap()
            .to_string();
        assert!(error.contains("scientific source"));
    }

    #[test]
    fn prepared_scalar_rejects_changed_native_assembly_definition() {
        let directory = tempfile::tempdir().unwrap();
        fixture().save_prepared(directory.path()).unwrap();
        edit_manifest(directory.path(), |value| {
            value["scientific_assembly_sha256"] = "stale".into()
        });
        assert!(
            LadderFixture::load_prepared(directory.path(), Projection::QuadraticTransverse)
                .is_err()
        );
    }

    #[test]
    fn prepared_scalar_rejects_missing_guard_even_with_updated_file_checksum() {
        let directory = tempfile::tempdir().unwrap();
        fixture().save_prepared(directory.path()).unwrap();
        edit_manifest(directory.path(), |value| {
            value["nonzero_conditions"].as_array_mut().unwrap().pop();
            let text = value["nonzero_conditions"]
                .as_array()
                .unwrap()
                .iter()
                .map(|condition| condition.as_str().unwrap())
                .collect::<Vec<_>>()
                .join("\n");
            std::fs::write(directory.path().join(CONDITIONS_FILE), &text).unwrap();
            value["nonzero_conditions_file"] =
                serde_json::to_value(PreparedFile::from_bytes(text.as_bytes())).unwrap();
        });
        let error = LadderFixture::load_prepared(directory.path(), Projection::QuadraticTransverse)
            .err()
            .unwrap()
            .to_string();
        assert!(error.contains("complete projector"));
    }

    #[test]
    fn prepared_scalar_does_not_replace_a_different_existing_bundle() {
        let directory = tempfile::tempdir().unwrap();
        let original = fixture();
        original.save_prepared(directory.path()).unwrap();
        let before = std::fs::read(directory.path().join(PREPARED_INPUT_MANIFEST)).unwrap();
        let mut changed = fixture();
        changed.input = Atom::num(2);
        assert!(changed.save_prepared(directory.path()).is_err());
        assert_eq!(
            std::fs::read(directory.path().join(PREPARED_INPUT_MANIFEST)).unwrap(),
            before
        );
        assert_eq!(
            std::fs::read_to_string(directory.path().join(INPUT_FILE)).unwrap(),
            original.input.to_canonical_string()
        );
    }
}
