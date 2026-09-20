//! Offline terminal values, persisted by RustRed's native Symbolica codec.
//!
//! This is a value catalog, not an IBP program or closure certificate. Embedded
//! native files have trusted generated provenance. Complete means producer-
//! declared finite coverage; the consumer still checks the exact terminal set.

use std::collections::{BTreeMap, BTreeSet};

use rustred::family::IntegralKey;
pub use rustred::persistence::TerminalCatalogCoverage as CatalogCoverage;
use rustred::persistence::{BinaryIoLimits, ExactTerminalCatalog};
use symbolica::atom::Atom;

/// Exact terminal values bound to one RustRed family identity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OfflineTerminalCatalog(ExactTerminalCatalog);

impl OfflineTerminalCatalog {
    pub fn from_terms(
        family_fingerprint: impl Into<String>,
        index_count: usize,
        terms: BTreeMap<IntegralKey, Atom>,
    ) -> Result<Self, String> {
        Self::from_terms_with_coverage(
            family_fingerprint,
            index_count,
            CatalogCoverage::Partial,
            terms,
        )
    }

    fn from_terms_with_coverage(
        family_fingerprint: impl Into<String>,
        index_count: usize,
        coverage: CatalogCoverage,
        terms: BTreeMap<IntegralKey, Atom>,
    ) -> Result<Self, String> {
        ExactTerminalCatalog::try_new(family_fingerprint, index_count, coverage, terms)
            .map(Self)
            .map_err(|error| error.to_string())
    }

    pub fn from_complete_terms(
        family_fingerprint: impl Into<String>,
        index_count: usize,
        terms: BTreeMap<IntegralKey, Atom>,
    ) -> Result<Self, String> {
        Self::from_terms_with_coverage(
            family_fingerprint,
            index_count,
            CatalogCoverage::Complete,
            terms,
        )
    }

    pub fn family_fingerprint(&self) -> &str {
        self.0.family_fingerprint()
    }
    pub fn index_count(&self) -> usize {
        self.0.index_count()
    }
    pub fn coverage(&self) -> CatalogCoverage {
        self.0.coverage()
    }
    pub fn is_complete(&self) -> bool {
        self.coverage() == CatalogCoverage::Complete
    }
    pub fn terms(&self) -> &BTreeMap<IntegralKey, Atom> {
        self.0.terms()
    }
    pub fn into_terms(self) -> BTreeMap<IntegralKey, Atom> {
        self.0.into_terms()
    }

    /// This flag does not prove coverage against a separate terminal owner.
    pub fn require_complete(&self) -> Result<(), String> {
        if self.is_complete() {
            Ok(())
        } else {
            Err(
                "offline terminal catalog is partial (it covers only finite probe residuals)"
                    .into(),
            )
        }
    }

    /// Resolve every declared key offline, never during ordinary evaluation.
    pub fn from_keys<I, F, E>(
        family_fingerprint: impl Into<String>,
        index_count: usize,
        keys: I,
        mut resolve: F,
    ) -> Result<Self, String>
    where
        I: IntoIterator<Item = IntegralKey>,
        F: FnMut(&IntegralKey) -> Result<Atom, E>,
        E: std::fmt::Display,
    {
        let mut terms = BTreeMap::new();
        for key in keys {
            let value = resolve(&key)
                .map_err(|error| format!("resolve terminal {:?}: {error}", key.powers()))?;
            if terms.insert(key.clone(), value).is_some() {
                return Err(format!("duplicate terminal key {:?}", key.powers()));
            }
        }
        Self::from_complete_terms(family_fingerprint, index_count, terms)
    }

    pub fn missing_keys<I>(&self, required: I) -> Vec<IntegralKey>
    where
        I: IntoIterator<Item = IntegralKey>,
    {
        required
            .into_iter()
            .filter(|key| !self.terms().contains_key(key))
            .collect()
    }

    pub fn require_keys<I>(&self, required: I) -> Result<(), String>
    where
        I: IntoIterator<Item = IntegralKey>,
    {
        let missing = self.missing_keys(required);
        if missing.is_empty() {
            return Ok(());
        }
        Err(format!(
            "offline terminal catalog is missing {} declared key(s), first {:?}",
            missing.len(),
            missing[0].powers()
        ))
    }

    /// Bind values to an installed output convention, with no unused raw keys.
    pub fn require_exact_keys(&self, required: &BTreeSet<IntegralKey>) -> Result<(), String> {
        self.require_keys(required.iter().cloned())?;
        if self.terms().len() != required.len() {
            return Err(
                "offline terminal catalog contains keys outside the installed output convention"
                    .into(),
            );
        }
        Ok(())
    }

    /// Native Atom/state persistence; no canonical strings or expression parse.
    pub fn encode(&self) -> Result<Vec<u8>, String> {
        self.0
            .encode_native(BinaryIoLimits::default())
            .map_err(|error| error.to_string())
    }

    /// Deterministic gzip around the unchanged native Atom/state catalog.
    pub fn encode_compressed(&self) -> Result<Vec<u8>, String> {
        super::compressed::encode(&self.encode()?)
    }

    /// Bound decompression before the native decoder and retain its trusted-
    /// generated provenance requirement. No text or legacy fallback is tried.
    pub fn decode_generated_compressed(
        input: &[u8],
        expected_family_fingerprint: &str,
        expected_index_count: usize,
        limits: BinaryIoLimits,
    ) -> Result<Self, String> {
        let bytes = super::compressed::decode(input, limits.max_program_bytes)?;
        ExactTerminalCatalog::decode_generated(
            &bytes,
            expected_family_fingerprint,
            expected_index_count,
            limits,
        )
        .map(Self)
        .map_err(|error| error.to_string())
    }

    /// Load trusted generated native data and bind its family and index arity.
    /// RustRed validates framing, exact values and metadata. Symbolica's native
    /// state reader is not a hostile-input allocation sandbox.
    pub fn decode_generated(
        input: &[u8],
        expected_family_fingerprint: &str,
        expected_index_count: usize,
    ) -> Result<Self, String> {
        ExactTerminalCatalog::decode_generated(
            input,
            expected_family_fingerprint,
            expected_index_count,
            BinaryIoLimits::default(),
        )
        .map(Self)
        .map_err(|error| error.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    fn sample() -> OfflineTerminalCatalog {
        OfflineTerminalCatalog::from_terms(
            "family-test-fingerprint",
            4,
            BTreeMap::from([(
                IntegralKey::try_new([1, 1, 0, 0]).unwrap(),
                vk_parse!("PR4+PR9/2").unwrap(),
            )]),
        )
        .unwrap()
    }

    #[test]
    fn encoding_is_deterministic_and_round_trips() {
        let catalog = sample();
        let encoded = catalog.encode().unwrap();
        assert_eq!(encoded, catalog.clone().encode().unwrap());
        let decoded =
            OfflineTerminalCatalog::decode_generated(&encoded, "family-test-fingerprint", 4)
                .unwrap();
        assert_eq!(decoded, catalog);
        assert_eq!(decoded.coverage(), CatalogCoverage::Partial);
        assert!(!decoded.is_complete());
        assert!(decoded.require_complete().is_err());
        assert!(
            OfflineTerminalCatalog::decode_generated(
                b"rustred-vakint-experimental-terminal-catalog\nschema=1\n",
                "family-test-fingerprint",
                4
            )
            .is_err()
        );
    }

    #[test]
    fn decode_rejects_wrong_family_arity_and_approximate_values() {
        let encoded = sample().encode().unwrap();
        assert!(
            OfflineTerminalCatalog::decode_generated(&encoded, "other", 4)
                .unwrap_err()
                .contains("fingerprint mismatch")
        );
        assert!(
            OfflineTerminalCatalog::decode_generated(&encoded, "family-test-fingerprint", 3)
                .unwrap_err()
                .contains("arity mismatch")
        );
        assert!(
            OfflineTerminalCatalog::from_terms(
                "family-test-fingerprint",
                4,
                BTreeMap::from([(
                    IntegralKey::try_new([1, 1, 0, 0]).unwrap(),
                    vk_parse!("1.25`40*PR4").unwrap()
                )])
            )
            .unwrap_err()
            .contains("approximate")
        );
    }

    #[test]
    fn from_keys_rejects_duplicate_keys() {
        let key = IntegralKey::try_new([1, 1, 0, 0]).unwrap();
        assert!(
            OfflineTerminalCatalog::from_keys(
                "family-test-fingerprint",
                4,
                [key.clone(), key],
                |_| Ok::<_, std::convert::Infallible>(Atom::num(1))
            )
            .unwrap_err()
            .contains("duplicate terminal key")
        );
    }

    #[test]
    fn compressed_catalog_preserves_exact_values_and_native_binding() {
        let catalog = sample();
        let bytes = catalog.encode_compressed().unwrap();
        assert_eq!(bytes, catalog.encode_compressed().unwrap());
        let decoded = OfflineTerminalCatalog::decode_generated_compressed(
            &bytes,
            catalog.family_fingerprint(),
            catalog.index_count(),
            BinaryIoLimits::default(),
        )
        .unwrap();
        assert_eq!(decoded, catalog);
        assert!(
            OfflineTerminalCatalog::decode_generated_compressed(
                &bytes,
                "foreign",
                catalog.index_count(),
                BinaryIoLimits::default()
            )
            .is_err()
        );
        assert!(
            OfflineTerminalCatalog::decode_generated_compressed(
                &catalog.encode().unwrap(),
                catalog.family_fingerprint(),
                catalog.index_count(),
                BinaryIoLimits::default()
            )
            .is_err()
        );
    }

    #[test]
    fn exact_output_coverage_rejects_missing_and_unused_raw_values() {
        let catalog = sample();
        let expected = catalog.terms().keys().cloned().collect::<BTreeSet<_>>();
        catalog.require_exact_keys(&expected).unwrap();
        assert!(
            catalog
                .require_exact_keys(&BTreeSet::new())
                .unwrap_err()
                .contains("outside")
        );
        let mut missing = expected;
        missing.insert(IntegralKey::try_new([2, 1, 0, 0]).unwrap());
        assert!(
            catalog
                .require_exact_keys(&missing)
                .unwrap_err()
                .contains("missing")
        );
    }

    #[test]
    fn from_keys_and_coverage_check_cover_declared_terminals() {
        let keys = vec![
            IntegralKey::try_new([2, 1, 0, 0]).unwrap(),
            IntegralKey::try_new([1, 1, 0, 0]).unwrap(),
        ];
        let catalog =
            OfflineTerminalCatalog::from_keys("family-test-fingerprint", 4, keys.clone(), |key| {
                Ok::<_, std::convert::Infallible>(
                    vk_parse!(format!("PR{}", key.powers()[0])).unwrap(),
                )
            })
            .unwrap();
        assert!(catalog.require_keys(keys.clone()).is_ok());
        assert_eq!(catalog.coverage(), CatalogCoverage::Complete);
        assert!(catalog.require_complete().is_ok());
        let absent = IntegralKey::try_new([3, 1, 0, 0]).unwrap();
        assert_eq!(catalog.missing_keys([absent.clone()]), vec![absent]);
        assert!(
            catalog
                .require_keys([IntegralKey::try_new([3, 1, 0, 0]).unwrap()])
                .is_err()
        );
        let decoded = OfflineTerminalCatalog::decode_generated(
            &catalog.encode().unwrap(),
            "family-test-fingerprint",
            4,
        )
        .unwrap();
        assert_eq!(decoded, catalog);
        assert_eq!(
            decoded.terms().keys().cloned().collect::<Vec<_>>(),
            vec![keys[1].clone(), keys[0].clone()]
        );
    }
}
