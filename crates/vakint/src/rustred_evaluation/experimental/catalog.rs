//! Deterministic persistence for the experimental four-loop terminal catalog.
//!
//! The catalog is deliberately only a value cache: it contains exact Symbolica
//! expressions for keys reached by an independently supplied candidate reducer.
//! It does not contain rules, certify a family, or grant an artifact authority.
//! Keeping this format small and line-oriented makes it possible to inspect and
//! review an offline FMFT export without adding a second serialization/CAS layer.

use std::collections::BTreeMap;

use rustred::family::IntegralKey;
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::utils::vakint_macros::vk_parse;

/// Version of the experimental offline terminal catalog format.
pub const SCHEMA_VERSION: u32 = 1;

const MAGIC: &str = "rustred-vakint-experimental-terminal-catalog";

/// Whether a catalog was imported for every declared candidate terminal or
/// only for the residuals reached by a finite probe matrix.
///
/// Partial catalogs remain useful for diagnostics and warm-cache experiments,
/// but must never be presented as a complete finite terminal declaration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CatalogCoverage {
    Partial,
    Complete,
}

/// Exact terminal values bound to one authenticated RustRed family identity.
///
/// This type is intentionally independent of `ClosedArtifact`: callers still
/// have to supply and validate the candidate reducer/artifact separately.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OfflineTerminalCatalog {
    family_fingerprint: String,
    index_count: usize,
    coverage: CatalogCoverage,
    terms: BTreeMap<IntegralKey, Atom>,
}

impl OfflineTerminalCatalog {
    /// Build a catalog from exact values obtained by an offline oracle.
    pub fn from_terms(
        family_fingerprint: impl Into<String>,
        index_count: usize,
        terms: BTreeMap<IntegralKey, Atom>,
    ) -> Result<Self, String> {
        let family_fingerprint = family_fingerprint.into();
        validate_identity(&family_fingerprint, index_count)?;
        let catalog = Self {
            family_fingerprint,
            index_count,
            coverage: CatalogCoverage::Partial,
            terms,
        };
        catalog.validate_terms()?;
        Ok(catalog)
    }

    fn from_terms_with_coverage(
        family_fingerprint: impl Into<String>,
        index_count: usize,
        coverage: CatalogCoverage,
        terms: BTreeMap<IntegralKey, Atom>,
    ) -> Result<Self, String> {
        let family_fingerprint = family_fingerprint.into();
        validate_identity(&family_fingerprint, index_count)?;
        let catalog = Self {
            family_fingerprint,
            index_count,
            coverage,
            terms,
        };
        catalog.validate_terms()?;
        Ok(catalog)
    }

    pub fn family_fingerprint(&self) -> &str {
        &self.family_fingerprint
    }

    pub fn index_count(&self) -> usize {
        self.index_count
    }

    pub fn coverage(&self) -> CatalogCoverage {
        self.coverage
    }

    pub fn is_complete(&self) -> bool {
        self.coverage == CatalogCoverage::Complete
    }

    /// Require declaration-complete coverage without making any closure claim.
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

    pub fn terms(&self) -> &BTreeMap<IntegralKey, Atom> {
        &self.terms
    }

    pub fn into_terms(self) -> BTreeMap<IntegralKey, Atom> {
        self.terms
    }

    /// Build a catalog for a complete, finite terminal declaration.
    ///
    /// `keys` is normally the immutable terminal set reported by a candidate
    /// reducer, rather than only the terminals reached by a small probe
    /// matrix.  Values are obtained through `resolve` at the offline import
    /// boundary (for example, from FMFT); the resulting catalog itself has no
    /// oracle dependency.  The map is sorted by `IntegralKey`, so encoding is
    /// deterministic even when the caller supplies an arbitrary iterator.
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
        Self::from_terms_with_coverage(
            family_fingerprint,
            index_count,
            CatalogCoverage::Complete,
            terms,
        )
    }

    /// Return the required keys that are absent from this catalog.
    ///
    /// This is intentionally a coverage check only.  It does not infer rules,
    /// promote a missing key to a master, or make any statement about family
    /// closure.
    pub fn missing_keys<I>(&self, required: I) -> Vec<IntegralKey>
    where
        I: IntoIterator<Item = IntegralKey>,
    {
        required
            .into_iter()
            .filter(|key| !self.terms.contains_key(key))
            .collect()
    }

    /// Fail with a concise coverage diagnostic if any declared terminal is
    /// absent.  Callers can use `missing_keys` when they need to report or
    /// recover the individual keys.
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

    /// Encode a deterministic, human-auditable representation.
    pub fn encode(&self) -> String {
        let mut output = String::new();
        output.push_str(MAGIC);
        output.push('\n');
        output.push_str(&format!("schema={}\n", SCHEMA_VERSION));
        output.push_str("family_fingerprint=");
        output.push_str(&self.family_fingerprint);
        output.push('\n');
        output.push_str(&format!("index_count={}\n", self.index_count));
        output.push_str("coverage=");
        output.push_str(match self.coverage {
            CatalogCoverage::Partial => "partial",
            CatalogCoverage::Complete => "complete",
        });
        output.push('\n');
        for (key, value) in &self.terms {
            output.push_str("terminal=");
            output.push_str(&format_powers(key));
            output.push('\t');
            output.push_str(&value.to_canonical_string());
            output.push('\n');
        }
        output
    }

    /// Decode and validate a catalog at the untrusted/offline load boundary.
    ///
    /// The expected family identity and arity are mandatory so a value file
    /// cannot accidentally be reused for a different integral family.
    pub fn decode(
        input: &str,
        expected_family_fingerprint: &str,
        expected_index_count: usize,
    ) -> Result<Self, String> {
        validate_identity(expected_family_fingerprint, expected_index_count)?;
        let mut lines = input.lines().peekable();
        if lines.next() != Some(MAGIC) {
            return Err("invalid terminal catalog magic".into());
        }
        let schema = required_field(&mut lines, "schema")?;
        let schema = schema
            .parse::<u32>()
            .map_err(|error| format!("invalid catalog schema {schema:?}: {error}"))?;
        if schema != SCHEMA_VERSION {
            return Err(format!(
                "unsupported terminal catalog schema {schema}; expected {SCHEMA_VERSION}"
            ));
        }
        let family_fingerprint = required_field(&mut lines, "family_fingerprint")?.to_owned();
        if family_fingerprint != expected_family_fingerprint {
            return Err(format!(
                "terminal catalog family fingerprint mismatch: got {family_fingerprint:?}, expected {expected_family_fingerprint:?}"
            ));
        }
        let index_count_field = required_field(&mut lines, "index_count")?;
        let index_count = index_count_field.parse::<usize>().map_err(|error| {
            format!("invalid terminal catalog index count {index_count_field:?}: {error}")
        })?;
        if index_count != expected_index_count {
            return Err(format!(
                "terminal catalog arity mismatch: got {index_count}, expected {expected_index_count}"
            ));
        }

        // Catalogs written before coverage metadata are explicitly treated as
        // partial.  This keeps old probe caches loadable while preventing them
        // from being mistaken for a complete finite terminal declaration.
        let coverage = match lines.peek().copied() {
            Some(line) if line.starts_with("coverage=") => {
                let value = lines.next().expect("peeked coverage line");
                match value.strip_prefix("coverage=") {
                    Some("partial") => CatalogCoverage::Partial,
                    Some("complete") => CatalogCoverage::Complete,
                    Some(other) => return Err(format!("invalid catalog coverage {other:?}")),
                    None => unreachable!("coverage prefix checked above"),
                }
            }
            _ => CatalogCoverage::Partial,
        };

        let mut terms = BTreeMap::new();
        for (line_number, line) in lines.enumerate() {
            let line_number = line_number + 5;
            let payload = line
                .strip_prefix("terminal=")
                .ok_or_else(|| format!("line {line_number} is not a terminal record"))?;
            let (powers, expression) = payload.split_once('\t').ok_or_else(|| {
                format!("line {line_number} terminal record is missing a TAB separator")
            })?;
            if expression.contains(['\n', '\r', '\t']) {
                return Err(format!(
                    "line {line_number} terminal expression contains a control character"
                ));
            }
            let powers = parse_powers(powers, expected_index_count)
                .map_err(|error| format!("line {line_number}: {error}"))?;
            let key = IntegralKey::try_new(powers)
                .map_err(|error| format!("line {line_number}: invalid terminal key: {error}"))?;
            if terms.contains_key(&key) {
                return Err(format!(
                    "line {line_number}: duplicate terminal key {:?}",
                    key.powers()
                ));
            }
            let value = vk_parse!(expression).map_err(|error| {
                format!("line {line_number}: invalid Symbolica expression: {error}")
            })?;
            validate_exact_value(&value, key.powers())
                .map_err(|error| format!("line {line_number}: {error}"))?;
            terms.insert(key, value);
        }
        Self::from_terms_with_coverage(family_fingerprint, index_count, coverage, terms)
    }

    fn validate_terms(&self) -> Result<(), String> {
        for (key, value) in &self.terms {
            if key.powers().len() != self.index_count {
                return Err(format!(
                    "terminal {:?} has arity {}, expected {}",
                    key.powers(),
                    key.powers().len(),
                    self.index_count
                ));
            }
            validate_exact_value(value, key.powers())?;
        }
        Ok(())
    }
}

fn validate_identity(fingerprint: &str, index_count: usize) -> Result<(), String> {
    if fingerprint.is_empty() || fingerprint.contains(['\n', '\r', '\t']) {
        return Err("family fingerprint must be nonempty and line-safe".into());
    }
    if index_count == 0 {
        return Err("terminal catalog index count must be positive".into());
    }
    Ok(())
}

fn validate_exact_value(value: &Atom, powers: &[i64]) -> Result<(), String> {
    let mut approximate = false;
    value.visitor(&mut |part| {
        if let AtomView::Num(number) = part {
            approximate |= number.get_coeff_view().is_float();
        }
        !approximate
    });
    if approximate {
        return Err(format!(
            "terminal {:?} contains an approximate coefficient; catalog values must be exact",
            powers
        ));
    }
    let rendered = value.to_canonical_string();
    if rendered.contains(['\n', '\r', '\t']) {
        return Err(format!(
            "terminal {:?} contains a control character",
            powers
        ));
    }
    Ok(())
}

fn format_powers(key: &IntegralKey) -> String {
    key.powers()
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn parse_powers(input: &str, expected_len: usize) -> Result<Vec<i64>, String> {
    if input.is_empty() {
        return Err("terminal key is empty".into());
    }
    let powers = input
        .split(',')
        .map(|power| {
            power
                .parse::<i64>()
                .map_err(|error| format!("invalid terminal power {power:?}: {error}"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if powers.len() != expected_len {
        return Err(format!(
            "terminal key has arity {}, expected {expected_len}",
            powers.len()
        ));
    }
    Ok(powers)
}

fn required_field<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    name: &str,
) -> Result<&'a str, String> {
    let line = lines
        .next()
        .ok_or_else(|| format!("terminal catalog is missing {name}"))?;
    line.strip_prefix(&format!("{name}="))
        .ok_or_else(|| format!("terminal catalog expected {name} field, got {line:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> OfflineTerminalCatalog {
        let key = IntegralKey::try_new([1, 1, 0, 0]).unwrap();
        OfflineTerminalCatalog::from_terms(
            "family-test-fingerprint",
            4,
            BTreeMap::from([(key, vk_parse!("PR4+PR9/2").unwrap())]),
        )
        .unwrap()
    }

    #[test]
    fn encoding_is_deterministic_and_round_trips() {
        let catalog = sample();
        let encoded = catalog.encode();
        assert_eq!(encoded, catalog.clone().encode());
        let decoded =
            OfflineTerminalCatalog::decode(&encoded, "family-test-fingerprint", 4).unwrap();
        assert_eq!(decoded, catalog);
        assert_eq!(decoded.coverage(), CatalogCoverage::Partial);
        assert!(!decoded.is_complete());
        assert!(decoded.require_complete().is_err());

        // Probe-era files had no coverage field.  They remain loadable but
        // are conservatively classified as partial rather than complete.
        let legacy = encoded
            .lines()
            .filter(|line| !line.starts_with("coverage="))
            .collect::<Vec<_>>()
            .join("\n");
        let legacy_decoded =
            OfflineTerminalCatalog::decode(&legacy, "family-test-fingerprint", 4).unwrap();
        assert_eq!(legacy_decoded.coverage(), CatalogCoverage::Partial);
    }

    #[test]
    fn decode_rejects_wrong_family_arity_and_approximate_values() {
        let encoded = sample().encode();
        assert!(
            OfflineTerminalCatalog::decode(&encoded, "other", 4)
                .unwrap_err()
                .contains("fingerprint mismatch")
        );
        assert!(
            OfflineTerminalCatalog::decode(&encoded, "family-test-fingerprint", 3)
                .unwrap_err()
                .contains("arity mismatch")
        );
        let header = encoded.lines().take(4).collect::<Vec<_>>().join("\n");
        let approximate = format!("{header}\nterminal=1,1,0,0\t1.25`40*PR4\n");
        assert!(
            OfflineTerminalCatalog::decode(&approximate, "family-test-fingerprint", 4)
                .unwrap_err()
                .contains("approximate")
        );
    }

    #[test]
    fn decode_rejects_duplicate_keys() {
        let catalog = sample();
        let mut encoded = catalog.encode();
        encoded.push_str("terminal=1,1,0,0\tPR1\n");
        assert!(
            OfflineTerminalCatalog::decode(&encoded, "family-test-fingerprint", 4)
                .unwrap_err()
                .contains("duplicate terminal key")
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
        // BTreeMap ownership keeps the serialized key order deterministic,
        // independent of the declaration order above.
        assert!(
            catalog.encode().find("terminal=1,1,0,0").unwrap()
                < catalog.encode().find("terminal=2,1,0,0").unwrap()
        );
    }
}
