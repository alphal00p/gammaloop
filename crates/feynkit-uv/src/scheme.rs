//! Renormalization schemes and the rules that assign them to UV subgraphs.

use std::collections::BTreeSet;
use std::fmt::{self, Display, Formatter};

use bincode_trait_derive::{Decode, Encode};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize, Serializer};

/// The approximation, or renormalization scheme, applied to one UV spinney:
/// `MUV` subtracts with an auxiliary UV mass, `PolePart` subtracts the
/// integrated pole only, `OS` subtracts on shell, `IR` applies the completed
/// soft prescription and `Unsubtracted` leaves the subgraph out of the wood.
///
/// ```
/// use feynkit_uv::ApproximationType;
///
/// let scheme: ApproximationType = serde_json::from_str("\"pole_part\"").unwrap();
/// assert_eq!(scheme, ApproximationType::PolePart);
/// assert_eq!(scheme.to_string(), "PolePart");
/// ```
#[derive(
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    Encode,
    Decode,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    JsonSchema,
    Default,
)]
#[cfg_attr(feature = "python", pyo3::pyclass(from_py_object))]
#[schemars(description = "Renormalization scheme applied to one UV subgraph.")]
#[serde(deny_unknown_fields)]
pub enum ApproximationType {
    #[default]
    #[serde(rename = "MUV", alias = "muv")]
    MUV,
    #[serde(rename = "PolePart", alias = "pole_part")]
    PolePart,
    #[serde(rename = "OS", alias = "os")]
    OS,
    #[serde(rename = "IR", alias = "ir")]
    IR,
    #[serde(rename = "Unsubtracted", alias = "unsubtracted")]
    Unsubtracted,
}

impl Display for ApproximationType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            ApproximationType::MUV => write!(f, "MUV"),
            ApproximationType::PolePart => write!(f, "PolePart"),
            ApproximationType::OS => write!(f, "OS"),
            ApproximationType::IR => write!(f, "IR"),
            ApproximationType::Unsubtracted => write!(f, "Unsubtracted"),
        }
    }
}

/// Particle content of a counterterm: the signed PDG codes on its boundary and,
/// optionally, the PDG codes of its internal propagators.
///
/// Boundary codes are negated on half-edges whose flow leaves the subgraph, so
/// a particle and its antiparticle are distinguished by direction.
///
/// ```
/// use std::collections::BTreeSet;
/// use feynkit_uv::CTIdentifier;
///
/// let quark_self_energy = CTIdentifier::new(BTreeSet::from([-1, 1]), Some(BTreeSet::from([1, 21])));
/// let any_internal = CTIdentifier::new(BTreeSet::from([-1, 1]), None);
/// assert!(any_internal.matches(&quark_self_energy));
/// assert!(!quark_self_energy.matches(&any_internal));
/// ```
#[cfg_attr(feature = "python", pyo3::pyclass(from_py_object, get_all, set_all))]
#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    Encode,
    Decode,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    JsonSchema,
)]
#[schemars(description = "Boundary and optional internal PDG codes that select counterterms.")]
#[serde(default, deny_unknown_fields)]
pub struct CTIdentifier {
    /// Signed PDG identifiers on the counterterm's external boundary.
    #[serde(skip_serializing_if = "BTreeSet::is_empty")]
    pub external_pdg_set: BTreeSet<isize>,
    /// Optional internal PDG identifiers; omitting them makes the rule match any internal set.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub internal_pdg_set: Option<BTreeSet<isize>>,
}

impl CTIdentifier {
    /// Identify a counterterm by its boundary and, optionally, internal particle content.
    ///
    /// See [`CTIdentifier`] for an example.
    pub fn new(
        external_pdg_set: BTreeSet<isize>,
        internal_pdg_set: Option<BTreeSet<isize>>,
    ) -> Self {
        Self {
            external_pdg_set,
            internal_pdg_set,
        }
    }

    /// Whether this pattern selects `candidate`: the boundaries agree and either
    /// this pattern leaves the internal content open or both internal sets agree.
    ///
    /// See [`CTIdentifier`] for an example.
    pub fn matches(&self, candidate: &Self) -> bool {
        self.external_pdg_set == candidate.external_pdg_set
            && self
                .internal_pdg_set
                .as_ref()
                .is_none_or(|internal| candidate.internal_pdg_set.as_ref() == Some(internal))
    }
}

/// A scheme override for the counterterms selected by a [`CTIdentifier`].
///
/// ```
/// use std::collections::BTreeSet;
/// use feynkit_uv::{ApproximationType, CTIdentifier, CTRenormalizationRule};
///
/// let rule = CTRenormalizationRule::new(
///     CTIdentifier::new(BTreeSet::from([21, -21]), None),
///     ApproximationType::OS,
/// );
/// assert_eq!(rule.prescription, ApproximationType::OS);
/// ```
#[cfg_attr(feature = "python", pyo3::pyclass(from_py_object, get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema)]
#[schemars(description = "Scheme override for the counterterms matching an identifier.")]
#[serde(deny_unknown_fields)]
pub struct CTRenormalizationRule {
    /// Particle-content pattern that selects the counterterms to which this rule applies.
    pub ct_identifier: CTIdentifier,
    /// Approximation that replaces the default prescription for matching counterterms.
    pub prescription: ApproximationType,
}

impl CTRenormalizationRule {
    /// Override the prescription of the counterterms matched by `ct_identifier`.
    ///
    /// See [`CTRenormalizationRule`] for an example.
    pub fn new(ct_identifier: CTIdentifier, prescription: ApproximationType) -> Self {
        Self {
            ct_identifier,
            prescription,
        }
    }
}

/// Default and counterterm-specific renormalization schemes.
///
/// Serialization omits fields equal to their defaults. Hosts that display every
/// setting use [`Self::serialize_with_defaults`] instead.
///
/// ```
/// use std::collections::BTreeSet;
/// use feynkit_uv::{ApproximationType, CTIdentifier, RenormalizationPrescriptionSettings};
///
/// let settings = RenormalizationPrescriptionSettings {
///     massless_power_divergent: ApproximationType::IR,
///     ..Default::default()
/// };
/// let gluon_self_energy = CTIdentifier::new(BTreeSet::from([21]), None);
/// assert_eq!(
///     settings.approximation_scheme_for(&gluon_self_energy, 2, false),
///     ApproximationType::IR
/// );
/// assert_eq!(
///     settings.approximation_scheme_for(&gluon_self_energy, 0, false),
///     ApproximationType::MUV
/// );
/// ```
#[cfg_attr(feature = "python", pyo3::pyclass(from_py_object, get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[schemars(description = "Default and counterterm-specific renormalization schemes.")]
#[serde(default, deny_unknown_fields)]
pub struct RenormalizationPrescriptionSettings {
    /// Approximation applied to logarithmically divergent ultraviolet counterterms.
    #[serde(skip_serializing_if = "is_default")]
    pub log_divergent: ApproximationType,
    /// Approximation applied to power-divergent counterterms with massive external structure.
    #[serde(skip_serializing_if = "is_default")]
    pub massive_power_divergent: ApproximationType,
    /// Approximation applied to power-divergent counterterms with massless external structure.
    #[serde(skip_serializing_if = "is_default")]
    pub massless_power_divergent: ApproximationType,
    /// Counterterm-specific prescription rules evaluated before the divergence-class defaults.
    #[serde(
        default,
        skip_serializing_if = "Vec::is_empty",
        with = "ct_renormalization_overrides_serde"
    )]
    #[schemars(with = "Vec<CTRenormalizationRule>")]
    pub overrides: Vec<CTRenormalizationRule>,
}

impl Default for RenormalizationPrescriptionSettings {
    fn default() -> Self {
        Self {
            log_divergent: ApproximationType::MUV,
            massive_power_divergent: ApproximationType::MUV,
            massless_power_divergent: ApproximationType::MUV,
            overrides: Vec::new(),
        }
    }
}

/// Every field of [`RenormalizationPrescriptionSettings`], including defaults.
#[derive(Serialize)]
struct PrescriptionWithDefaults<'a> {
    log_divergent: ApproximationType,
    massive_power_divergent: ApproximationType,
    massless_power_divergent: ApproximationType,
    overrides: &'a [CTRenormalizationRule],
}

impl RenormalizationPrescriptionSettings {
    /// Select the scheme of a divergent counterterm.
    ///
    /// An override with exactly this identifier wins; otherwise the first override
    /// whose internal content is open and whose boundary matches. Without a matching
    /// override, logarithmic divergences (`dod == 0`) use [`Self::log_divergent`],
    /// and power divergences choose by whether the counterterm has a massive
    /// boundary particle.
    ///
    /// See [`RenormalizationPrescriptionSettings`] for an example.
    pub fn approximation_scheme_for(
        &self,
        ct_identifier: &CTIdentifier,
        dod: i32,
        has_massive_externals: bool,
    ) -> ApproximationType {
        if let Some(rule) = self
            .overrides
            .iter()
            .find(|rule| &rule.ct_identifier == ct_identifier)
        {
            return rule.prescription;
        }

        if let Some(rule) = self.overrides.iter().find(|rule| {
            rule.ct_identifier.internal_pdg_set.is_none()
                && rule.ct_identifier.matches(ct_identifier)
        }) {
            return rule.prescription;
        }

        if dod == 0 {
            self.log_divergent
        } else if has_massive_externals {
            self.massive_power_divergent
        } else {
            self.massless_power_divergent
        }
    }

    /// Serialize every field, including those equal to their defaults.
    ///
    /// # Errors
    ///
    /// Returns the serializer's error.
    ///
    /// ```
    /// use feynkit_uv::RenormalizationPrescriptionSettings;
    ///
    /// let defaults = RenormalizationPrescriptionSettings::default();
    /// assert_eq!(serde_json::to_string(&defaults).unwrap(), "{}");
    /// let mut full = Vec::new();
    /// defaults
    ///     .serialize_with_defaults(&mut serde_json::Serializer::new(&mut full))
    ///     .unwrap();
    /// assert!(String::from_utf8(full).unwrap().contains("\"log_divergent\":\"MUV\""));
    /// ```
    pub fn serialize_with_defaults<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        PrescriptionWithDefaults {
            log_divergent: self.log_divergent,
            massive_power_divergent: self.massive_power_divergent,
            massless_power_divergent: self.massless_power_divergent,
            overrides: &self.overrides,
        }
        .serialize(serializer)
    }
}

fn is_default<T: Default + PartialEq>(value: &T) -> bool {
    value == &T::default()
}

mod ct_renormalization_overrides_serde {
    use std::collections::BTreeSet;

    use serde::de::Error;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    use super::CTRenormalizationRule;

    pub fn serialize<S>(
        overrides: &[CTRenormalizationRule],
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        overrides.serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Vec<CTRenormalizationRule>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let rules = Vec::<CTRenormalizationRule>::deserialize(deserializer)?;
        let mut seen = BTreeSet::new();

        for rule in &rules {
            if !seen.insert(rule.ct_identifier.clone()) {
                return Err(D::Error::custom(format!(
                    "duplicate CTIdentifier in renormalization_prescription.overrides: {:?}",
                    rule.ct_identifier
                )));
            }
        }

        Ok(rules)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_suppression_and_full_serialization_agree_on_values() {
        let settings = RenormalizationPrescriptionSettings {
            log_divergent: ApproximationType::PolePart,
            ..Default::default()
        };
        let compact = toml::to_string(&settings).unwrap();
        assert_eq!(compact.trim(), "log_divergent = \"PolePart\"");

        let mut full = Vec::new();
        settings
            .serialize_with_defaults(&mut serde_json::Serializer::new(&mut full))
            .unwrap();
        let full = String::from_utf8(full).unwrap();
        assert_eq!(
            full,
            r#"{"log_divergent":"PolePart","massive_power_divergent":"MUV","massless_power_divergent":"MUV","overrides":[]}"#
        );
        assert_eq!(
            serde_json::from_str::<RenormalizationPrescriptionSettings>(&full).unwrap(),
            settings
        );
    }
}
