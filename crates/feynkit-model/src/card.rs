use std::{
    collections::BTreeMap,
    fs,
    ops::{Deref, DerefMut},
    path::Path,
};

use serde::{Deserialize, Serialize};
use symbolica::domains::float::Complex;

use crate::ModelError;

/// The complex scalar used by model parameters and couplings.
pub type ComplexValue = Complex<f64>;

pub(crate) mod optional_complex {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    use super::ComplexValue;

    pub fn serialize<S>(value: &Option<ComplexValue>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        value
            .map(|value| (value.re, value.im))
            .serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Option<ComplexValue>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Ok(Option::<(f64, f64)>::deserialize(deserializer)?
            .map(|(re, im)| ComplexValue::new(re, im)))
    }
}

mod complex_map {
    use std::collections::BTreeMap;

    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    use super::ComplexValue;

    pub fn serialize<S>(
        values: &BTreeMap<String, ComplexValue>,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        values
            .iter()
            .map(|(name, value)| (name, (value.re, value.im)))
            .collect::<BTreeMap<_, _>>()
            .serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<BTreeMap<String, ComplexValue>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Ok(BTreeMap::<String, (f64, f64)>::deserialize(deserializer)?
            .into_iter()
            .map(|(name, (re, im))| (name, ComplexValue::new(re, im)))
            .collect())
    }
}

/// Values supplied for named model parameters.
///
/// Normalized restriction cards may include internal parameters whose values
/// were fixed while simplifying a model.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct ParameterCard(#[serde(with = "complex_map")] BTreeMap<String, ComplexValue>);

impl ParameterCard {
    /// Create an empty card.
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_json(json: &str) -> Result<Self, ModelError> {
        Ok(serde_json::from_str(json)?)
    }

    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, ModelError> {
        let path = path.as_ref();
        let json = fs::read_to_string(path).map_err(|source| ModelError::Read {
            path: path.to_path_buf(),
            source,
        })?;
        Self::from_json(&json)
    }

    pub fn to_json(&self) -> Result<String, ModelError> {
        Ok(serde_json::to_string(self)?)
    }

    pub fn to_json_pretty(&self) -> Result<String, ModelError> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<(), ModelError> {
        let path = path.as_ref();
        fs::write(path, self.to_json_pretty()?).map_err(|source| ModelError::Write {
            path: path.to_path_buf(),
            source,
        })
    }
}

impl Deref for ParameterCard {
    type Target = BTreeMap<String, ComplexValue>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for ParameterCard {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl IntoIterator for ParameterCard {
    type Item = (String, ComplexValue);
    type IntoIter = std::collections::btree_map::IntoIter<String, ComplexValue>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl FromIterator<(String, ComplexValue)> for ParameterCard {
    fn from_iter<T: IntoIterator<Item = (String, ComplexValue)>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}
