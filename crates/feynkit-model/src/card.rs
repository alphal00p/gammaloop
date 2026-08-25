use std::{
    collections::BTreeMap,
    fs,
    ops::{Deref, DerefMut},
    path::Path,
};

use serde::{Deserialize, Serialize};

use crate::ModelError;

/// A complex number encoded by UFO JSON as `[real, imaginary]`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(from = "(f64, f64)", into = "(f64, f64)")]
pub struct ComplexValue {
    pub re: f64,
    pub im: f64,
}

impl ComplexValue {
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };

    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn is_zero(self) -> bool {
        self == Self::ZERO
    }
}

impl From<(f64, f64)> for ComplexValue {
    fn from((re, im): (f64, f64)) -> Self {
        Self { re, im }
    }
}

impl From<ComplexValue> for (f64, f64) {
    fn from(value: ComplexValue) -> Self {
        (value.re, value.im)
    }
}

/// Values supplied for named model parameters.
///
/// Normalized restriction cards may include internal parameters whose values
/// were fixed while simplifying a model.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct ParameterCard(BTreeMap<String, ComplexValue>);

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
