use serde::{Deserialize, Serialize};
use thiserror::Error;

#[cfg(feature = "shadowing")]
use symbolica::atom::{Symbol, SymbolAttribute};

/// Prefix reserved for versioned Young-tableau metadata on tensor symbols.
pub const YOUNG_TABLEAU_TAG_PREFIX: &str = "spenso::young_tableau:";
const YOUNG_TABLEAU_TAG_VERSION: &str = "v1";

/// The intrinsic permutation behavior that Symbolica can represent directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum YoungTableauClass {
    Symmetric,
    Antisymmetric,
    General,
}

/// A Young diagram together with the tensor-slot order assigned to its boxes.
///
/// Rows are ordered from longest to shortest. `slot_order` lists the original
/// zero-based tensor slot placed in each box, reading rows from left to right
/// and top to bottom.
///
/// Version 1 metadata uses the normalized manifest Young symmetrizer
/// `P_T = C_T R_T / h_T`: `R_T` symmetrizes within the listed rows, `C_T`
/// antisymmetrizes within the corresponding manifest columns, and `h_T` is the
/// product of the diagram's hook lengths.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(try_from = "UncheckedYoungTableau")]
pub struct YoungTableau {
    shape: Vec<usize>,
    slot_order: Vec<usize>,
}

#[derive(Deserialize)]
struct UncheckedYoungTableau {
    shape: Vec<usize>,
    slot_order: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum YoungTableauError {
    #[error("a Young tableau shape must contain at least one row")]
    EmptyShape,
    #[error("Young tableau row {row} has length zero")]
    EmptyRow { row: usize },
    #[error(
        "Young tableau row lengths must be non-increasing, but row {row} has length {length} after length {previous}"
    )]
    NonPartition {
        row: usize,
        previous: usize,
        length: usize,
    },
    #[error("Young tableau rank overflows usize")]
    RankOverflow,
    #[error("cannot allocate storage for Young tableau rank {rank}")]
    Allocation { rank: usize },
    #[error("Young tableau has rank {expected}, but its slot order contains {actual} slots")]
    SlotCount { expected: usize, actual: usize },
    #[error(
        "Young tableau slot {slot} at position {position} is outside the rank-{rank} slot range"
    )]
    SlotOutOfRange {
        position: usize,
        slot: usize,
        rank: usize,
    },
    #[error("Young tableau slot {slot} occurs more than once")]
    DuplicateSlot { slot: usize },
    #[error("tag is not Young-tableau metadata")]
    NotYoungTableauTag,
    #[error("unsupported Young-tableau tag version `{version}`")]
    UnsupportedTagVersion { version: String },
    #[error("malformed Young-tableau tag")]
    MalformedTag,
    #[error("invalid integer `{value}` in Young-tableau tag {section}")]
    InvalidTagInteger {
        section: &'static str,
        value: String,
    },
    #[error("a symbol has {count} Young-tableau metadata tags; expected at most one")]
    MultipleSymbolTags { count: usize },
}

impl YoungTableau {
    pub fn new(shape: Vec<usize>, slot_order: Vec<usize>) -> Result<Self, YoungTableauError> {
        let rank = Self::validate_shape(&shape)?;
        if slot_order.len() != rank {
            return Err(YoungTableauError::SlotCount {
                expected: rank,
                actual: slot_order.len(),
            });
        }

        let mut seen = Vec::new();
        seen.try_reserve_exact(rank)
            .map_err(|_| YoungTableauError::Allocation { rank })?;
        seen.resize(rank, false);
        for (position, &slot) in slot_order.iter().enumerate() {
            let Some(entry) = seen.get_mut(slot) else {
                return Err(YoungTableauError::SlotOutOfRange {
                    position,
                    slot,
                    rank,
                });
            };
            if *entry {
                return Err(YoungTableauError::DuplicateSlot { slot });
            }
            *entry = true;
        }

        Ok(Self { shape, slot_order })
    }

    /// Construct a tableau whose boxes follow the tensor's existing slot order.
    pub fn canonical(shape: Vec<usize>) -> Result<Self, YoungTableauError> {
        let rank = Self::validate_shape(&shape)?;
        let mut slot_order = Vec::new();
        slot_order
            .try_reserve_exact(rank)
            .map_err(|_| YoungTableauError::Allocation { rank })?;
        slot_order.extend(0..rank);
        Self::new(shape, slot_order)
    }

    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    pub fn slot_order(&self) -> &[usize] {
        &self.slot_order
    }

    pub fn rank(&self) -> usize {
        self.slot_order.len()
    }

    /// Iterate over the tensor slots assigned to each row of the diagram.
    pub fn rows(&self) -> impl ExactSizeIterator<Item = &[usize]> {
        let mut start = 0;
        self.shape.iter().map(move |&length| {
            let row = &self.slot_order[start..start + length];
            start += length;
            row
        })
    }

    /// Iterate over the manifest columns, with slots ordered from top to bottom.
    pub fn columns(&self) -> impl ExactSizeIterator<Item = Vec<usize>> + '_ {
        (0..self.shape[0]).map(|column| {
            self.rows()
                .filter_map(|row| row.get(column).copied())
                .collect()
        })
    }

    /// Classify tableaux whose complete symmetry is representable by a
    /// Symbolica intrinsic attribute. The one-box tableau is treated as a row.
    pub fn class(&self) -> YoungTableauClass {
        if self.shape.len() == 1 {
            YoungTableauClass::Symmetric
        } else if self.shape.iter().all(|&length| length == 1) {
            YoungTableauClass::Antisymmetric
        } else {
            YoungTableauClass::General
        }
    }

    pub fn into_parts(self) -> (Vec<usize>, Vec<usize>) {
        (self.shape, self.slot_order)
    }

    /// Encode this tableau as deterministic, versioned symbol metadata.
    pub fn to_tag(&self) -> String {
        format!(
            "{YOUNG_TABLEAU_TAG_PREFIX}{YOUNG_TABLEAU_TAG_VERSION}:{}:{}",
            Self::encode_sequence(&self.shape),
            Self::encode_sequence(&self.slot_order)
        )
    }

    pub fn from_tag(tag: &str) -> Result<Self, YoungTableauError> {
        let payload = tag
            .strip_prefix(YOUNG_TABLEAU_TAG_PREFIX)
            .ok_or(YoungTableauError::NotYoungTableauTag)?;
        let mut fields = payload.split(':');
        let version = fields.next().ok_or(YoungTableauError::MalformedTag)?;
        if version != YOUNG_TABLEAU_TAG_VERSION {
            return Err(YoungTableauError::UnsupportedTagVersion {
                version: version.to_owned(),
            });
        }
        let shape = Self::decode_sequence(
            fields.next().ok_or(YoungTableauError::MalformedTag)?,
            "shape",
        )?;
        let slot_order = Self::decode_sequence(
            fields.next().ok_or(YoungTableauError::MalformedTag)?,
            "slot order",
        )?;
        if fields.next().is_some() {
            return Err(YoungTableauError::MalformedTag);
        }
        Self::new(shape, slot_order)
    }

    #[cfg(feature = "shadowing")]
    pub fn from_symbol(symbol: Symbol) -> Result<Option<Self>, YoungTableauError> {
        let tags = symbol
            .get_tags()
            .iter()
            .filter(|tag| tag.starts_with(YOUNG_TABLEAU_TAG_PREFIX))
            .collect::<Vec<_>>();
        match tags.as_slice() {
            [] => Ok(None),
            [tag] => Self::from_tag(tag).map(Some),
            tags => Err(YoungTableauError::MultipleSymbolTags { count: tags.len() }),
        }
    }

    #[cfg(feature = "shadowing")]
    pub fn symbol_attribute(&self) -> Option<SymbolAttribute> {
        match self.class() {
            YoungTableauClass::Symmetric => Some(SymbolAttribute::Symmetric),
            YoungTableauClass::Antisymmetric => Some(SymbolAttribute::Antisymmetric),
            YoungTableauClass::General => None,
        }
    }

    fn validate_shape(shape: &[usize]) -> Result<usize, YoungTableauError> {
        if shape.is_empty() {
            return Err(YoungTableauError::EmptyShape);
        }

        let mut rank = 0usize;
        for (row, &length) in shape.iter().enumerate() {
            if length == 0 {
                return Err(YoungTableauError::EmptyRow { row });
            }
            if row > 0 && length > shape[row - 1] {
                return Err(YoungTableauError::NonPartition {
                    row,
                    previous: shape[row - 1],
                    length,
                });
            }
            rank = rank
                .checked_add(length)
                .ok_or(YoungTableauError::RankOverflow)?;
        }
        Ok(rank)
    }

    fn encode_sequence(values: &[usize]) -> String {
        values
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(".")
    }

    fn decode_sequence(
        encoded: &str,
        section: &'static str,
    ) -> Result<Vec<usize>, YoungTableauError> {
        if encoded.is_empty() {
            return Err(YoungTableauError::MalformedTag);
        }
        encoded
            .split('.')
            .map(|value| {
                value
                    .parse()
                    .map_err(|_| YoungTableauError::InvalidTagInteger {
                        section,
                        value: value.to_owned(),
                    })
            })
            .collect()
    }
}

impl TryFrom<UncheckedYoungTableau> for YoungTableau {
    type Error = YoungTableauError;

    fn try_from(value: UncheckedYoungTableau) -> Result<Self, Self::Error> {
        Self::new(value.shape, value.slot_order)
    }
}

#[cfg(test)]
mod tests {
    use super::{YoungTableau, YoungTableauClass, YoungTableauError};

    #[test]
    fn validates_partition_and_slot_permutation() {
        assert!(matches!(
            YoungTableau::canonical(vec![]),
            Err(YoungTableauError::EmptyShape)
        ));
        assert!(matches!(
            YoungTableau::canonical(vec![1, 2]),
            Err(YoungTableauError::NonPartition { row: 1, .. })
        ));
        assert!(matches!(
            YoungTableau::new(vec![2, 1], vec![0, 1, 1]),
            Err(YoungTableauError::DuplicateSlot { slot: 1 })
        ));
        assert!(matches!(
            YoungTableau::new(vec![2, 1], vec![0, 1, 3]),
            Err(YoungTableauError::SlotOutOfRange { slot: 3, .. })
        ));
    }

    #[test]
    fn reports_unallocatable_rank() {
        assert_eq!(
            YoungTableau::canonical(vec![usize::MAX]),
            Err(YoungTableauError::Allocation { rank: usize::MAX })
        );
    }

    #[test]
    fn classifies_rows_columns_and_general_shapes() {
        assert_eq!(
            YoungTableau::canonical(vec![3]).unwrap().class(),
            YoungTableauClass::Symmetric
        );
        assert_eq!(
            YoungTableau::canonical(vec![1, 1, 1]).unwrap().class(),
            YoungTableauClass::Antisymmetric
        );
        assert_eq!(
            YoungTableau::canonical(vec![2, 1]).unwrap().class(),
            YoungTableauClass::General
        );
    }

    #[test]
    fn versioned_tag_round_trips_slot_order() {
        let tableau = YoungTableau::new(vec![2, 1], vec![2, 0, 1]).unwrap();
        let tag = tableau.to_tag();

        assert_eq!(tag, "spenso::young_tableau:v1:2.1:2.0.1");
        assert_eq!(YoungTableau::from_tag(&tag).unwrap(), tableau);
        assert_eq!(
            tableau.rows().collect::<Vec<_>>(),
            vec![&[2, 0][..], &[1][..]]
        );
        assert_eq!(
            tableau.columns().collect::<Vec<_>>(),
            vec![vec![2, 1], vec![0]]
        );
    }

    #[test]
    fn serde_round_trip_revalidates_tableau_invariants() {
        let tableau = YoungTableau::new(vec![2, 1], vec![2, 0, 1]).unwrap();
        let encoded = serde_json::to_string(&tableau).unwrap();
        assert_eq!(
            serde_json::from_str::<YoungTableau>(&encoded).unwrap(),
            tableau
        );

        let invalid = r#"{"shape":[1,2],"slot_order":[0,1,2]}"#;
        assert!(serde_json::from_str::<YoungTableau>(invalid).is_err());
    }
}
