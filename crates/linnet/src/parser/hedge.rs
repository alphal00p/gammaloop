use std::fmt::Display;

use crate::half_edge::involution::{Flow, Hedge};

#[cfg(feature = "rkyv")]
use rkyv::{
    with::{ArchiveWith, DeserializeWith, SerializeWith},
    Archive, Archived, Deserialize, Fallible, Resolver, Serialize,
};

use super::{syntax::Port, CompassPt, DotParseError, ExplicitIdKind};

/// Archives a compass point as its discriminant byte.
#[cfg(feature = "rkyv")]
struct CompassPtRkyv;

#[cfg(feature = "rkyv")]
impl ArchiveWith<CompassPt> for CompassPtRkyv {
    type Archived = Archived<u8>;
    type Resolver = Resolver<u8>;

    unsafe fn resolve_with(
        field: &CompassPt,
        pos: usize,
        _: Self::Resolver,
        out: *mut Self::Archived,
    ) {
        (*field as u8).resolve(pos, (), out);
    }
}

#[cfg(feature = "rkyv")]
impl<S: Fallible + ?Sized> SerializeWith<CompassPt, S> for CompassPtRkyv
where
    u8: Serialize<S>,
{
    fn serialize_with(field: &CompassPt, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        (*field as u8).serialize(serializer)
    }
}

#[cfg(feature = "rkyv")]
impl<D: Fallible + ?Sized> DeserializeWith<Archived<u8>, CompassPt, D> for CompassPtRkyv
where
    Archived<u8>: Deserialize<u8, D>,
{
    fn deserialize_with(field: &Archived<u8>, deserializer: &mut D) -> Result<CompassPt, D::Error> {
        let byte: u8 = field.deserialize(deserializer)?;
        Ok(CompassPt::try_from(byte).expect("invalid archived compass point"))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub struct DotHedgeData {
    pub statement: Option<String>,
    pub id: Option<Hedge>,
    pub payload: Option<Vec<u8>>,
    pub port_label: Option<String>,
    #[cfg_attr(feature = "rkyv", with(rkyv::with::Map<CompassPtRkyv>))]
    pub compasspt: Option<CompassPt>,
}

impl DotHedgeData {
    /// The half-edge data of an endpoint port. A port name that parses as a
    /// nonnegative integer requests that half-edge index; any other integer
    /// spelling is an invalid index.
    pub(crate) fn from_parser(port: Option<&Port>) -> Result<Self, DotParseError> {
        let Some(port) = port else {
            return Ok(Self::default());
        };
        let name = port.name();
        let data = Self {
            port_label: name.map(str::to_owned),
            compasspt: port.compass(),
            id: name.and_then(|name| name.parse().ok()),
            ..Self::default()
        };
        if let (None, Some(value)) = (data.id, name) {
            let digits = value
                .strip_prefix('+')
                .or_else(|| value.strip_prefix('-'))
                .unwrap_or(value);
            if !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit()) {
                if let Err(source) = value.parse::<usize>() {
                    return Err(DotParseError::InvalidExplicitId {
                        kind: ExplicitIdKind::Hedge,
                        value: value.to_string(),
                        source,
                    });
                }
            }
        }
        Ok(data)
    }

    pub fn is_none(&self) -> bool {
        self.statement.is_none()
            && self.id.is_none()
            && self.payload.is_none()
            && self.port_label.is_none()
            && self.compasspt.is_none()
    }

    #[allow(clippy::useless_asref)]
    pub fn dot_serialize(&self) -> Option<String> {
        let mut out = String::new();
        let mut info = false;
        if let Some(statement) = &self.statement {
            info = true;
            out.push_str(statement);
        }

        if let Some(id) = &self.id {
            info = true;
            out.push_str(&format!(" [id={id}]"));
        }
        if info {
            Some(out)
        } else {
            None
        }
    }

    pub fn with_statement(mut self, statement: String) -> Self {
        self.statement = Some(statement);
        self
    }

    pub fn with_id(mut self, id: Hedge) -> Self {
        self.id = Some(id);
        self
    }
}

impl Display for DotHedgeData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(statement) = &self.statement {
            write!(f, "{statement}")?;
        }
        Ok(())
    }
}

impl From<Option<String>> for DotHedgeData {
    fn from(statement: Option<String>) -> Self {
        DotHedgeData {
            statement,
            ..Default::default()
        }
    }
}

pub enum ParsingHedgePair {
    Unpaired {
        hedge: Option<Hedge>,
        flow: Flow,
        data: DotHedgeData,
    },
    Paired {
        source: Option<Hedge>,
        source_data: DotHedgeData,
        sink: Option<Hedge>,
        sink_data: DotHedgeData,
    },
}
