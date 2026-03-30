use std::fmt::Display;

use dot_parser::ast::{CompassPt, Port};

use crate::half_edge::involution::{Flow, Hedge};

#[cfg(feature = "rkyv")]
use rkyv::{
    with::{ArchiveWith, DeserializeWith, SerializeWith},
    Archive, Archived, Deserialize, Fallible, Resolver, Serialize,
};

use super::{strip_quotes, subgraph_free::PortExt};

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
        compass_pt_to_u8(*field).resolve(pos, (), out);
    }
}

#[cfg(feature = "rkyv")]
impl<S: Fallible + ?Sized> SerializeWith<CompassPt, S> for CompassPtRkyv
where
    u8: Serialize<S>,
{
    fn serialize_with(field: &CompassPt, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
        compass_pt_to_u8(*field).serialize(serializer)
    }
}

#[cfg(feature = "rkyv")]
impl<D: Fallible + ?Sized> DeserializeWith<Archived<u8>, CompassPt, D> for CompassPtRkyv
where
    Archived<u8>: Deserialize<u8, D>,
{
    fn deserialize_with(field: &Archived<u8>, deserializer: &mut D) -> Result<CompassPt, D::Error> {
        Ok(u8_to_compass_pt(field.deserialize(deserializer)?))
    }
}

#[cfg(feature = "rkyv")]
pub(crate) const fn compass_pt_to_u8(value: CompassPt) -> u8 {
    match value {
        CompassPt::N => 0,
        CompassPt::NE => 1,
        CompassPt::E => 2,
        CompassPt::SE => 3,
        CompassPt::S => 4,
        CompassPt::SW => 5,
        CompassPt::W => 6,
        CompassPt::NW => 7,
        CompassPt::C => 8,
        CompassPt::Underscore => 9,
    }
}

#[cfg(feature = "rkyv")]
pub(crate) const fn u8_to_compass_pt(value: u8) -> CompassPt {
    match value {
        0 => CompassPt::N,
        1 => CompassPt::NE,
        2 => CompassPt::E,
        3 => CompassPt::SE,
        4 => CompassPt::S,
        5 => CompassPt::SW,
        6 => CompassPt::W,
        7 => CompassPt::NW,
        8 => CompassPt::C,
        9 => CompassPt::Underscore,
        _ => panic!("invalid archived compass point"),
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct DotHedgeData {
    pub statement: Option<String>,
    pub id: Option<Hedge>,
    pub port_label: Option<String>,
    #[cfg_attr(feature = "rkyv", with(rkyv::with::Map<CompassPtRkyv>))]
    pub compasspt: Option<CompassPt>,
}

impl DotHedgeData {
    pub fn is_none(&self) -> bool {
        self.statement.is_none()
            && self.id.is_none()
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

    pub fn with_port(mut self, port: &Port) -> Self {
        self.port_label = port.id().map(|a| a.to_string());
        self.compasspt = port.compass().cloned();
        self.id = port.hedge();
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
            statement: statement.map(|s| strip_quotes(&s).to_string()),
            ..Default::default()
        }
    }
}

impl From<Option<&Port>> for DotHedgeData {
    fn from(port: Option<&Port>) -> Self {
        if let Some(port) = port {
            DotHedgeData::default().with_port(port)
        } else {
            DotHedgeData::default()
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
