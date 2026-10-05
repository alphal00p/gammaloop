//! Read an exact JSON token, or a scalar supplied by another Serde format.

use std::fmt;

use serde::de::{Error, MapAccess, Visitor};
use serde::{Deserialize, Deserializer, forward_to_deserialize_any};
use serde_json::value::RawValue;

pub(super) struct RawScalar(pub Box<RawValue>);

impl<'de> Deserialize<'de> for RawScalar {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Box::<RawValue>::deserialize(ScalarDeserializer(d)).map(Self)
    }
}

// Ask serde_json for its raw token through RawValue's public implementation.
// Other human-readable formats typically supply a transparent newtype or scalar.
struct ScalarDeserializer<D>(D);

impl<'de, D: Deserializer<'de>> Deserializer<'de> for ScalarDeserializer<D> {
    type Error = D::Error;

    fn deserialize_any<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value, D::Error> {
        self.0.deserialize_any(visitor)
    }

    fn deserialize_newtype_struct<V: Visitor<'de>>(
        self,
        name: &'static str,
        visitor: V,
    ) -> Result<V::Value, D::Error> {
        self.0
            .deserialize_newtype_struct(name, RawOrScalar { name, visitor })
    }

    forward_to_deserialize_any! {
        bool i8 i16 i32 i64 i128 u8 u16 u32 u64 u128 f32 f64 char str string
        bytes byte_buf option unit unit_struct seq tuple tuple_struct map struct
        enum identifier ignored_any
    }
}

struct RawOrScalar<V> {
    name: &'static str,
    visitor: V,
}

impl<'de, V: Visitor<'de>> RawOrScalar<V> {
    fn scalar<E: Error>(self, value: impl Into<serde_json::Value>) -> Result<V::Value, E> {
        value
            .into()
            .deserialize_newtype_struct(self.name, self.visitor)
            .map_err(E::custom)
    }
}

impl<'de, V: Visitor<'de>> Visitor<'de> for RawOrScalar<V> {
    type Value = V::Value;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("a number or an exact fraction string")
    }

    fn visit_map<A: MapAccess<'de>>(self, map: A) -> Result<Self::Value, A::Error> {
        self.visitor.visit_map(map)
    }

    fn visit_newtype_struct<D: Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        self.scalar(serde_json::Value::deserialize(d)?)
    }

    fn visit_i64<E: Error>(self, value: i64) -> Result<Self::Value, E> {
        self.scalar(value)
    }

    fn visit_u64<E: Error>(self, value: u64) -> Result<Self::Value, E> {
        self.scalar(value)
    }

    fn visit_f64<E: Error>(self, value: f64) -> Result<Self::Value, E> {
        self.scalar(value)
    }

    fn visit_str<E: Error>(self, value: &str) -> Result<Self::Value, E> {
        self.scalar(value)
    }
}
