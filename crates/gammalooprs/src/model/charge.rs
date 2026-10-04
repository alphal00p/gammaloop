use serde::{Deserialize, Deserializer, Serialize, Serializer, de::Error};
use symbolica::domains::{integer::Integer, rational::Rational};

/// UFO 1.0 exports integral charges as numbers and fractional charges as strings.
/// Keep them exact through model import/export; older numeric JSON remains valid.
#[derive(Debug, Clone)]
pub(super) struct SerializedCharge(pub Rational);

impl Serialize for SerializedCharge {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        if self.0.denominator_ref().is_one()
            && let Some(value) = self.0.numerator_ref().to_i64()
        {
            serializer.serialize_i64(value)
        } else {
            serializer.serialize_str(&self.0.to_string())
        }
    }
}

impl<'de> Deserialize<'de> for SerializedCharge {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Input {
            Rational(String),
            Number(serde_json::Number),
        }

        let value = match Input::deserialize(deserializer)? {
            Input::Rational(text) => {
                let (num, den) = text.split_once('/').unwrap_or((&text, "1"));
                let num = num.parse::<Integer>().map_err(D::Error::custom)?;
                let den = den.parse::<Integer>().map_err(D::Error::custom)?;
                if den.is_zero() {
                    return Err(D::Error::custom("particle charge has a zero denominator"));
                }
                Rational::new(num, den)
            }
            Input::Number(number) => {
                if let Some(integer) = number.as_i64() {
                    integer.into()
                } else if let Some(integer) = number.as_u64() {
                    integer.into()
                } else {
                    Rational::try_from(number.as_f64().ok_or_else(|| {
                        D::Error::custom("particle charge must be a finite number")
                    })?)
                    .map_err(D::Error::custom)?
                }
            }
        };
        Ok(Self(value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_and_numeric_particle_charges_roundtrip() {
        for (input, expected) in [
            ("0", Rational::zero()),
            ("-1", Rational::from(-1)),
            ("0.5", Rational::from((1, 2))),
            (r#""2/3""#, Rational::from((2, 3))),
            (r#""-4/3""#, Rational::from((-4, 3))),
        ] {
            let charge: SerializedCharge = serde_json::from_str(input).unwrap();
            assert_eq!(charge.0, expected);
            let encoded = serde_json::to_string(&charge).unwrap();
            let decoded: SerializedCharge = serde_json::from_str(&encoded).unwrap();
            assert_eq!(decoded.0, expected);
        }
        for input in [r#""1/0""#, r#""Q""#, "null", "true"] {
            assert!(
                serde_json::from_str::<SerializedCharge>(input).is_err(),
                "{input}"
            );
        }
    }
}
