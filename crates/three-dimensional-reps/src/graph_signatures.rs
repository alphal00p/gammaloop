use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct MomentumSignature {
    pub loop_signature: Vec<i32>,
    pub external_signature: Vec<i32>,
}

impl MomentumSignature {
    pub fn negated(&self) -> Self {
        Self {
            loop_signature: self.loop_signature.iter().map(|x| -*x).collect(),
            external_signature: self.external_signature.iter().map(|x| -*x).collect(),
        }
    }

    pub fn canonical_up_to_sign(&self) -> (Self, i32) {
        let negated = self.negated();
        if self <= &negated {
            (self.clone(), 1)
        } else {
            (negated, -1)
        }
    }
}
