use bincode_trait_derive::{Decode, Encode};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq)]
pub struct UVgenerationSettings {
    pub generate_integrated: bool,
}

impl Default for UVgenerationSettings {
    fn default() -> Self {
        UVgenerationSettings {
            generate_integrated: true,
        }
    }
}
