use clap::Args;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Eq)]
pub struct StartCommandsBlock {
    /// Name under which the recorded command sequence will be stored.
    #[arg(value_name = "COMMAND_BLOCK_NAME")]
    pub name: String,
}
