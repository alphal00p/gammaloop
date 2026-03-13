use clap::Args;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Args, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq, Eq)]
pub struct StartCommandsBlock {
    #[arg(value_name = "COMMAND_BLOCK_NAME")]
    pub name: String,
}
