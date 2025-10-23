use crate::utils::serde_utils::{
    is_default_form_path, is_default_python_path, is_default_vakint_evaluation_methods, is_true,
    is_usize, IsDefault,
};
use bincode_trait_derive::{Decode, Encode};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
pub struct VakintSettings {
    #[serde(skip_serializing_if = "is_default_form_path")]
    pub form_exe_path: String,
    #[serde(skip_serializing_if = "is_default_python_path")]
    pub python_exe_path: String,
    #[serde(skip_serializing_if = "is_default_vakint_evaluation_methods")]
    pub evaluation_methods: Vec<String>,
    #[serde(skip_serializing_if = "is_usize::<100>")]
    pub run_time_decimal_precision: usize,
    #[serde(skip_serializing_if = "is_true")]
    pub clean_tmp_dir: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub temporary_directory: Option<String>,
}

impl Default for VakintSettings {
    fn default() -> Self {
        Self {
            form_exe_path: "form".to_string(),
            python_exe_path: "python3".to_string(),
            // Supported evaluation methods: "alphaloop", "matad", "fmft" and "pysecdec"
            evaluation_methods: vec![
                "alphaloop".to_string(),
                "matad".to_string(),
                "fmft".to_string(),
            ],
            run_time_decimal_precision: 100,
            clean_tmp_dir: true,
            temporary_directory: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[serde(default, deny_unknown_fields)]
pub struct UVgenerationSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub generate_integrated: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub vakint: VakintSettings,
}

impl Default for UVgenerationSettings {
    fn default() -> Self {
        UVgenerationSettings {
            generate_integrated: true,
            vakint: VakintSettings::default(),
        }
    }
}
