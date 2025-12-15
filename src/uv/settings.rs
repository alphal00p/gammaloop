use crate::utils::serde_utils::{
    IsDefault, is_default_form_path, is_default_pysecdec_relative_precision,
    is_default_python_path, is_default_vakint_evaluation_methods, is_default_vakint_normalization,
    is_true, is_usize,
};
use bincode_trait_derive::{Decode, Encode};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MATADSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub expand_masters: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub susbstitute_masters: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub substitute_hpls: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub direct_numerical_substition: bool,
}

impl Default for MATADSettings {
    fn default() -> Self {
        Self {
            expand_masters: true,
            susbstitute_masters: true,
            substitute_hpls: true,
            direct_numerical_substition: true,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FMFTSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub expand_masters: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub susbstitute_masters: bool,
}

impl Default for FMFTSettings {
    fn default() -> Self {
        Self {
            expand_masters: true,
            susbstitute_masters: true,
        }
    }
}
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PySecDecSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub quiet: bool,
    #[serde(skip_serializing_if = "is_default_pysecdec_relative_precision")]
    pub relative_precision: f64,
    #[serde(skip_serializing_if = "is_usize::<10_000>")]
    pub min_n_evals: usize,
    #[serde(skip_serializing_if = "is_usize::<1_000_000_000_000>")]
    pub max_n_evals: usize,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub reuse_existing_output: Option<String>,
}

impl Default for PySecDecSettings {
    fn default() -> Self {
        Self {
            quiet: true,
            relative_precision: 1.0e-7,
            min_n_evals: 10_000,
            max_n_evals: 1_000_000_000_000,
            reuse_existing_output: None,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct VakintSettings {
    #[serde(skip_serializing_if = "is_default_form_path")]
    pub form_exe_path: String,
    #[serde(skip_serializing_if = "is_default_python_path")]
    pub python_exe_path: String,
    #[serde(skip_serializing_if = "is_default_vakint_evaluation_methods")]
    pub evaluation_methods: Vec<String>,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub matad: MATADSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub fmft: FMFTSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub pysecdec: PySecDecSettings,
    #[serde(skip_serializing_if = "is_usize::<100>")]
    pub run_time_decimal_precision: usize,
    #[serde(skip_serializing_if = "is_true")]
    pub clean_tmp_dir: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub temporary_directory: Option<String>,
    #[serde(skip_serializing_if = "is_default_vakint_normalization")]
    pub normalization: String,
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
            matad: MATADSettings::default(),
            fmft: FMFTSettings::default(),
            pysecdec: PySecDecSettings::default(),
            run_time_decimal_precision: 100,
            clean_tmp_dir: true,
            temporary_directory: None,
            normalization: "MSbar".to_string(),
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
