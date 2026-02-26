use crate::utils::{
    GS,
    serde_utils::{
        IsDefault, is_default_form_path, is_default_pysecdec_relative_precision,
        is_default_python_path, is_default_vakint_evaluation_methods,
        is_default_vakint_normalization, is_false, is_one_string, is_true, is_usize,
    },
};
use bincode_trait_derive::{Decode, Encode};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use vakint::{EvaluationMethod, FMFTOptions, MATADOptions, PySecDecOptions, AlphaLoopOptions};

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
pub struct AlphaLoopSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub susbstitute_masters: bool,
}

impl Default for AlphaLoopSettings {
    fn default() -> Self {
        Self {
            susbstitute_masters: true,
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
    pub alphaloop: AlphaLoopSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub fmft: FMFTSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub pysecdec: PySecDecSettings,
    #[serde(skip_serializing_if = "is_usize::<16>")]
    pub run_time_decimal_precision: usize,
    #[serde(skip_serializing_if = "is_true")]
    pub clean_tmp_dir: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub temporary_directory: Option<String>,
    #[serde(skip_serializing_if = "is_default_vakint_normalization")]
    pub normalization: String,
    #[serde(skip_serializing_if = "is_one_string")]
    pub additional_normalization: String,
}

impl VakintSettings {
    pub fn true_settings(&self) -> vakint::VakintSettings {
        vakint::VakintSettings {
            form_exe_path: self.form_exe_path.clone(),
            python_exe_path: self.python_exe_path.clone(),
            verify_numerator_identification: false,
            run_time_decimal_precision: self.run_time_decimal_precision as u32,
            allow_unknown_integrals: false,
            clean_tmp_dir: self.clean_tmp_dir,
            evaluation_order: vakint::EvaluationOrder(
                self.evaluation_methods
                    .iter()
                    .map(|a| match a.as_str() {
                        "alphaloop" => EvaluationMethod::AlphaLoop(AlphaLoopOptions {
                            susbstitute_masters: self.alphaloop.susbstitute_masters,
                        }),
                        "matad" => EvaluationMethod::MATAD(MATADOptions {
                            expand_masters: self.matad.expand_masters,
                            susbstitute_masters: self.matad.susbstitute_masters,
                            substitute_hpls: self.matad.substitute_hpls,
                            direct_numerical_substition: self.matad.direct_numerical_substition,
                        }),
                        "fmft" => EvaluationMethod::FMFT(FMFTOptions {
                            expand_masters: self.fmft.expand_masters,
                            susbstitute_masters: self.fmft.susbstitute_masters,
                        }),
                        "pysecdec" => EvaluationMethod::PySecDec(PySecDecOptions {
                            quiet: self.pysecdec.quiet,
                            relative_precision: self.pysecdec.relative_precision,
                            min_n_evals: self.pysecdec.min_n_evals as u64,
                            max_n_evals: self.pysecdec.max_n_evals as u64,
                            reuse_existing_output: self.pysecdec.reuse_existing_output.clone(),
                            ..Default::default()
                        }),
                        _ => panic!("Unknown vakint evaluation method: {}", a),
                    })
                    .collect(),
            ),
            use_dot_product_notation: false,
            temporary_directory: self.temporary_directory.clone(),
            epsilon_symbol: GS.dim_epsilon.get_name().into(),
            mu_r_sq_symbol: GS.mu_r_sq.get_name().into(),
            integral_normalization_factor: match self.normalization.as_str() {
                "MSbar" => vakint::LoopNormalizationFactor::Custom("1".to_string()),
                "FMFTandMATAD" => vakint::LoopNormalizationFactor::FMFTandMATAD,
                "pySecDec" => vakint::LoopNormalizationFactor::pySecDec,
                _ => vakint::LoopNormalizationFactor::Custom(self.normalization.clone()),
            },
            //Custom("1".to_string()),
            number_of_terms_in_epsilon_expansion: 5,
            // ..Default::default()
        }
    }
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
            alphaloop: AlphaLoopSettings::default(),
            fmft: FMFTSettings::default(),
            pysecdec: PySecDecSettings::default(),
            run_time_decimal_precision: 100,
            clean_tmp_dir: true,
            temporary_directory: None,
            normalization: "MSbar".to_string(),
            additional_normalization: "1".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[serde(default, deny_unknown_fields)]
pub struct UVgenerationSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub softct: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub generate_integrated: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub only_integrated: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub add_sigma: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub inner_products: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub vakint: VakintSettings,
}

impl Default for UVgenerationSettings {
    fn default() -> Self {
        UVgenerationSettings {
            softct: true,
            generate_integrated: true,
            only_integrated: false,
            inner_products: true,
            add_sigma: false,
            vakint: VakintSettings::default(),
        }
    }
}
