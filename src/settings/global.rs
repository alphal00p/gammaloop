use std::{fs::File, path::Path};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Result, Section};
use eyre::Context;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    evaluate::{CompileOptions, InlineASM},
    function,
};

use crate::_default_true;
use crate::{
    cff::expression::GraphOrientation,
    processes::EvaluatorSettings,
    utils::{serde_utils::IsDefault, symbolica_ext::StringSerializedAtom, GS, W_},
    uv::UVgenerationSettings,
    GammaLoopContext,
};

use super::GlobalSettings;

impl GlobalSettings {
    pub(crate) fn from_file(filename: impl AsRef<Path>) -> Result<Self> {
        let filename = filename.as_ref();
        let f = File::open(filename)
            .wrap_err_with(|| {
                format!(
                    "Could not open generation settings file {}",
                    filename.display()
                )
            })
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .wrap_err("Could not parse generation  settings file")
            .suggestion("Is it a correct yaml file")
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GenerationSettings {
    // Generation Time settings
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub evaluator_settings: EvaluatorSettings,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub orientation_pattern: OrientationPattern,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub gammaloop_compile_options: GammaloopCompileOptions,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub tropical_subgraph_table_settings: TropicalSubgraphTableSettings,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")] // default false
    pub enable_thresholds: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub uv_settings: UVgenerationSettings,
}

impl Default for GenerationSettings {
    fn default() -> Self {
        Self {
            evaluator_settings: EvaluatorSettings::default(),
            orientation_pattern: OrientationPattern::default(),
            gammaloop_compile_options: GammaloopCompileOptions::default(),
            tropical_subgraph_table_settings: TropicalSubgraphTableSettings::default(),
            enable_thresholds: false,
            uv_settings: UVgenerationSettings::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Copy, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
pub enum CompilationOptimizationLevel {
    O0,
    O1,
    O2,
    O3,
}

impl Default for CompilationOptimizationLevel {
    fn default() -> Self {
        CompilationOptimizationLevel::O3
    }
}

impl From<CompilationOptimizationLevel> for usize {
    fn from(value: CompilationOptimizationLevel) -> Self {
        match value {
            CompilationOptimizationLevel::O0 => 0,
            CompilationOptimizationLevel::O1 => 1,
            CompilationOptimizationLevel::O2 => 2,
            CompilationOptimizationLevel::O3 => 3,
        }
    }
}

fn _default_compiler() -> String {
    "g++".to_owned()
}

pub const fn yes() -> bool {
    true
}

pub fn gpp() -> String {
    "g++".to_owned()
}

pub fn is_gpp(compiler: &str) -> bool {
    "g++" == compiler
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
pub struct GammaloopCompileOptions {
    #[serde(default = "yes", skip_serializing_if = "std::clone::Clone::clone")]
    pub inline_asm: bool,

    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub optimization_level: CompilationOptimizationLevel,

    #[serde(default = "yes", skip_serializing_if = "std::clone::Clone::clone")] // default true
    pub fast_math: bool,

    #[serde(default = "yes", skip_serializing_if = "std::clone::Clone::clone")] // default true
    pub unsafe_math: bool,

    #[serde(default = "gpp", skip_serializing_if = "is_gpp")] // default g++
    pub compiler: String,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub custom: Vec<String>,
}

impl Default for GammaloopCompileOptions {
    fn default() -> Self {
        Self {
            inline_asm: true,
            optimization_level: CompilationOptimizationLevel::O3,
            fast_math: true,
            unsafe_math: true,
            compiler: gpp(),
            custom: vec![],
        }
    }
}

impl GammaloopCompileOptions {
    pub(crate) fn inline_asm(&self) -> InlineASM {
        if self.inline_asm {
            InlineASM::default()
        } else {
            InlineASM::None
        }
    }

    #[allow(clippy::needless_update)]
    pub(crate) fn to_symbolica_compile_options(&self) -> CompileOptions {
        CompileOptions {
            optimization_level: self.optimization_level.into(),
            fast_math: self.fast_math,
            unsafe_math: self.unsafe_math,
            compiler: self.compiler.clone(),
            custom: self.custom.clone(),
            ..CompileOptions::default()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
pub struct TropicalSubgraphTableSettings {
    pub panic_on_fail: bool,
    pub target_omega: f64,
    #[serde(default = "_default_true", skip_serializing_if = "std::ops::Not::not")] // default true
    pub disable_tropical_generation: bool,
}

impl Default for TropicalSubgraphTableSettings {
    fn default() -> Self {
        Self {
            panic_on_fail: false,
            target_omega: 1.0,
            disable_tropical_generation: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default, PartialEq, JsonSchema)]
#[trait_decode(trait = GammaLoopContext)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
pub struct OrientationPattern {
    pub pat: Option<StringSerializedAtom>,
}

impl From<Atom> for OrientationPattern {
    fn from(value: Atom) -> Self {
        OrientationPattern {
            pat: Some(StringSerializedAtom(value)),
        }
    }
}

impl OrientationPattern {
    pub fn from_orientation<O: GraphOrientation>(orientation: &O) -> Self {
        orientation.orientation_delta().into()
    }

    pub fn select_pattern<'a>(&self, atom: impl AtomCore) -> Option<Atom> {
        Some(
            atom.replace(self.pat.as_ref()?.as_view().to_pattern())
                .with(function!(GS.selected, &self.pat.as_ref()?.0))
                .replace(function!(GS.orientation_delta, W_.a___))
                .level_range((0, Some(0)))
                .with(Atom::Zero),
        )
    }

    pub fn filter<O: GraphOrientation>(&self, orientation: &O) -> bool {
        if let Some(pat) = &self.pat {
            let a = orientation.orientation_delta();

            // println!("{a}vs{}", pat.0);
            let a = a
                .pattern_match(&pat.to_pattern(), None, None)
                .next()
                .is_some();
            // println!("{a}");
            a
        } else {
            true
        }
    }
}
