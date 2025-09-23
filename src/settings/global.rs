use bincode_trait_derive::{Decode, Encode};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    evaluate::{CompileOptions, ExportSettings, InlineASM},
    function,
};

use crate::{
    cff::expression::GraphOrientation,
    processes::EvaluatorSettings,
    utils::{
        serde_utils::{is_false, is_float, is_true, IsDefault},
        symbolica_ext::StringSerializedAtom,
        GS, W_,
    },
    uv::UVgenerationSettings,
    GammaLoopContext,
};

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[trait_decode(trait = GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct GenerationSettings {
    // Generation Time settings
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub evaluator: EvaluatorSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub orientation_pattern: OrientationPattern,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub compile: GammaloopCompileOptions,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub tropical_subgraph_table: TropicalSubgraphTableSettings,
    #[serde(skip_serializing_if = "is_true")]
    pub enable_thresholds: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub uv: UVgenerationSettings,
}

impl Default for GenerationSettings {
    fn default() -> Self {
        Self {
            evaluator: EvaluatorSettings::default(),
            orientation_pattern: OrientationPattern::default(),
            compile: GammaloopCompileOptions::default(),
            tropical_subgraph_table: TropicalSubgraphTableSettings::default(),
            enable_thresholds: true,
            uv: UVgenerationSettings::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Copy, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[serde(deny_unknown_fields)]
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
#[serde(default, deny_unknown_fields)]
pub struct GammaloopCompileOptions {
    #[serde(skip_serializing_if = "is_true")]
    pub inline_asm: bool,

    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub optimization_level: CompilationOptimizationLevel,

    #[serde(skip_serializing_if = "is_true")]
    pub fast_math: bool,

    #[serde(skip_serializing_if = "is_true")] // default true
    pub unsafe_math: bool,

    #[serde(skip_serializing_if = "is_gpp")] // default g++
    pub compiler: String,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
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
    pub(crate) fn export_settings(&self) -> ExportSettings {
        ExportSettings {
            inline_asm: if self.inline_asm {
                InlineASM::default()
            } else {
                InlineASM::None
            },
            ..Default::default()
        }
    }

    #[allow(clippy::needless_update)]
    pub fn to_symbolica_compile_options(&self) -> CompileOptions {
        CompileOptions {
            optimization_level: self.optimization_level.into(),
            fast_math: self.fast_math,
            unsafe_math: self.unsafe_math,
            compiler: self.compiler.clone(),
            // custom: self.custom.clone(),
            ..CompileOptions::default()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[serde(default, deny_unknown_fields)]
pub struct TropicalSubgraphTableSettings {
    #[serde(skip_serializing_if = "is_false")]
    pub panic_on_fail: bool,
    #[serde(skip_serializing_if = "is_float::<1>")] // default 1.0
    pub target_omega: f64,
    #[serde(skip_serializing_if = "is_false")]
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
#[serde(default, deny_unknown_fields)]
pub struct OrientationPattern {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
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
