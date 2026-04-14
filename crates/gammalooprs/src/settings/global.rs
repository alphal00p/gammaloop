use std::fmt;

use bincode_trait_derive::{Decode, Encode};
use eyre::{Result as EyreResult, eyre};
use schemars::JsonSchema;
use serde::{Deserialize, Deserializer, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    evaluate::{CompileOptions, ExportSettings, InlineASM},
    function, try_parse,
};

use crate::{
    GammaLoopContext,
    cff::expression::GraphOrientation,
    processes::EvaluatorSettings,
    utils::{
        GS, W_,
        serde_utils::{IsDefault, is_false, is_float, is_true, is_usize, show_defaults_helper},
        symbolica_ext::StringSerializedAtom,
    },
    uv::UVgenerationSettings,
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[trait_decode(trait = GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
#[derive(Default)]
pub struct GenerationSettings {
    // Generation Time settings
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub evaluator: EvaluatorSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub feyngen: FeyGenSettings,

    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub orientation_pattern: OrientationPattern,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub compile: GammaloopCompileOptions,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub tropical_subgraph_table: TropicalSubgraphTableSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub threshold_subtraction: ThresholdSubtractionSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub vector_polarization_sum_gauge: VectorPolarizationSumGauge,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub uv: UVgenerationSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub force_cuts: Vec<Vec<String>>,
    #[serde(skip_serializing_if = "is_false")]
    pub override_lmb_heuristics: bool,
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[trait_decode(trait = GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct ThresholdSubtractionSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub enable_thresholds: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub check_esurface_at_generation: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub skip_thresholds_that_are_cuts: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub disable_integrated_ct: bool,
}

impl Default for ThresholdSubtractionSettings {
    fn default() -> Self {
        Self {
            enable_thresholds: true,
            check_esurface_at_generation: false,
            skip_thresholds_that_are_cuts: true,
            disable_integrated_ct: false,
        }
    }
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema, Default,
)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum VectorPolarizationSumGauge {
    #[serde(rename = "Feynman", alias = "feynman")]
    Feynman,
    #[default]
    #[serde(rename = "LightLikeAxial", alias = "light_like_axial")]
    LightLikeAxial,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, Copy, JsonSchema)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
#[derive(Default)]
pub enum CompilationOptimizationLevel {
    O0,
    O1,
    O2,
    #[default]
    O3,
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

impl fmt::Display for CompilationOptimizationLevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CompilationOptimizationLevel::O0 => f.write_str("O0"),
            CompilationOptimizationLevel::O1 => f.write_str("O1"),
            CompilationOptimizationLevel::O2 => f.write_str("O2"),
            CompilationOptimizationLevel::O3 => f.write_str("O3"),
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
    show_defaults_helper("g++" == compiler)
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Default)]
pub enum CompilationMode {
    #[serde(rename = "c++", alias = "cpp")]
    Cpp,
    #[serde(rename = "assembly")]
    Assembly,
    #[default]
    #[serde(rename = "symjit")]
    Symjit,
}

impl fmt::Display for CompilationMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CompilationMode::Cpp => f.write_str("c++"),
            CompilationMode::Assembly => f.write_str("assembly"),
            CompilationMode::Symjit => f.write_str("symjit"),
        }
    }
}

#[derive(
    Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema, Default,
)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[serde(default, deny_unknown_fields)]
pub struct ExternalCompilationOptionsSnapshot {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub optimization_level: CompilationOptimizationLevel,
    #[serde(skip_serializing_if = "is_true")]
    pub fast_math: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub unsafe_math: bool,
    #[serde(skip_serializing_if = "is_gpp")]
    pub compiler: String,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub custom: Vec<String>,
}

impl fmt::Display for ExternalCompilationOptionsSnapshot {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} fast_math={} unsafe_math={} compiler={} custom={}",
            self.optimization_level,
            self.fast_math,
            self.unsafe_math,
            self.compiler,
            if self.custom.is_empty() {
                "[]".to_string()
            } else {
                format!("[{}]", self.custom.join(", "))
            }
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema)]
pub enum FrozenCompilationMode {
    Eager,
    Symjit,
    Cpp(ExternalCompilationOptionsSnapshot),
    Assembly(ExternalCompilationOptionsSnapshot),
}

impl FrozenCompilationMode {
    pub fn compile_enabled(&self) -> bool {
        !matches!(self, FrozenCompilationMode::Eager)
    }

    pub fn active_backend_name(&self) -> &'static str {
        match self {
            FrozenCompilationMode::Eager => "eager",
            FrozenCompilationMode::Symjit => "symjit",
            FrozenCompilationMode::Cpp(_) => "c++",
            FrozenCompilationMode::Assembly(_) => "assembly",
        }
    }

    pub fn external_options(&self) -> Option<&ExternalCompilationOptionsSnapshot> {
        match self {
            FrozenCompilationMode::Cpp(options) | FrozenCompilationMode::Assembly(options) => {
                Some(options)
            }
            FrozenCompilationMode::Eager | FrozenCompilationMode::Symjit => None,
        }
    }

    pub fn requires_external_compilation(&self) -> bool {
        matches!(
            self,
            FrozenCompilationMode::Cpp(_) | FrozenCompilationMode::Assembly(_)
        )
    }

    pub(crate) fn export_settings(&self) -> ExportSettings {
        ExportSettings {
            inline_asm: match self {
                FrozenCompilationMode::Assembly(_) => InlineASM::default(),
                FrozenCompilationMode::Cpp(_)
                | FrozenCompilationMode::Symjit
                | FrozenCompilationMode::Eager => InlineASM::None,
            },
            ..Default::default()
        }
    }

    pub fn to_symbolica_compile_options(&self) -> Option<CompileOptions> {
        let options = self.external_options()?;
        Some(CompileOptions {
            optimization_level: options.optimization_level.into(),
            fast_math: options.fast_math,
            unsafe_math: options.unsafe_math,
            compiler: options.compiler.clone(),
            args: options.custom.clone(),
            ..CompileOptions::default()
        })
    }
}

impl fmt::Display for FrozenCompilationMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FrozenCompilationMode::Eager => f.write_str("eager"),
            FrozenCompilationMode::Symjit => f.write_str("symjit"),
            FrozenCompilationMode::Cpp(options) => write!(f, "c++ ({options})"),
            FrozenCompilationMode::Assembly(options) => write!(f, "assembly ({options})"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[serde(default, deny_unknown_fields)]
pub struct GammaloopCompileOptions {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub compilation_mode: CompilationMode,

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
            compilation_mode: CompilationMode::Symjit,
            optimization_level: CompilationOptimizationLevel::O3,
            fast_math: true,
            unsafe_math: true,
            compiler: gpp(),
            custom: vec![],
        }
    }
}

impl GammaloopCompileOptions {
    pub fn external_options_snapshot(&self) -> ExternalCompilationOptionsSnapshot {
        ExternalCompilationOptionsSnapshot {
            optimization_level: self.optimization_level,
            fast_math: self.fast_math,
            unsafe_math: self.unsafe_math,
            compiler: self.compiler.clone(),
            custom: self.custom.clone(),
        }
    }

    pub fn requires_external_compilation(&self) -> bool {
        matches!(
            self.compilation_mode,
            CompilationMode::Cpp | CompilationMode::Assembly
        )
    }

    pub fn frozen_mode(&self, evaluator_settings: &EvaluatorSettings) -> FrozenCompilationMode {
        if !evaluator_settings.compile {
            return FrozenCompilationMode::Eager;
        }

        match self.compilation_mode {
            CompilationMode::Cpp => FrozenCompilationMode::Cpp(self.external_options_snapshot()),
            CompilationMode::Assembly => {
                FrozenCompilationMode::Assembly(self.external_options_snapshot())
            }
            CompilationMode::Symjit => FrozenCompilationMode::Symjit,
        }
    }

    pub fn to_symbolica_compile_options(&self) -> CompileOptions {
        CompileOptions {
            optimization_level: self.optimization_level.into(),
            fast_math: self.fast_math,
            unsafe_math: self.unsafe_math,
            compiler: self.compiler.clone(),
            args: self.custom.clone(),
            ..CompileOptions::default()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[serde(default, deny_unknown_fields)]
#[derive(Default)]
pub struct FeyGenSettings {
    #[serde(skip_serializing_if = "is_false")]
    pub gamma_simplification_closure_check: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
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
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(default, deny_unknown_fields)]
pub struct OrientationPattern {
    #[serde(
        default,
        skip_serializing_if = "IsDefault::is_default",
        deserialize_with = "deserialize_orientation_pattern_atom"
    )]
    pub pat: Option<StringSerializedAtom>,
}

fn deserialize_orientation_pattern_atom<'de, D>(
    deserializer: D,
) -> std::result::Result<Option<StringSerializedAtom>, D::Error>
where
    D: Deserializer<'de>,
{
    let raw = Option::<String>::deserialize(deserializer)?;
    raw.map(|value| OrientationPattern::parse_user_pattern(&value))
        .transpose()
        .map(|parsed| parsed.map(StringSerializedAtom))
        .map_err(serde::de::Error::custom)
}

impl From<Atom> for OrientationPattern {
    fn from(value: Atom) -> Self {
        OrientationPattern {
            pat: Some(StringSerializedAtom(value)),
        }
    }
}

impl OrientationPattern {
    fn is_orientation_delta_atom(atom: &Atom) -> bool {
        matches!(
            atom.as_view(),
            AtomView::Fun(function)
                if function.get_symbol().get_stripped_name() == "orientation_delta"
        )
    }

    fn split_top_level_args(input: &str) -> EyreResult<Vec<String>> {
        let mut args = Vec::new();
        let mut start = 0usize;
        let mut depth = 0usize;

        for (index, ch) in input.char_indices() {
            match ch {
                '(' | '[' | '{' => depth += 1,
                ')' | ']' | '}' => {
                    if depth == 0 {
                        return Err(eyre!(
                            "Unbalanced delimiter in orientation pattern: {input}"
                        ));
                    }
                    depth -= 1;
                }
                ',' if depth == 0 => {
                    let arg = input[start..index].trim();
                    if arg.is_empty() {
                        return Err(eyre!("Empty orientation-pattern entry in pattern: {input}"));
                    }
                    args.push(arg.to_string());
                    start = index + ch.len_utf8();
                }
                _ => {}
            }
        }

        if depth != 0 {
            return Err(eyre!(
                "Unbalanced delimiter in orientation pattern: {input}"
            ));
        }

        let tail = input[start..].trim();
        if tail.is_empty() {
            if args.is_empty() {
                return Err(eyre!("Orientation pattern cannot be empty"));
            }
            return Err(eyre!("Empty orientation-pattern entry in pattern: {input}"));
        }
        args.push(tail.to_string());

        Ok(args)
    }

    fn normalize_user_pattern(pattern: &str) -> EyreResult<String> {
        let trimmed = pattern.trim();
        if trimmed.is_empty() {
            return Err(eyre!("Orientation pattern cannot be empty"));
        }

        let args = if let Some(rest) = trimmed.strip_prefix("orientation_delta") {
            let rest = rest.trim();
            if !(rest.starts_with('(') && rest.ends_with(')')) {
                return Err(eyre!(
                    "orientation_delta patterns must use parentheses, got: {pattern}"
                ));
            }
            Self::split_top_level_args(&rest[1..rest.len() - 1])?
        } else if trimmed.starts_with('(') && trimmed.ends_with(')') {
            Self::split_top_level_args(&trimmed[1..trimmed.len() - 1])?
        } else {
            Self::split_top_level_args(trimmed)?
        };

        let normalized_args = args
            .into_iter()
            .map(|arg| match arg.as_str() {
                "+" | "+1" => "1".to_string(),
                "-" | "-1" => "-1".to_string(),
                _ => arg,
            })
            .collect::<Vec<_>>()
            .join(",");

        Ok(format!("orientation_delta({normalized_args})"))
    }

    pub fn parse_user_pattern(pattern: &str) -> EyreResult<Atom> {
        let trimmed = pattern.trim();
        if trimmed.is_empty() {
            return Err(eyre!("Orientation pattern cannot be empty"));
        }

        if let Ok(parsed) = try_parse!(trimmed)
            && Self::is_orientation_delta_atom(&parsed)
        {
            return Ok(parsed);
        }

        let normalized = Self::normalize_user_pattern(pattern)?;
        try_parse!(normalized.as_str())
            .map_err(|error| eyre!("Symbolica parsing error for orientation pattern: {error}"))
    }

    pub fn from_user_pattern(pattern: &str) -> EyreResult<Self> {
        Ok(Self {
            pat: Some(StringSerializedAtom(Self::parse_user_pattern(pattern)?)),
        })
    }

    pub fn from_orientation<O: GraphOrientation>(orientation: &O) -> Self {
        orientation.orientation_delta().into()
    }

    pub fn select_pattern(&self, atom: impl AtomCore) -> Option<Atom> {
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

            // println!("{a}");
            a.pattern_match(&pat.to_pattern(), None, None)
                .next()
                .is_some()
        } else {
            true
        }
    }

    pub fn alt_filter<O: GraphOrientation>(&self, orientation: &O) -> bool {
        self.filter(orientation)
    }
}

#[cfg(test)]
mod orientation_pattern_tests {
    use super::OrientationPattern;
    use linnet::half_edge::involution::{EdgeVec, Orientation};
    use symbolica::atom::AtomCore;

    fn orientation(value: i8) -> Orientation {
        match value {
            1 => Orientation::Default,
            -1 => Orientation::Reversed,
            0 => Orientation::Undirected,
            _ => panic!("invalid orientation encoding"),
        }
    }

    fn edgevec(values: impl IntoIterator<Item = i8>) -> EdgeVec<Orientation> {
        EdgeVec::from_iter(values.into_iter().map(orientation))
    }

    #[test]
    fn orientation_pattern_deserialization_supports_shorthand_and_missing_wrapper() {
        let wrapped: OrientationPattern =
            toml::from_str(r#"pat = "orientation_delta(+,-,0)""#).unwrap();
        let tuple_only: OrientationPattern = toml::from_str(r#"pat = "(+,-,0)""#).unwrap();
        let bare_args: OrientationPattern = toml::from_str(r#"pat = "+,-,0""#).unwrap();

        assert_eq!(wrapped, tuple_only);
        assert_eq!(wrapped, bare_args);
        assert_eq!(
            wrapped.pat.as_ref().unwrap().to_string(),
            "orientation_delta(1,-1,0)"
        );
        assert!(wrapped.filter(&edgevec([1, -1, 0])));
        assert!(!wrapped.filter(&edgevec([1, 1, 0])));
    }

    #[test]
    fn orientation_pattern_repeated_wildcards_enforce_identical_bindings() {
        let pattern: OrientationPattern =
            toml::from_str(r#"pat = "(+,+,-,x_,-,x_,+,0,+,y_,-,+)""#).unwrap();

        assert!(pattern.filter(&edgevec([1, 1, -1, 1, -1, 1, 1, 0, 1, 0, -1, 1])));
        assert!(pattern.filter(&edgevec([1, 1, -1, 0, -1, 0, 1, 0, 1, -1, -1, 1])));
        assert!(!pattern.filter(&edgevec([1, 1, -1, 1, -1, 0, 1, 0, 1, 0, -1, 1])));
        assert!(pattern.alt_filter(&edgevec([1, 1, -1, 1, -1, 1, 1, 0, 1, 0, -1, 1])));
    }

    #[test]
    fn orientation_pattern_accepts_canonical_namespaced_roundtrip() {
        let original = OrientationPattern::from_user_pattern("(+,-,0)").unwrap();
        let canonical = original.pat.as_ref().unwrap().0.to_canonical_string();
        let reparsed = OrientationPattern::from_user_pattern(&canonical).unwrap();

        assert_eq!(original, reparsed);
        assert!(reparsed.filter(&edgevec([1, -1, 0])));
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[serde(default, deny_unknown_fields)]
pub struct Parallelisation {
    #[serde(skip_serializing_if = "is_usize::<1>")]
    pub feyngen: usize,
    #[serde(skip_serializing_if = "is_usize::<1>")]
    pub generate: usize,
    #[serde(skip_serializing_if = "is_usize::<1>")]
    pub compile: usize,
    #[serde(skip_serializing_if = "is_usize::<1>")]
    pub integrate: usize,
}

impl Default for Parallelisation {
    fn default() -> Self {
        Self {
            feyngen: 1,
            generate: 1,
            compile: 1,
            integrate: 1,
        }
    }
}
