use std::{fs::File, path::Path};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Result, Section};
use eyre::Context;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    evaluate::{CompileOptions, InlineASM},
    function,
};

use crate::{
    cff::expression::GraphOrientation,
    numerator::NumeratorSettings,
    utils::{symbolica_ext::StringSerializedAtom, GS, W_},
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
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GenerationSettings {
    // Generation Time settings
    pub compile_cff: bool,
    pub numerator_settings: NumeratorSettings,
    pub cpe_rounds_cff: Option<usize>,
    pub orientation_pattern: OrientationPattern,
    pub compile_separate_orientations: bool,
    pub gammaloop_compile_options: GammaloopCompileOptions,
    pub tropical_subgraph_table_settings: TropicalSubgraphTableSettings,
    pub enable_thresholds: bool,
}

impl Default for GenerationSettings {
    fn default() -> Self {
        Self {
            compile_cff: true,
            numerator_settings: NumeratorSettings::default(),
            cpe_rounds_cff: Some(1),
            orientation_pattern: OrientationPattern::default(),
            compile_separate_orientations: true,
            gammaloop_compile_options: GammaloopCompileOptions::default(),
            tropical_subgraph_table_settings: TropicalSubgraphTableSettings::default(),
            enable_thresholds: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct GammaloopCompileOptions {
    pub inline_asm: bool,
    pub optimization_level: usize,
    pub fast_math: bool,
    pub unsafe_math: bool,
    pub compiler: String,
    pub custom: Vec<String>,
}

impl Default for GammaloopCompileOptions {
    fn default() -> Self {
        Self {
            inline_asm: false,
            optimization_level: 3,
            fast_math: true,
            unsafe_math: false,
            compiler: "g++".to_owned(),
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
            optimization_level: self.optimization_level,
            fast_math: self.fast_math,
            unsafe_math: self.unsafe_math,
            compiler: self.compiler.clone(),
            custom: self.custom.clone(),
            ..CompileOptions::default()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct TropicalSubgraphTableSettings {
    pub panic_on_fail: bool,
    pub target_omega: f64,
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

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
#[trait_decode(trait = GammaLoopContext)]
pub struct OrientationPattern {
    pat: Option<StringSerializedAtom>,
}

impl OrientationPattern {
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
