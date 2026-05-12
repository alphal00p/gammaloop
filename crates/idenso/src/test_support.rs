use std::sync::LazyLock;

use spenso::{
    s,
    structure::representation::{Minkowski, RepName, Representation},
};

use crate::{
    representations::{Bispinor, ColorAdjoint, ColorFundamental},
    symbol_set,
};

#[derive(Clone)]
pub(crate) struct TestReps {
    pub mink4: Representation<Minkowski>,
    pub mink_d: Representation<Minkowski>,
    pub bis4: Representation<Bispinor>,
    pub bis_d: Representation<Bispinor>,
    pub cof_nc: Representation<ColorFundamental>,
    pub coad_na: Representation<ColorAdjoint>,
}

// Generate TestSymbols with all alphabet characters and some multi-character symbols
symbol_set!(TestSymbols, TS;
    a b c d e f h i j k l m n o p q r s t u v w x y z
    A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
    mul nu mu
);

symbol_set!(SpensoTestSymbols, SPENSO_TS, namespace = "spenso";
    l_1 l_2 l_3 l_4 l_5 l_6 l_7 l_8 l_9 l_10 l_20
    EMRID G dummy_ss l r dim
    vbar
    ebar
    edge_1_1 edge_2_1 edge_3_1 edge_4_1 edge_5_1 edge_6_1 edge_7_1 edge_8_1 edge_9_1 edge_10_1 edge_11_1 edge_12_1 edge_13_1 edge_14_1 edge_15_1 edge_16_1 edge_17_1 edge_18_1 edge_19_1 edge_20_1
    hedge0 hedge1 hedge2 hedge3 hedge4 hedge5 hedge6 hedge7 hedge8 hedge9 hedge10 hedge11 hedge12 hedge13 hedge14 hedge15 hedge16 hedge17 hedge18 hedge19 hedge20
    hedge_0 hedge_1 hedge_2 hedge_3 hedge_4 hedge_5 hedge_6 hedge_7 hedge_8 hedge_9 hedge_10 hedge_11 hedge_12 hedge_13 hedge_14 hedge_15 hedge_16 hedge_17 hedge_18 hedge_19 hedge_20
);

pub(crate) static TEST_REPS: LazyLock<TestReps> = LazyLock::new(|| {
    TestReps::initialize_symbols();
    TestReps::build()
});

pub(crate) fn test_initialize() -> &'static TestReps {
    &TEST_REPS
}

impl TestReps {
    pub(crate) fn new() -> Self {
        test_initialize().clone()
    }

    fn build() -> Self {
        Self {
            mink4: Minkowski {}.new_rep(4),
            mink_d: Minkowski {}.new_rep(s!(d)),
            bis4: Bispinor {}.new_rep(4),
            bis_d: Bispinor {}.new_rep(s!(d)),
            cof_nc: ColorFundamental {}.new_rep(s!(Nc)),
            coad_na: ColorAdjoint {}.new_rep(s!(NA)),
        }
    }

    fn initialize_symbols() {
        crate::representations::initialize();

        let tags = &spenso::network::tags::SPENSO_TAG;
        let _ = tags.rank_one_tensor_symbol("P");
        let _ = tags.rank_one_tensor_symbol("Q");
        let _ = tags.rank_one_tensor_symbol("K");

        let _ = spenso::p!();
        let _ = spenso::vector_symbol!(q);
        let _ = spenso::vector_symbol!(P);
        let _ = spenso::vector_symbol!(Q);
        let _ = spenso::vector_symbol!(K);

        let _ = TS.A;
        let _ = SPENSO_TS.G;
    }
}
