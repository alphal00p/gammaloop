use spenso::{
    s,
    structure::representation::{Minkowski, RepName, Representation},
};

use crate::representations::{Bispinor, ColorAdjoint, ColorFundamental, initialize};

pub(crate) struct TestReps {
    pub mink4: Representation<Minkowski>,
    pub mink_d: Representation<Minkowski>,
    pub bis4: Representation<Bispinor>,
    pub bis_d: Representation<Bispinor>,
    pub cof_nc: Representation<ColorFundamental>,
    pub coad_na: Representation<ColorAdjoint>,
}

impl TestReps {
    pub(crate) fn new() -> Self {
        initialize();
        Self {
            mink4: Minkowski {}.new_rep(4),
            mink_d: Minkowski {}.new_rep(s!(d)),
            bis4: Bispinor {}.new_rep(4),
            bis_d: Bispinor {}.new_rep(s!(d)),
            cof_nc: ColorFundamental {}.new_rep(s!(Nc)),
            coad_na: ColorAdjoint {}.new_rep(s!(NA)),
        }
    }
}
