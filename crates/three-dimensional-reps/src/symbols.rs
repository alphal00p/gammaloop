use std::{collections::HashSet, sync::LazyLock};

use linnet::half_edge::involution::{EdgeIndex, Orientation};
use symbolica::{
    atom::{Atom, AtomOrView, FunctionBuilder, Symbol},
    function, symbol,
};

use crate::expression::GraphOrientation;

pub const NAMESPACE: &str = "three_dimensional_reps";
pub const GAMMALOOP_NAMESPACE: &str = "gammalooprs";
pub const SPENSO_NAMESPACE: &str = "spenso";

#[derive(Debug)]
pub struct ThreeDimensionalRepSymbols {
    pub sign: Symbol,
    pub theta: Symbol,
    pub orientation_delta: Symbol,
    pub ose: Symbol,
    pub cut_energy: Symbol,
    pub emr_momentum: Symbol,
    pub concrete_index: Symbol,
    pub numerator_sampling_scale: Symbol,
    pub affine_parameter: Symbol,
    pub coefficient: Symbol,
    pub loop_energy: Symbol,
}

impl ThreeDimensionalRepSymbols {
    pub fn sign(&self, edge: EdgeIndex) -> Atom {
        function!(self.sign, usize::from(edge) as i64)
    }

    pub fn theta<'a>(&self, arg: impl Into<AtomOrView<'a>>) -> Atom {
        let arg = arg.into();
        function!(self.theta, arg.as_view())
    }

    pub fn orientation_delta<O: GraphOrientation>(&self, orientation: &O) -> Atom {
        let args: Vec<i32> = orientation
            .orientation()
            .iter()
            .map(|(_, orientation)| match orientation {
                Orientation::Default => 1,
                Orientation::Reversed => -1,
                Orientation::Undirected => 0,
            })
            .collect();
        FunctionBuilder::new(self.orientation_delta)
            .add_args(&args)
            .finish()
    }

    pub fn ose(&self, edge: EdgeIndex) -> Atom {
        function!(self.ose, usize::from(edge) as i64)
    }

    pub fn cut_energy(&self, edge: EdgeIndex) -> Atom {
        function!(self.cut_energy, usize::from(edge) as i64)
    }

    pub fn concrete_index(&self, index: usize) -> Atom {
        function!(self.concrete_index, Atom::num(index as i64))
    }

    pub fn emr_energy(&self, edge: EdgeIndex) -> Atom {
        function!(
            self.emr_momentum,
            usize::from(edge) as i64,
            self.concrete_index(0)
        )
    }

    pub fn numerator_sampling_scale(&self) -> Atom {
        Atom::var(self.numerator_sampling_scale)
    }
}

pub static S: LazyLock<ThreeDimensionalRepSymbols> = LazyLock::new(|| ThreeDimensionalRepSymbols {
    sign: symbol!("gammalooprs::σ"),
    theta: symbol!("gammalooprs::θ"),
    orientation_delta: symbol!("gammalooprs::orientation_delta"),
    ose: symbol!("gammalooprs::OSE"),
    cut_energy: symbol!("gammalooprs::E"),
    emr_momentum: symbol!("gammalooprs::Q"),
    concrete_index: symbol!("spenso::cind"),
    numerator_sampling_scale: symbol!("gammalooprs::M"),
    affine_parameter: symbol!("three_dimensional_reps::a"),
    coefficient: symbol!("three_dimensional_reps::c"),
    loop_energy: symbol!("three_dimensional_reps::ell0"),
});

pub static SYMBOL_REGISTRY: LazyLock<HashSet<Symbol>> = LazyLock::new(|| {
    let s = &*S;
    [
        s.sign,
        s.theta,
        s.orientation_delta,
        s.ose,
        s.cut_energy,
        s.emr_momentum,
        s.concrete_index,
        s.numerator_sampling_scale,
        s.affine_parameter,
        s.coefficient,
        s.loop_energy,
    ]
    .into_iter()
    .collect()
});

pub fn sign(edge: EdgeIndex) -> Atom {
    S.sign(edge)
}

pub fn sign_theta<'a>(arg: impl Into<AtomOrView<'a>>) -> Atom {
    S.theta(arg)
}

pub fn orientation_delta<O: GraphOrientation>(orientation: &O) -> Atom {
    S.orientation_delta(orientation)
}

pub fn ose_atom_from_index(index: EdgeIndex) -> Atom {
    S.ose(index)
}

pub fn cut_energy(index: EdgeIndex) -> Atom {
    S.cut_energy(index)
}

pub fn external_energy_atom_from_index(index: EdgeIndex) -> Atom {
    S.emr_energy(index)
}

pub fn numerator_sampling_scale() -> Atom {
    S.numerator_sampling_scale()
}

pub fn sign_symbol() -> Symbol {
    S.sign
}

pub fn theta_symbol() -> Symbol {
    S.theta
}

pub fn orientation_delta_symbol() -> Symbol {
    S.orientation_delta
}

pub fn ose_symbol() -> Symbol {
    S.ose
}

pub fn energy_symbol() -> Symbol {
    S.cut_energy
}

pub fn emr_momentum_symbol() -> Symbol {
    S.emr_momentum
}

pub fn concrete_index_symbol() -> Symbol {
    S.concrete_index
}

pub fn numerator_sampling_scale_symbol() -> Symbol {
    S.numerator_sampling_scale
}
