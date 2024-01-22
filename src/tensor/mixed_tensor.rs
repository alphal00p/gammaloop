use num::Complex;
use std::{collections::BTreeMap, process::Output};
use symbolica::{
    representations::{
        default::Linear, number::Number, AsAtomView, Atom, FunctionBuilder, Identifier,
    },
    state::{ResettableBuffer, State, Workspace},
};

use super::{
    DenseTensor, Expr, HasTensorStructure, SmallestUpgrade, SparseTensor, TensorStructure,
    VecSlotExtension,
};

pub trait SmallestUpgradeSymbolic<T> {
    type LCMS<'a>;
    fn upgrade_sym<'a>(self, ws: &'a Workspace, state: &'a State) -> Option<Self::LCMS<'a>>;
}

// impl<'a> SmallestUpgrade<f64> for Expr<'a> {
//     type LCM = Expr<'a>;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// }

// impl<'a> SmallestUpgrade<Complex<f64>> for Expr<'a> {
//     type LCM = Expr<'a>;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// }

// impl SmallestUpgrade<f64> for Atom {
//     type LCM = Atom;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// }

// impl SmallestUpgrade<Complex<f64>> for Atom {
//     type LCM = Atom;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// }

impl<T, U> SmallestUpgradeSymbolic<T> for U
where
    U: SmallestUpgrade<T>,
    T: SmallestUpgrade<U, LCM = U::LCM>,
{
    type LCMS<'a> = U::LCM;
    fn upgrade_sym<'a>(self, _ws: &'a Workspace, _state: &'a State) -> Option<Self::LCMS<'a>> {
        Some(self.upgrade())
    }
}

// impl<T> SmallestUpgradeSymbolic<T> for T {
//     type LCM = T;
//     fn upgrade_sym(self, ws: &Workspace, state: &State) -> Option<Self::LCM> {
//         Some(self)
//     }
// }

// impl SmallestUpgradeSymbolic<f64> for Complex<f64> {
//     type LCM = Complex<f64>;
//     fn upgrade_sym(self, ws: &Workspace, state: &State) -> Option<Self::LCM> {
//         Some(self)
//     }
// }

// impl SmallestUpgradeSymbolic<Complex<f64>> for f64 {
//     type LCM = Complex<f64>;
//     fn upgrade_sym(self, ws: &Workspace, state: &State) -> Option<Self::LCM> {
//         Some(Complex::new(self, 0.0))
//     }
// }

// impl<T> SmallestUpgradeSymbolic<T> for Atom {
//     type LCM = Atom;
//     fn upgrade_sym(self, ws: &Workspace, state: &State) -> Option<Self::LCM> {
//         Some(self)
//     }
// }

impl SmallestUpgradeSymbolic<Atom> for f64 {
    type LCMS<'a> = Atom;
    fn upgrade_sym<'a>(self, ws: &'a Workspace, _state: &'a State) -> Option<Self::LCMS<'a>> {
        let rugrat = rug::Rational::from_f64(self)?;
        let natrat = symbolica::rings::rational::Rational::from_large(rugrat);
        let symrat = Atom::new_from_view(&ws.new_num(Number::from(natrat)).as_view());

        Some(symrat)
    }
}

impl SmallestUpgradeSymbolic<Atom> for &f64 {
    type LCMS<'a> = Atom;
    fn upgrade_sym<'a>(self, ws: &'a Workspace, _state: &'a State) -> Option<Self::LCMS<'a>> {
        let rugrat = rug::Rational::from_f64(*self)?;
        let natrat = symbolica::rings::rational::Rational::from_large(rugrat);
        let symrat = Atom::new_from_view(&ws.new_num(Number::from(natrat)).as_view());

        Some(symrat)
    }
}

impl SmallestUpgradeSymbolic<Atom> for num::Complex<f64> {
    type LCMS<'a> = Atom;
    fn upgrade_sym<'a>(self, ws: &'a Workspace, state: &'a State) -> Option<Self::LCMS<'a>> {
        return Some(Atom::new_from_view(
            &SmallestUpgradeSymbolic::<Expr>::upgrade_sym(self, ws, state)?.as_atom_view(),
        ));
    }
}

impl SmallestUpgradeSymbolic<Atom> for &num::Complex<f64> {
    type LCMS<'a> = Atom;
    fn upgrade_sym<'a>(self, ws: &'a Workspace, state: &'a State) -> Option<Self::LCMS<'a>> {
        return Some(Atom::new_from_view(
            &SmallestUpgradeSymbolic::<Expr>::upgrade_sym(*self, ws, state)?.as_atom_view(),
        ));
    }
}

impl<'a> SmallestUpgradeSymbolic<Atom> for Expr<'a> {
    type LCMS<'b> = Expr<'a>;
    fn upgrade_sym<'b>(self, _ws: &'b Workspace, _state: &'b State) -> Option<Self::LCMS<'b>> {
        Some(self)
    }
}

// impl<'a, T> SmallestUpgradeSymbolic<T> for Expr<'a> {
//     type LCM = Expr<'a>;
//     fn upgrade_sym(self, ws: &Workspace, state: &State) -> Option<Self::LCM> {
//         Some(self)
//     }
// }

impl<'a> SmallestUpgradeSymbolic<Expr<'a>> for f64 {
    type LCMS<'b> = Expr<'b>;
    fn upgrade_sym<'b>(self, ws: &'b Workspace, state: &'b State) -> Option<Self::LCMS<'b>> {
        let a: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state).unwrap();
        Some(a.builder(state, ws))
    }
}

impl<'a> SmallestUpgradeSymbolic<Expr<'a>> for &f64 {
    type LCMS<'b> = Expr<'b>;
    fn upgrade_sym<'b>(self, ws: &'b Workspace, state: &'b State) -> Option<Self::LCMS<'b>> {
        let a: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state).unwrap();
        Some(a.builder(state, ws))
    }
}

impl<'a> SmallestUpgradeSymbolic<Expr<'a>> for num::Complex<f64> {
    type LCMS<'b> = Expr<'b>;
    fn upgrade_sym<'c>(self, ws: &'c Workspace, state: &'c State) -> Option<Self::LCMS<'c>> {
        let real = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self.re, ws, state)?;
        let imag = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self.im, ws, state)?;
        let i = Atom::new_var(State::I);
        let symrat = (i.builder(state, ws) * &imag) + &real;

        Some(symrat)
    }
}

impl<'a> SmallestUpgradeSymbolic<Expr<'a>> for Atom {
    type LCMS<'b> = Atom;
    fn upgrade_sym<'c>(self, ws: &'c Workspace, state: &'c State) -> Option<Self::LCMS<'c>> {
        Some(self)
    }
}

// trait SymbolicMul<T> {
//     type Output<'a>;
//     fn mul_sym<'a>(self, other: T, ws: &'a Workspace, state: &'a State)
//         -> Option<Self::Output<'a>>;
// }

// impl<T> SymbolicMul<T> for T
// where
//     T: std::ops::Mul<T>,
// {
//     type Output<'b> = T::Output;
//     fn mul_sym<'b>(
//         self,
//         other: T,
//         _ws: &'b Workspace,
//         _state: &'b State,
//     ) -> Option<Self::Output<'b>> {
//         Some(self * other)
//     }
// }

// impl<'a> SymbolicMul<f64> for Expr<'a> {
//     type Output<'b> = Expr<'a>;
//     fn mul_sym<'b>(
//         self,
//         other: f64,
//         ws: &'b Workspace,
//         state: &'b State,
//     ) -> Option<Self::Output<'b>> {
//         let other: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(other, ws, state)?;
//         Some(self * &other)
//     }
// }

// impl<'a> SymbolicMul<Expr<'a>> for f64 {
//     type Output<'b> = Expr<'b>;
//     fn mul_sym<'b>(
//         self,
//         other: Expr<'b>,
//         ws: &'b Workspace,
//         state: &'b State,
//     ) -> Option<Self::Output<'b>> {
//         let atomself: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state)?;
//         Some(other * &atomself)
//     }
// }

// impl SymbolicMul<Atom> for f64 {
//     type Output<'a> = Atom;
//     fn mul_sym<'a>(
//         self,
//         other: Atom,
//         ws: &'a Workspace,
//         state: &'a State,
//     ) -> Option<Self::Output<'a>> {
//         let atomself: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state)?;
//         let mut out: Atom<Linear> = Atom::new();
//         other.mul(state, ws, &atomself, &mut out);
//         Some(out)
//     }
// }

// impl SymbolicMul<f64> for Atom {
//     type Output<'a> = Atom;
//     fn mul_sym<'a>(
//         self,
//         other: f64,
//         ws: &'a Workspace,
//         state: &'a State,
//     ) -> Option<Self::Output<'a>> {
//         let otheratom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(other, ws, state)?;
//         let mut out = Atom::new();
//         self.mul(state, ws, &otheratom, &mut out);
//         Some(out)
//     }
// }

trait SymbolicMul<T> {
    type Output;
    fn mul_sym<'b>(&self, other: &T, ws: &'b Workspace, state: &'b State) -> Option<Self::Output>;
}

impl<T, Out> SymbolicMul<T> for T
where
    for<'a> &'a T: std::ops::Mul<&'a T, Output = Out>,
{
    type Output = Out;
    fn mul_sym<'c>(
        &self,
        other: &T,
        _ws: &'c Workspace,
        _state: &'c State,
    ) -> Option<Self::Output> {
        Some(self * other)
    }
}

impl SymbolicMul<Atom> for f64 {
    type Output = Atom;
    fn mul_sym<'a>(
        &self,
        other: &Atom,
        ws: &'a Workspace,
        state: &'a State,
    ) -> Option<Self::Output> {
        let atomself: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state)?;
        let mut out: Atom<Linear> = Atom::new();
        other.mul(state, ws, &atomself, &mut out);
        Some(out)
    }
}

impl SymbolicMul<f64> for Atom {
    type Output = Atom;
    fn mul_sym<'a>(
        &self,
        other: &f64,
        ws: &'a Workspace,
        state: &'a State,
    ) -> Option<Self::Output> {
        let otheratom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(other, ws, state)?;
        let mut out = Atom::new();
        self.mul(state, ws, &otheratom, &mut out);
        Some(out)
    }
}

trait SymbolicAdd<T> {
    type Output;
    fn add_sym<'b>(&self, other: &T, ws: &'b Workspace, state: &'b State) -> Option<Self::Output>;
}

impl<T, Out> SymbolicAdd<T> for T
where
    for<'a> &'a T: std::ops::Add<&'a T, Output = Out>,
{
    type Output = Out;
    fn add_sym<'c>(
        &self,
        other: &T,
        _ws: &'c Workspace,
        _state: &'c State,
    ) -> Option<Self::Output> {
        Some(self + other)
    }
}

impl SymbolicAdd<Atom> for f64 {
    type Output = Atom;
    fn add_sym<'a>(
        &self,
        other: &Atom,
        ws: &'a Workspace,
        state: &'a State,
    ) -> Option<Self::Output> {
        let atomself: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state)?;
        let mut out: Atom<Linear> = Atom::new();
        other.add(state, ws, &atomself, &mut out);
        Some(out)
    }
}

impl SymbolicAdd<f64> for Atom {
    type Output = Atom;
    fn add_sym<'a>(
        &self,
        other: &f64,
        ws: &'a Workspace,
        state: &'a State,
    ) -> Option<Self::Output> {
        let otheratom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(other, ws, state)?;
        let mut out = Atom::new();
        self.add(state, ws, &otheratom, &mut out);
        Some(out)
    }
}

trait SymbolicSub<T> {
    type Output;
    fn sub_sym<'b>(&self, other: &T, ws: &'b Workspace, state: &'b State) -> Option<Self::Output>;
}

impl<T, Out> SymbolicSub<T> for T
where
    for<'a> &'a T: std::ops::Sub<&'a T, Output = Out>,
{
    type Output = Out;
    fn sub_sym<'c>(
        &self,
        other: &T,
        _ws: &'c Workspace,
        _state: &'c State,
    ) -> Option<Self::Output> {
        Some(self - other)
    }
}

impl SymbolicSub<Atom> for f64 {
    type Output = Atom;
    fn sub_sym<'a>(
        &self,
        other: &Atom,
        ws: &'a Workspace,
        state: &'a State,
    ) -> Option<Self::Output> {
        let atomself: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state)?;
        let mut negother: Atom<Linear> = Atom::new();
        other.neg(state, ws, &mut negother);
        let mut out = Atom::new();
        atomself.add(state, ws, &negother, &mut out);
        Some(out)
    }
}

impl SymbolicSub<f64> for Atom {
    type Output = Atom;
    fn sub_sym<'a>(
        &self,
        other: &f64,
        ws: &'a Workspace,
        state: &'a State,
    ) -> Option<Self::Output> {
        let otheratom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(other, ws, state)?;
        let mut negother: Atom<Linear> = Atom::new();
        otheratom.neg(state, ws, &mut negother);
        let mut out = Atom::new();
        self.add(state, ws, &negother, &mut out);
        Some(out)
    }
}
// trait SymbolicSub<T> {
//     type Output;
//     fn sub_sym<'a>(self, other: T, ws: &'a Workspace, state: &'a State)
//         -> Option<Self::Output<'a>>;
// }

// impl<T> SymbolicSub<T> for T
// where
//     T: std::ops::Sub<T>,
// {
//     type Output<'b> = T::Output;
//     fn sub_sym<'b>(
//         self,
//         other: T,
//         _ws: &'b Workspace,
//         _state: &'b State,
//     ) -> Option<Self::Output<'b>> {
//         Some(self - other)
//     }
// }

// impl SymbolicSub<Atom> for f64 {
//     type Output<'a> = Atom;
//     fn sub_sym<'a>(
//         self,
//         other: Atom,
//         ws: &'a Workspace,
//         state: &'a State,
//     ) -> Option<Self::Output<'a>> {
//         let atomself: Atom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(self, ws, state)?;
//         let mut negother: Atom<Linear> = Atom::new();
//         other.neg(state, ws, &mut negother);
//         let mut out = Atom::new();
//         atomself.add(state, ws, &negother, &mut out);
//         Some(out)
//     }
// }

// impl SymbolicSub<f64> for Atom {
//     type Output<'a> = Atom;
//     fn sub_sym<'a>(
//         self,
//         other: f64,
//         ws: &'a Workspace,
//         state: &'a State,
//     ) -> Option<Self::Output<'a>> {
//         let otheratom = SmallestUpgradeSymbolic::<Atom>::upgrade_sym(other, ws, state)?;
//         let mut negother: Atom<Linear> = Atom::new();
//         otheratom.neg(state, ws, &mut negother);
//         let mut out = Atom::new();
//         self.add(state, ws, &negother, &mut out);
//         Some(out)
//     }
// }

// pub fn mul_sym<'a, T, U>(
//     right: T,
//     left: U,
//     ws: &'a Workspace,
//     state: &'a State,
// ) -> Option<U::LCMS<'a>>
// where
//     T: SmallestUpgradeSymbolic<U>,
//     U: SmallestUpgradeSymbolic<T>,
//     U::LCMS<'a>: std::ops::Mul<T::LCMS<'a>, Output = U::LCMS<'a>>,
// {
//     let right = right.upgrade_sym(ws, state)?;
//     let left = left.upgrade_sym(ws, state)?;
//     Some(left * right)
// }

// impl<'a> SmallestUpgrade<Expr<'a>> for f64 {
//     type Output = Expr<'a>;
//     fn upgrade(self) -> Self::Output {
//         Atom::new_num(self)
//     }
// }

trait SymbolicZero {
    fn zero(state: &State, ws: &Workspace) -> Self;
}

// impl<T> SymbolicZero for T
// where
//     T: num::traits::Zero,
// {
//     fn zero(_state: &State, _ws: &Workspace) -> Self {
//         Self::zero()
//     }
// }

impl SymbolicZero for f64 {
    fn zero(_state: &State, _ws: &Workspace) -> Self {
        0.
    }
}

impl SymbolicZero for Atom {
    fn zero(_state: &State, _ws: &Workspace) -> Self {
        Atom::new_num(0)
    }
}

pub trait SymbolicContract<T> {
    type LCM;
    fn contract_sym(&self, other: &T, ws: &Workspace, state: &State) -> Option<Self::LCM>;
}

impl<T, U, Out> SymbolicContract<DenseTensor<T>> for DenseTensor<U>
where
    U: SymbolicAdd<T, Output = Out> + SymbolicMul<T, Output = Out> + SymbolicSub<T, Output = Out>,
    Out: SymbolicZero
        + Clone
        + SymbolicAdd<Out, Output = Out>
        + SymbolicSub<Out, Output = Out>
        + std::fmt::Debug
        + for<'a> std::ops::AddAssign<&'a Out>
        + std::ops::Neg<Output = Out>
        + for<'b> std::ops::SubAssign<&'b Out>,
{
    type LCM = DenseTensor<Out>;
    fn contract_sym(
        &self,
        other: &DenseTensor<T>,
        ws: &Workspace,
        state: &State,
    ) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            let dimension_of_contraction = self.shape()[i];
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let metric = self.structure()[i].representation.negative();
            for (index_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();
                    for i in 0..dimension_of_contraction {
                        if metric[i] {
                            result_data[result_index] = result_data[result_index].sub_sym(
                                &fiber_a[i].mul_sym(fiber_b[i], ws, state)?,
                                ws,
                                state,
                            )?;
                        } else {
                            result_data[result_index] = result_data[result_index].add_sym(
                                &fiber_a[i].mul_sym(fiber_b[i], ws, state)?,
                                ws,
                                state,
                            )?;
                        }
                    }
                }
            }

            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<T, U, Out> SymbolicContract<DenseTensor<T>> for SparseTensor<U>
where
    U: SymbolicAdd<T, Output = Out> + SymbolicMul<T, Output = Out> + SymbolicSub<T, Output = Out>,
    Out: SymbolicZero
        + Clone
        + SymbolicAdd<Out, Output = Out>
        + SymbolicSub<Out, Output = Out>
        + std::fmt::Debug
        + for<'a> std::ops::AddAssign<&'a Out>
        + std::ops::Neg<Output = Out>
        + for<'b> std::ops::SubAssign<&'b Out>,
{
    type LCM = DenseTensor<Out>;
    fn contract_sym(
        &self,
        other: &DenseTensor<T>,
        ws: &Workspace,
        state: &State,
    ) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let metric = self.structure()[i].representation.negative();
            for (index_a, nonzeros, fiber_a) in self.iter_fibers(i) {
                for (index_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();
                    for (i, k) in nonzeros.iter().enumerate() {
                        if metric[*k] {
                            result_data[result_index] = result_data[result_index].sub_sym(
                                &fiber_a[i].mul_sym(fiber_b[*k], ws, state)?,
                                ws,
                                state,
                            )?;
                        } else {
                            result_data[result_index] = result_data[result_index].add_sym(
                                &fiber_a[i].mul_sym(fiber_b[*k], ws, state)?,
                                ws,
                                state,
                            )?;
                        }
                    }
                }
            }

            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<T, U, Out> SymbolicContract<SparseTensor<T>> for SparseTensor<U>
where
    U: SymbolicAdd<T, Output = Out> + SymbolicMul<T, Output = Out> + SymbolicSub<T, Output = Out>,
    Out: SymbolicZero
        + Clone
        + SymbolicAdd<Out, Output = Out>
        + SymbolicSub<Out, Output = Out>
        + std::fmt::Debug
        + for<'a> std::ops::AddAssign<&'a Out>
        + std::ops::Neg<Output = Out>
        + for<'b> std::ops::SubAssign<&'b Out>,
{
    type LCM = SparseTensor<Out>;
    fn contract_sym(
        &self,
        other: &SparseTensor<T>,
        ws: &Workspace,
        state: &State,
    ) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = BTreeMap::new();
            let metric = self.structure()[i].representation.negative();
            for (index_a, nonzeros_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, nonzeros_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = index_a[..i]
                        .iter()
                        .chain(&index_a[i + 1..])
                        .chain(&index_b[..j])
                        .chain(&index_b[j + 1..])
                        .cloned()
                        .collect::<Vec<_>>();
                    let mut value = Out::zero(state, ws);
                    let mut nonzero = false;
                    for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                        nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x)) // Only store the positions
                    }) {
                        nonzero = true;
                        if metric[x] {
                            value = value.sub_sym(
                                &fiber_a[i].mul_sym(fiber_b[j], ws, state)?,
                                ws,
                                state,
                            )?;
                        } else {
                            value = value.add_sym(
                                &fiber_a[i].mul_sym(fiber_b[j], ws, state)?,
                                ws,
                                state,
                            )?;
                        }
                    }

                    if nonzero
                    // && value.as_atom_view() != zero.as_atom_view() TODO impl generic
                    {
                        result_data.insert(result_index, value);
                    }
                }
            }

            let result = SparseTensor {
                elements: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<T, U, Out> SymbolicContract<SparseTensor<T>> for DenseTensor<U>
where
    T: SymbolicAdd<U, Output = Out> + SymbolicMul<U, Output = Out> + SymbolicSub<U, Output = Out>,
    Out: SymbolicZero
        + Clone
        + SymbolicAdd<Out, Output = Out>
        + SymbolicSub<Out, Output = Out>
        + std::fmt::Debug
        + for<'a> std::ops::AddAssign<&'a Out>
        + std::ops::Neg<Output = Out>
        + for<'b> std::ops::SubAssign<&'b Out>,
{
    type LCM = DenseTensor<Out>;
    fn contract_sym(
        &self,
        other: &SparseTensor<T>,
        ws: &Workspace,
        state: &State,
    ) -> Option<Self::LCM> {
        other.contract_sym(self, ws, state)
    }
}

impl SparseTensor<Atom> {
    pub fn builder<'a>(&self, state: &'a State, ws: &'a Workspace) -> SparseTensor<Expr<'a>> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(index, value.builder(state, ws)).unwrap();
        }
        result
    }
}

impl DenseTensor<Atom> {
    pub fn builder<'a>(&self, state: &'a State, ws: &'a Workspace) -> DenseTensor<Expr<'a>> {
        let mut result = DenseTensor::neutral_builder(self.structure.clone(), ws, state);
        for (index, value) in self.iter() {
            result.set(&index, value.builder(state, ws));
        }
        result
    }
}

// impl<T: ConvertableToSymbolic> SparseTensor<T> {
//     pub fn to_symbolic<'a>(&self, ws: &'a Workspace, state: &'a State) -> SparseTensor<Atom> {
//         let mut result = SparseTensor::empty(self.structure.clone());
//         for (index, value) in self.iter() {
//             result
//                 .set(index, value.to_symbolic(ws, state).unwrap())
//                 .unwrap();
//         }
//         result
//     }

//     pub fn to_symbolic_builder<'a>(
//         &self,
//         ws: &'a Workspace,
//         state: &'a State,
//     ) -> SparseTensor<Expr<'a>> {
//         let mut result = SparseTensor::empty(self.structure.clone());
//         for (index, value) in self.iter() {
//             result
//                 .set(
//                     index,
//                     value.to_symbolic(ws, state).unwrap().builder(state, ws),
//                 )
//                 .unwrap();
//         }
//         result
//     }
// }

impl<'a, T> SparseTensor<T>
where
    for<'d> &'d T: SmallestUpgradeSymbolic<Atom, LCMS<'a> = Atom>,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> SparseTensor<Atom> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(
                &index,
                SmallestUpgradeSymbolic::<Atom>::upgrade_sym(value, ws, state).unwrap(),
            );
        }
        result
    }

    pub fn to_symbolic_builder<'b>(
        &'b self,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> SparseTensor<Expr<'a>> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(
                &index,
                SmallestUpgradeSymbolic::<Atom>::upgrade_sym(value, ws, state)
                    .unwrap()
                    .builder(state, ws),
            );
        }
        result
    }
}

impl<'a> DenseTensor<Expr<'a>> {
    pub fn symbolic_zeros(structure: TensorStructure) -> DenseTensor<Atom> {
        let result_data = vec![0; structure.size()];

        DenseTensor {
            data: result_data.iter().map(|&x| Atom::new_num(x)).collect(),
            structure,
        }
    }

    pub fn neutral_builder(
        structure: TensorStructure,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'a>> {
        let zero = Atom::new_num(0).builder(state, ws);
        let result_data = vec![zero; structure.size()];
        DenseTensor {
            data: result_data,
            structure,
        }
    }
    pub fn symbolic_labels(
        label: &str,
        structure: TensorStructure,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> DenseTensor<Atom> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let indices_str = index
                .into_iter()
                .map(|index| index.to_string())
                .collect::<Vec<String>>()
                .join("_");

            let value = Atom::parse(&format!("{}_{}", label, indices_str), state, ws).unwrap();

            data.push(value);
        }
        DenseTensor { data, structure }
    }

    pub fn numbered_labeled(
        number: usize,
        label: Identifier,
        structure: TensorStructure,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Atom> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let mut value_builder = FunctionBuilder::new(label, state, ws);
            value_builder = value_builder.add_arg(Atom::new_num(number as i64).as_atom_view());

            for i in index {
                value_builder = value_builder.add_arg(Atom::new_num(i as i64).as_atom_view());
            }
            // Atom::parse(&format!("{}_{}_{}", label, indices_str, i), state, ws).unwrap();

            let value = value_builder.finish().into_atom();

            data.push(value);
        }
        DenseTensor { data, structure }
    }

    pub fn numbered_labeled_builder(
        number: usize,
        label: Identifier,
        structure: TensorStructure,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'a>> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let mut value_builder = FunctionBuilder::new(label, state, ws);
            value_builder = value_builder.add_arg(Atom::new_num(number as i64).as_atom_view());

            for i in index {
                value_builder = value_builder.add_arg(Atom::new_num(i as i64).as_atom_view());
            }
            // Atom::parse(&format!("{}_{}_{}", label, indices_str, i), state, ws).unwrap();

            let value = value_builder.finish();

            data.push(value);
        }
        DenseTensor { data, structure }
    }

    pub fn finish(self) -> DenseTensor<Atom> {
        DenseTensor {
            data: self.data.into_iter().map(|x| x.into_atom()).collect(),
            structure: self.structure,
        }
    }
}

impl<'a, T> DenseTensor<T>
where
    for<'d> &'d T: SmallestUpgradeSymbolic<Atom, LCMS<'a> = Atom>,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> DenseTensor<Atom> {
        let mut result = DenseTensor::symbolic_zeros(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(
                &index,
                SmallestUpgradeSymbolic::<Atom>::upgrade_sym(value, ws, state).unwrap(),
            );
        }
        result
    }

    pub fn to_symbolic_builder<'b>(
        &'b self,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> DenseTensor<Expr<'a>> {
        let mut result = DenseTensor::neutral_builder(self.structure.clone(), ws, state);
        for (index, value) in self.iter() {
            result.set(
                &index,
                SmallestUpgradeSymbolic::<Atom>::upgrade_sym(value, ws, state)
                    .unwrap()
                    .builder(state, ws),
            );
        }
        result
    }
}
