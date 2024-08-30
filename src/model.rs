use crate::momentum::{FourMomentum, Helicity, Polarization};
use crate::utils::{self, FloatLike, F};

use ahash::{AHashMap, RandomState};
use color_eyre::{Help, Report};
use eyre::{eyre, Context};
use itertools::Itertools;

// use log::{info, trace};
use serde::{Deserialize, Serialize};
use serde_yaml::Error;
use smartstring::{LazyCompact, SmartString};
use spenso::parametric::{ExpandedCoefficent, TensorCoefficient};
use spenso::structure::{Dimension, ExpandedIndex, FlatIndex, VecStructure, CONCRETEIND};
use spenso::{
    contraction::IsZero,
    structure::{
        BaseRepName, Bispinor, ColorAdjoint, ColorFundamental, ColorSextet, Dual, DualSlotTo,
        IsAbstractSlot, Lorentz, PhysReps, RepName, Representation, Slot, ABSTRACTIND,
    },
};
use std::fs;

use eyre::Result;
use std::ops::{Index, Neg};
use std::path::Path;
use symbolica::domains::rational::Rational;
use symbolica::evaluate::FunctionMap;
use symbolica::id::Pattern;
// use std::str::pattern::Pattern;
use std::sync::Arc;
use std::{collections::HashMap, fs::File};
use symbolica::atom::{Atom, AtomView, FunctionBuilder, Symbol};

use spenso::complex::Complex;
use symbolica::domains::float::NumericalFloatLike;

use symbolica::fun;
use symbolica::printer::{AtomPrinter, PrintOptions};
use symbolica::state::State;
use utils::approx_eq_cmplx_vec;

#[allow(unused)]
fn normalise_complex(atom: &Atom) -> Atom {
    let re = Atom::parse("re_").unwrap();
    let im = Atom::parse("im_").unwrap();

    let comp_id = State::get_symbol("complex");

    let complexfn = fun!(comp_id, re, im).into_pattern();

    let i = Atom::new_var(State::I);
    let complexpanded = &re + i * &im;

    complexfn.replace_all(atom.as_view(), &complexpanded.into_pattern(), None, None)
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub enum ParameterNature {
    #[default]
    #[serde(rename = "external")]
    External,
    #[serde(rename = "internal")]
    Internal,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub enum ParameterType {
    #[default]
    #[serde(rename = "real")]
    Real,
    #[serde(rename = "complex")]
    Imaginary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableVertexRule {
    pub name: SmartString<LazyCompact>,
    pub particles: Vec<SmartString<LazyCompact>>,
    pub color_structures: Vec<SmartString<LazyCompact>>,
    pub lorentz_structures: Vec<SmartString<LazyCompact>>,
    pub couplings: Vec<Vec<Option<SmartString<LazyCompact>>>>,
}

impl SerializableVertexRule {
    pub fn from_vertex_rule(vertex_rule: &VertexRule) -> SerializableVertexRule {
        SerializableVertexRule {
            name: vertex_rule.name.clone(),
            particles: vertex_rule
                .particles
                .iter()
                .map(|particle| particle.name.clone())
                .collect(),
            color_structures: vertex_rule
                .color_structures
                .iter()
                .map(utils::to_str_expression)
                .map(SmartString::from)
                .collect(),
            lorentz_structures: vertex_rule
                .lorentz_structures
                .iter()
                .map(|lorentz_structure| lorentz_structure.name.clone())
                .collect(),
            couplings: vertex_rule
                .couplings
                .iter()
                .map(|couplings| {
                    couplings
                        .iter()
                        .map(|coupling| coupling.as_ref().map(|cpl| cpl.name.clone()))
                        .collect()
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ColorStructure {
    color_structure: Vec<Atom>,
}

impl ColorStructure {
    pub fn new(color_structure: Vec<Atom>) -> Self {
        ColorStructure { color_structure }
    }

    pub fn iter(&self) -> std::slice::Iter<Atom> {
        self.color_structure.iter()
    }
}

impl FromIterator<Atom> for ColorStructure {
    fn from_iter<T: IntoIterator<Item = Atom>>(iter: T) -> Self {
        ColorStructure {
            color_structure: iter.into_iter().collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct VertexRule {
    pub name: SmartString<LazyCompact>,
    pub particles: Vec<Arc<Particle>>,
    pub color_structures: ColorStructure,
    pub lorentz_structures: Vec<Arc<LorentzStructure>>,
    pub couplings: Vec<Vec<Option<Arc<Coupling>>>>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VertexSlots {
    edge_slots: Vec<EdgeSlots<Dual<Lorentz>>>,
}

impl Index<usize> for VertexSlots {
    type Output = EdgeSlots<Dual<Lorentz>>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.edge_slots[index]
    }
}

impl From<EdgeSlots<Dual<Lorentz>>> for VertexSlots {
    fn from(value: EdgeSlots<Dual<Lorentz>>) -> Self {
        VertexSlots {
            edge_slots: vec![value],
        }
    }
}

impl VertexRule {
    pub fn generate_vertex_slots(
        &self,
        mut shifts: (usize, usize, usize),
    ) -> (VertexSlots, (usize, usize, usize)) {
        let mut edge_slots = vec![];
        for p in &self.particles {
            let (e, s) = p.slots(shifts);
            edge_slots.push(e);
            shifts = s;
        }
        (VertexSlots { edge_slots }, shifts)
    }
    pub fn from_serializable_vertex_rule(
        model: &Model,
        vertex_rule: &SerializableVertexRule,
    ) -> VertexRule {
        VertexRule {
            name: vertex_rule.name.clone(),
            particles: vertex_rule
                .particles
                .iter()
                .map(|particle_name| model.get_particle(particle_name).clone())
                .collect(),
            color_structures: vertex_rule
                .color_structures
                .iter()
                .map(|color_structure_name| {
                    utils::parse_python_expression(color_structure_name.as_str())
                })
                .collect(),
            lorentz_structures: vertex_rule
                .lorentz_structures
                .iter()
                .map(|lorentz_structure_name| {
                    model.get_lorentz_structure(lorentz_structure_name).clone()
                })
                .collect(),
            couplings: vertex_rule
                .couplings
                .iter()
                .map(|coupling_names| {
                    coupling_names
                        .iter()
                        .map(|coupling_name| {
                            coupling_name
                                .as_ref()
                                .map(|cpl_name| model.get_coupling(cpl_name))
                        })
                        .collect()
                })
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializablePropagator {
    pub name: SmartString<LazyCompact>,
    pub particle: SmartString<LazyCompact>,
    pub numerator: SmartString<LazyCompact>,
    pub denominator: SmartString<LazyCompact>,
}

impl SerializablePropagator {
    pub fn from_propagator(propagator: &Propagator) -> SerializablePropagator {
        SerializablePropagator {
            name: propagator.name.clone(),
            particle: propagator.particle.name.clone(),
            numerator: utils::to_str_expression(&propagator.numerator).into(),
            denominator: utils::to_str_expression(&propagator.denominator).into(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Propagator {
    pub name: SmartString<LazyCompact>,
    pub particle: Arc<Particle>,
    pub numerator: Atom,
    pub denominator: Atom,
}

impl Propagator {
    pub fn from_serializable_propagator(
        model: &Model,
        propagator: &SerializablePropagator,
    ) -> Propagator {
        Propagator {
            name: propagator.name.clone(),
            particle: model.get_particle(&propagator.particle).clone(),
            numerator: utils::parse_python_expression(propagator.numerator.as_str()),
            denominator: utils::parse_python_expression(propagator.denominator.as_str()),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableCoupling {
    name: SmartString<LazyCompact>,
    expression: SmartString<LazyCompact>,
    #[serde(with = "vectorize")]
    orders: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    value: Option<(f64, f64)>,
}

impl SerializableCoupling {
    pub fn from_coupling(coupling: &Coupling) -> SerializableCoupling {
        SerializableCoupling {
            name: coupling.name.clone(),
            expression: utils::to_str_expression(&coupling.expression).into(),
            orders: coupling.orders.clone(),
            value: coupling.value.map(|value| (value.re, value.im)),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Coupling {
    pub name: SmartString<LazyCompact>,
    pub expression: Atom,
    pub orders: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub value: Option<Complex<f64>>,
}

impl Coupling {
    pub fn from_serializable_coupling(coupling: &SerializableCoupling) -> Coupling {
        Coupling {
            name: coupling.name.clone(),
            expression: utils::parse_python_expression(coupling.expression.as_str()),
            orders: coupling.orders.clone(),
            value: coupling.value.map(|value| Complex::new(value.0, value.1)),
        }
    }

    pub fn rep_rule(&self) -> [Atom; 2] {
        let lhs = Atom::parse(&self.name).unwrap();
        //let rhs = normalise_complex(&self.expression);
        let rhs = self.expression.clone();

        [lhs, rhs]
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableParticle {
    pdg_code: isize,
    name: SmartString<LazyCompact>,
    antiname: SmartString<LazyCompact>,
    spin: isize,
    color: isize,
    mass: SmartString<LazyCompact>,
    width: SmartString<LazyCompact>,
    texname: SmartString<LazyCompact>,
    antitexname: SmartString<LazyCompact>,
    charge: f64,
    ghost_number: isize,
    lepton_number: isize,
    y_charge: isize,
}

impl SerializableParticle {
    pub fn from_particle(particle: &Particle) -> SerializableParticle {
        SerializableParticle {
            pdg_code: particle.pdg_code,
            name: particle.name.clone(),
            antiname: particle.antiname.clone(),
            spin: particle.spin,
            color: particle.color,
            mass: particle.mass.name.clone(),
            width: particle.width.name.clone(),
            texname: particle.texname.clone(),
            antitexname: particle.antitexname.clone(),
            charge: particle.charge,
            ghost_number: particle.ghost_number,
            lepton_number: particle.lepton_number,
            y_charge: particle.y_charge,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Particle {
    pub pdg_code: isize,
    pub name: SmartString<LazyCompact>,
    pub antiname: SmartString<LazyCompact>,
    pub spin: isize,
    pub color: isize,
    pub mass: Arc<Parameter>,
    pub width: Arc<Parameter>,
    pub texname: SmartString<LazyCompact>,
    pub antitexname: SmartString<LazyCompact>,
    pub charge: f64,
    pub ghost_number: isize,
    pub lepton_number: isize,
    pub y_charge: isize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InOutIndex {
    incoming: Slot<PhysReps>,
    outgoing: Slot<PhysReps>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeSlots<LorRep: RepName> {
    lorentz: Vec<Slot<LorRep>>,
    spin: Vec<Slot<Bispinor>>,
    pub color: Vec<Slot<PhysReps>>,
}

impl From<EdgeSlots<Lorentz>> for VecStructure {
    fn from(value: EdgeSlots<Lorentz>) -> Self {
        VecStructure {
            structure: value
                .lorentz
                .into_iter()
                .map(|x| x.into())
                .chain(value.spin.into_iter().map(|a| a.into()))
                .chain(value.color)
                .collect_vec(),
        }
    }
}

impl<LorRep: BaseRepName> EdgeSlots<LorRep>
where
    PhysReps: From<LorRep>,
{
    pub fn expanded_index(&self, flat_index: FlatIndex) -> Result<ExpandedIndex> {
        let mut indices = vec![];
        let mut index: usize = flat_index.into();
        for &stride in &self.strides_row_major()? {
            indices.push(index / stride);
            index %= stride;
        }
        if usize::from(flat_index) < self.size() {
            Ok(indices.into())
        } else {
            Err(eyre!("Index {flat_index} out of bounds"))
        }
    }

    fn order(&self) -> usize {
        self.color.len() + self.lorentz.len() + self.spin.len()
    }

    fn shape(&self) -> Vec<Dimension> {
        let mut dims = vec![];
        dims.extend(self.lorentz.iter().map(|s| s.dim()));
        dims.extend(self.spin.iter().map(|s| s.dim()));
        dims.extend(self.color.iter().map(|s| s.dim()));
        dims
    }

    fn strides_row_major(&self) -> Result<Vec<usize>> {
        let mut strides = vec![1; self.order()];
        if self.order() == 0 {
            return Ok(strides);
        }

        for i in (0..self.order() - 1).rev() {
            strides[i] = strides[i + 1] * usize::try_from(self.shape()[i + 1])?;
        }

        Ok(strides)
    }
    pub fn to_dense_labels<T>(&self, index_to_atom: impl Fn(&Self, FlatIndex) -> T) -> Vec<Atom>
    where
        Self: Sized,
        T: TensorCoefficient,
    {
        let mut data = vec![];
        for index in 0..self.size() {
            data.push(index_to_atom(self, index.into()).to_atom().unwrap());
        }
        data
    }
    pub fn size(&self) -> usize {
        if self.spin.is_empty() && self.lorentz.is_empty() && self.color.is_empty() {
            0
        } else {
            self.spin_size() * self.color_size() * self.lorentz_size()
        }
    }

    pub fn lorentz_size(&self) -> usize {
        self.lorentz
            .iter()
            .map(|s| usize::try_from(s.dim()).unwrap())
            .product()
    }

    pub fn spin_size(&self) -> usize {
        self.spin
            .iter()
            .map(|s| usize::try_from(s.dim()).unwrap())
            .product()
    }

    pub fn color_size(&self) -> usize {
        self.color
            .iter()
            .map(|s| usize::try_from(s.dim()).unwrap())
            .product()
    }

    pub fn dual(&self) -> EdgeSlots<LorRep::Dual> {
        EdgeSlots {
            lorentz: self.lorentz.iter().map(|l| l.dual()).collect(),
            spin: self.spin.iter().map(|s| s.dual()).collect(),
            color: self.color.iter().map(|c| c.dual()).collect(),
        }
    }
    pub fn replacements(&self, id: usize) -> Vec<(Pattern, Pattern)> {
        let rhs_lor = LorRep::new_slot_selfless(4, id)
            .to_symbolic_wrapped()
            .into_pattern();

        let rhs_spin = Bispinor::new_slot_selfless(4, id);

        let rhs_spin = rhs_spin.to_symbolic_wrapped().into_pattern();

        let mut reps = vec![];
        for l in &self.lorentz {
            reps.push((rhs_lor.clone(), l.to_symbolic().into_pattern()));
        }

        for s in &self.spin {
            reps.push((rhs_spin.clone(), s.to_symbolic().into_pattern()));
        }

        for c in &self.color {
            let mut rhs_color = *c;
            rhs_color.aind = id.into();
            let rhs_color = rhs_color.to_symbolic_wrapped().into_pattern();

            reps.push((rhs_color.clone(), c.to_symbolic().into_pattern()));
        }

        reps
    }
    pub fn to_aind_atom(&self) -> Atom {
        let mut builder = FunctionBuilder::new(State::get_symbol(ABSTRACTIND));

        for l in &self.lorentz {
            builder = builder.add_arg(l.to_symbolic().as_view());
        }

        for s in &self.spin {
            builder = builder.add_arg(s.to_symbolic().as_view());
        }

        for c in &self.color {
            builder = builder.add_arg(c.to_symbolic().as_view());
        }

        builder.finish()
    }

    pub fn to_cind_atom(&self) -> Atom {
        let mut builder = FunctionBuilder::new(State::get_symbol(CONCRETEIND));

        for l in &self.lorentz {
            builder = builder.add_arg(l.to_symbolic().as_view());
        }

        for s in &self.spin {
            builder = builder.add_arg(s.to_symbolic().as_view());
        }

        for c in &self.color {
            builder = builder.add_arg(c.to_symbolic().as_view());
        }

        builder.finish()
    }
}

impl Particle {
    fn lorentz_slots<LR: BaseRepName>(&self, shift: usize) -> (Vec<Slot<LR>>, usize) {
        let fourd_lor = LR::new_dimed_rep_selfless(4);

        match self.spin {
            3 => (vec![Representation::new_slot(&fourd_lor, shift)], shift + 1),
            _ => (vec![], shift),
        }
    }

    fn spin_slots(&self, shift: usize) -> (Vec<Slot<Bispinor>>, usize) {
        let fourd_bis: Representation<_> = Bispinor::new_dimed_rep_selfless(4);

        match self.spin {
            2 => (vec![fourd_bis.new_slot(shift)], shift + 1),
            _ => (vec![], shift),
        }
    }

    fn color_slots(&self, shift: usize) -> (Vec<Slot<PhysReps>>, usize) {
        let rep: Representation<PhysReps> = match self.color {
            3 => ColorFundamental::new_dimed_rep_selfless(3).cast(),
            -3 => Dual::<ColorFundamental>::new_dimed_rep_selfless(3).cast(),
            6 => ColorSextet::new_dimed_rep_selfless(6).cast(),
            -6 => Dual::<ColorSextet>::new_dimed_rep_selfless(6).cast(),
            8 => ColorAdjoint::new_dimed_rep_selfless(8).cast(),
            _ => return (vec![], shift),
        };

        (vec![rep.new_slot(shift)], shift + 1)
    }

    pub fn slots<LR: BaseRepName>(
        &self,
        shifts: (usize, usize, usize),
    ) -> (EdgeSlots<LR>, (usize, usize, usize)) {
        let (lorentz, shift_lor) = self.lorentz_slots(shifts.0);
        let (spin, shift_spin) = self.spin_slots(shifts.1);
        let (color, shift_color) = self.color_slots(shifts.2);

        (
            EdgeSlots {
                lorentz,
                spin,
                color,
            },
            (shift_lor, shift_spin, shift_color),
        )
    }

    pub fn from_serializable_particle(model: &Model, particle: &SerializableParticle) -> Particle {
        Particle {
            pdg_code: particle.pdg_code,
            name: particle.name.clone(),
            antiname: particle.antiname.clone(),
            spin: particle.spin,
            color: particle.color,
            mass: model.get_parameter(&particle.mass),
            width: model.get_parameter(&particle.width),
            texname: particle.texname.clone(),
            antitexname: particle.antitexname.clone(),
            charge: particle.charge,
            ghost_number: particle.ghost_number,
            lepton_number: particle.lepton_number,
            y_charge: particle.y_charge,
        }
    }

    pub fn incoming_polarization_atom(&self, edge_slots: &EdgeSlots<Lorentz>, num: usize) -> Atom {
        match self.spin {
            1 => Atom::parse("1").unwrap(),
            2 => {
                if self.pdg_code > 0 {
                    let mut u = FunctionBuilder::new(State::get_symbol("u"));
                    u = u.add_arg(&Atom::new_num(num as i64));
                    u = u.add_arg(&edge_slots.to_aind_atom());
                    u.finish()
                } else {
                    let mut vbar = FunctionBuilder::new(State::get_symbol("vbar"));
                    vbar = vbar.add_arg(&Atom::new_num(num as i64));
                    vbar = vbar.add_arg(&edge_slots.to_aind_atom());
                    vbar.finish()
                }
            }
            3 => {
                let mut e = FunctionBuilder::new(State::get_symbol("ϵ"));
                e = e.add_arg(&Atom::new_num(num as i64));
                e = e.add_arg(&edge_slots.to_aind_atom());
                e.finish()
            }
            _ => panic!("higher spin not supported"), //Atom::parse("1").unwrap(),
        }
    }

    pub fn in_pol_symbol(&self) -> Option<Symbol> {
        match self.spin {
            2 => {
                if self.pdg_code > 0 {
                    Some(State::get_symbol("u"))
                } else {
                    Some(State::get_symbol("vbar"))
                }
            }
            3 => Some(State::get_symbol("ϵ")),
            _ => None,
        }
    }

    pub fn out_pol_symbol(&self) -> Option<Symbol> {
        match self.spin {
            2 => {
                if self.pdg_code > 0 {
                    Some(State::get_symbol("ubar"))
                } else {
                    Some(State::get_symbol("v"))
                }
            }
            3 => Some(State::get_symbol("ϵbar")),
            _ => None,
        }
    }

    pub fn incoming_polarization_atom_concrete(
        &self,
        edge_slots: &EdgeSlots<Lorentz>,
        num: usize,
    ) -> Vec<Atom> {
        edge_slots.to_dense_labels(|v, i| ExpandedCoefficent::<usize> {
            index: v.expanded_index(i).unwrap(),
            name: self.in_pol_symbol(),
            args: Some(num),
        })
    }

    pub fn outgoing_polarization_atom_concrete(
        &self,
        edge_slots: &EdgeSlots<Lorentz>,
        num: usize,
    ) -> Vec<Atom> {
        edge_slots.to_dense_labels(|v, i| ExpandedCoefficent::<usize> {
            index: v.expanded_index(i).unwrap(),
            name: self.out_pol_symbol(),
            args: Some(num),
        })
    }

    pub fn incoming_polarization_match<T: FloatLike>(
        &self,
        num: usize,
        mom: &FourMomentum<F<T>>,
        helicity: Helicity,
    ) -> Vec<(Atom, Complex<F<T>>)> {
        let mut out = vec![];

        match self.spin {
            2 => {
                if self.pdg_code > 0 {
                    let pol = self.incoming_polarization(mom, helicity);
                    let (u1, u2, u3, u4) = (
                        pol[0].clone(),
                        pol[1].clone(),
                        pol[2].clone(),
                        pol[3].clone(),
                    );
                    out.push((Atom::parse(&format!("u({num},cind(0))")).unwrap(), u1));
                    out.push((Atom::parse(&format!("u({num},cind(1))")).unwrap(), u2));
                    out.push((Atom::parse(&format!("u({num},cind(2))")).unwrap(), u3));
                    out.push((Atom::parse(&format!("u({num},cind(3))")).unwrap(), u4));
                } else {
                    let pol = self.incoming_polarization(mom, helicity);
                    let (v1, v2, v3, v4) = (
                        pol[0].clone(),
                        pol[1].clone(),
                        pol[2].clone(),
                        pol[3].clone(),
                    );
                    out.push((Atom::parse(&format!("vbar({num},cind(0))")).unwrap(), v1));
                    out.push((Atom::parse(&format!("vbar({num},cind(1))")).unwrap(), v2));
                    out.push((Atom::parse(&format!("vbar({num},cind(2))")).unwrap(), v3));
                    out.push((Atom::parse(&format!("vbar({num},cind(3))")).unwrap(), v4));
                }
            }
            3 => {
                let pol = self.incoming_polarization(mom, helicity);
                let (e1, e2, e3, e4) = (
                    pol[0].clone(),
                    pol[1].clone(),
                    pol[2].clone(),
                    pol[3].clone(),
                );
                out.push((Atom::parse(&format!("ϵ({num},cind(0))")).unwrap(), e1));
                out.push((Atom::parse(&format!("ϵ({num},cind(1))")).unwrap(), e2));
                out.push((Atom::parse(&format!("ϵ({num},cind(2))")).unwrap(), e3));
                out.push((Atom::parse(&format!("ϵ({num},cind(3))")).unwrap(), e4));
            }
            _ => {}
        }
        out
    }

    pub fn incoming_polarization<T: FloatLike>(
        &self,
        mom: &FourMomentum<F<T>>,
        helicity: Helicity,
    ) -> Polarization<Complex<F<T>>> {
        Self::incoming_polarization_impl(self.spin, self.pdg_code, mom, helicity)
    }

    pub fn incoming_polarization_impl<T: FloatLike>(
        spin: isize,
        pdg_code: isize,
        mom: &FourMomentum<F<T>>,
        helicity: Helicity,
    ) -> Polarization<Complex<F<T>>> {
        let one: Complex<F<T>> = mom.temporal.value.one().into();
        match spin {
            1 => Polarization::scalar(one),
            2 => {
                if pdg_code > 0 {
                    mom.u(helicity.try_into().unwrap())
                } else {
                    mom.v(helicity.try_into().unwrap()).bar()
                }
            }
            3 => mom.pol(helicity),
            i => panic!("Spin {}/2 not implemented", i - 1),
        }
    }

    pub fn outgoing_polarization_atom(&self, edge_slots: &EdgeSlots<Lorentz>, num: usize) -> Atom {
        match self.spin {
            1 => Atom::parse("1").unwrap(),
            2 => {
                if self.pdg_code > 0 {
                    let mut ubar = FunctionBuilder::new(State::get_symbol("ubar"));
                    ubar = ubar.add_arg(&Atom::new_num(num as i64));
                    ubar = ubar.add_arg(&edge_slots.to_aind_atom());
                    ubar.finish()
                } else {
                    let mut v = FunctionBuilder::new(State::get_symbol("v"));
                    v = v.add_arg(&Atom::new_num(num as i64));
                    v = v.add_arg(&edge_slots.to_aind_atom());
                    v.finish()
                }
            }
            3 => {
                let mut ebar = FunctionBuilder::new(State::get_symbol("ϵbar"));
                ebar = ebar.add_arg(&Atom::new_num(num as i64));
                ebar = ebar.add_arg(&edge_slots.to_aind_atom());
                ebar.finish()
            }
            _ => Atom::parse("1").unwrap(),
        }
    }

    pub fn outgoing_polarization<T: FloatLike>(
        &self,
        mom: &FourMomentum<F<T>>,
        helicity: Helicity,
    ) -> Polarization<Complex<F<T>>> {
        Self::outgoing_polarization_impl(self.spin, self.pdg_code, mom, helicity)
    }

    pub fn outgoing_polarization_impl<T: FloatLike>(
        spin: isize,
        pdg_code: isize,
        mom: &FourMomentum<F<T>>,
        helicity: Helicity,
    ) -> Polarization<Complex<F<T>>> {
        let one: Complex<F<T>> = mom.temporal.value.one().into();

        match spin {
            1 => Polarization::scalar(one),
            2 => {
                if pdg_code > 0 {
                    mom.u(helicity.try_into().unwrap()).bar()
                } else {
                    mom.v(helicity.try_into().unwrap())
                }
            }
            3 => {
                let mut mom = mom.clone();
                mom.temporal = -mom.temporal;
                mom.spatial = -mom.spatial;
                mom.pol(helicity)
            }
            i => panic!("Spin {}/2 not implemented", i - 1),
        }
    }

    pub fn outgoing_polarization_match<T: FloatLike>(
        &self,
        num: usize,
        mom: &FourMomentum<F<T>>,
        helicity: Helicity,
    ) -> Vec<(Atom, Complex<F<T>>)> {
        let mut out = vec![];

        match self.spin {
            2 => {
                if self.pdg_code > 0 {
                    let pol = self.outgoing_polarization(mom, helicity);
                    let (ubar1, ubar2, ubar3, ubar4) = (
                        pol[0].clone(),
                        pol[1].clone(),
                        pol[2].clone(),
                        pol[3].clone(),
                    );
                    out.push((Atom::parse(&format!("ubar({num},cind(0))")).unwrap(), ubar1));
                    out.push((Atom::parse(&format!("ubar({num},cind(1))")).unwrap(), ubar2));
                    out.push((Atom::parse(&format!("ubar({num},cind(2))")).unwrap(), ubar3));
                    out.push((Atom::parse(&format!("ubar({num},cind(3))")).unwrap(), ubar4));
                } else {
                    let pol = self.outgoing_polarization(mom, helicity);
                    let (v1, v2, v3, v4) = (
                        pol[0].clone(),
                        pol[1].clone(),
                        pol[2].clone(),
                        pol[3].clone(),
                    );
                    out.push((Atom::parse(&format!("v({num},cind(0))")).unwrap(), v1));
                    out.push((Atom::parse(&format!("v({num},cind(1))")).unwrap(), v2));
                    out.push((Atom::parse(&format!("v({num},cind(2))")).unwrap(), v3));
                    out.push((Atom::parse(&format!("v({num},cind(3))")).unwrap(), v4));
                }
            }
            3 => {
                let pol = self.outgoing_polarization(mom, helicity);
                let (e1, e2, e3, e4) = (
                    pol[0].clone(),
                    pol[1].clone(),
                    pol[2].clone(),
                    pol[3].clone(),
                );
                out.push((Atom::parse(&format!("ϵbar({num},cind(0))")).unwrap(), e1));
                out.push((Atom::parse(&format!("ϵbar({num},cind(1))")).unwrap(), e2));
                out.push((Atom::parse(&format!("ϵbar({num},cind(2))")).unwrap(), e3));
                out.push((Atom::parse(&format!("ϵbar({num},cind(3))")).unwrap(), e4));
            }
            _ => {}
        }
        out
    }
}

#[test]
fn test_polarization() {
    let mom = FourMomentum::from_args(F(1.0), F(0.0), F(0.0), F(1.0));

    let sqrt_2_inv = F(f64::sqrt(2.)).inv();
    let isqrt_inv = Complex::new_i() * sqrt_2_inv;
    let zero = Complex::new_zero();
    let sqrt_inv = Complex::new_re(sqrt_2_inv);

    let pol_in_plus = Particle::incoming_polarization_impl(3, 1, &mom, Helicity::Plus);
    let pol_out_plus = Particle::outgoing_polarization_impl(3, 1, &mom, Helicity::Plus);
    let pol_in_minus = Particle::incoming_polarization_impl(3, 1, &mom, Helicity::Minus);
    let pol_out_minus = Particle::outgoing_polarization_impl(3, 1, &mom, Helicity::Minus);
    assert_eq!(
        vec![zero, sqrt_inv, isqrt_inv, zero],
        pol_out_minus.tensor.data
    );
    assert_eq!(
        vec![zero, sqrt_inv, -isqrt_inv, zero],
        pol_in_minus.tensor.data
    );
    assert_eq!(
        vec![zero, -sqrt_inv, isqrt_inv, zero],
        pol_out_plus.tensor.data
    );
    assert_eq!(
        vec![zero, -sqrt_inv, -isqrt_inv, zero],
        pol_in_plus.tensor.data
    );
}

#[test]
fn test_polarization_with_transverse_momentum() {
    let mom = FourMomentum::from_args(F(13.0), F(3.0), F(4.0), F(12.0));

    let cmplx_zero = Complex::new_zero();
    let zero = mom.temporal.value.zero();
    let threshold = F(1e-12);

    let pol_in_plus = Particle::incoming_polarization_impl(3, 1, &mom, Helicity::Plus);
    let pol_out_plus = Particle::outgoing_polarization_impl(3, 1, &mom, Helicity::Plus);
    let pol_in_minus = Particle::incoming_polarization_impl(3, 1, &mom, Helicity::Minus);
    let pol_out_minus = Particle::outgoing_polarization_impl(3, 1, &mom, Helicity::Minus);

    let target_pol_in_minus = [
        &cmplx_zero,
        &Complex::new(F(0.3916283711187033), F(0.5656854249492381)),
        &Complex::new(F(0.5221711614916043), F(-0.4242640687119285)),
        &Complex::new(F(-0.2719641466102106), zero),
    ];
    assert!(
        approx_eq_cmplx_vec(
            &target_pol_in_minus,
            &pol_in_minus.tensor.data.iter().collect_vec(),
            &threshold
        ),
        "pol(in,-) does not match:\n{:?}\n!=\n{:?}",
        pol_in_minus.tensor.data,
        target_pol_in_minus
    );

    let target_pol_out_minus = [
        &cmplx_zero,
        &Complex::new(F(0.3916283711187033), F(-0.5656854249492381)),
        &Complex::new(F(0.5221711614916043), F(0.4242640687119285)),
        &Complex::new(F(-0.2719641466102106), zero),
    ];
    assert!(
        approx_eq_cmplx_vec(
            &target_pol_out_minus,
            &pol_out_minus.tensor.data.iter().collect_vec(),
            &threshold
        ),
        "pol(out,-) does not match:\n{:?}\n!=\n{:?}",
        pol_out_minus.tensor.data,
        target_pol_out_minus
    );

    let target_pol_in_plus = [
        &cmplx_zero,
        &Complex::new(F(-0.3916283711187033), F(0.5656854249492381)),
        &Complex::new(F(-0.5221711614916043), F(-0.4242640687119285)),
        &Complex::new(F(0.2719641466102106), zero),
    ];
    assert!(
        approx_eq_cmplx_vec(
            &target_pol_in_plus,
            &pol_in_plus.tensor.data.iter().collect_vec(),
            &threshold
        ),
        "pol(out,-) does not match:\n{:?}\n!=\n{:?}",
        pol_in_plus.tensor.data,
        target_pol_in_plus
    );

    let target_pol_out_plus = [
        &cmplx_zero,
        &Complex::new(F(-0.3916283711187033), F(-0.5656854249492381)),
        &Complex::new(F(-0.5221711614916043), F(0.4242640687119285)),
        &Complex::new(F(0.2719641466102106), zero),
    ];
    assert!(
        approx_eq_cmplx_vec(
            &target_pol_out_plus,
            &pol_out_plus.tensor.data.iter().collect_vec(),
            &threshold
        ),
        "pol(out,-) does not match:\n{:?}\n!=\n{:?}",
        pol_out_plus.tensor.data,
        target_pol_out_plus
    );
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableLorentzStructure {
    name: SmartString<LazyCompact>,
    spins: Vec<isize>,
    structure: SmartString<LazyCompact>,
}

impl SerializableLorentzStructure {
    pub fn from_lorentz_structure(ls: &LorentzStructure) -> SerializableLorentzStructure {
        SerializableLorentzStructure {
            name: ls.name.clone(),
            spins: ls.spins.clone(),
            structure: utils::to_str_expression(&ls.structure).into(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct LorentzStructure {
    pub name: SmartString<LazyCompact>,
    pub spins: Vec<isize>,
    pub structure: Atom,
}

impl LorentzStructure {
    pub fn from_serializable_lorentz_structure(
        ls: &SerializableLorentzStructure,
    ) -> LorentzStructure {
        LorentzStructure {
            name: ls.name.clone(),
            spins: ls.spins.clone(),
            structure: utils::parse_python_expression(ls.structure.as_str()),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableParameter {
    name: SmartString<LazyCompact>,
    lhablock: Option<SmartString<LazyCompact>>,
    lhacode: Option<Vec<usize>>,
    nature: ParameterNature,
    parameter_type: ParameterType,
    value: Option<(F<f64>, F<f64>)>,
    expression: Option<SmartString<LazyCompact>>,
}

impl SerializableParameter {
    pub fn from_parameter(param: &Parameter) -> SerializableParameter {
        SerializableParameter {
            name: param.name.clone(),
            lhablock: param.lhablock.clone(),
            lhacode: param.lhacode.clone(),
            nature: param.nature.clone(),
            parameter_type: param.parameter_type.clone(),
            value: param.value.map(|value| (value.re, value.im)),
            expression: param
                .expression
                .as_ref()
                .map(utils::to_str_expression)
                .map(SmartString::from),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Parameter {
    pub name: SmartString<LazyCompact>,
    pub lhablock: Option<SmartString<LazyCompact>>,
    pub lhacode: Option<Vec<usize>>,
    pub nature: ParameterNature,
    pub parameter_type: ParameterType,
    pub value: Option<Complex<F<f64>>>,
    pub expression: Option<Atom>,
}

impl Parameter {
    pub fn from_serializable_parameter(param: &SerializableParameter) -> Parameter {
        Parameter {
            name: param.name.clone(),
            lhablock: param.lhablock.clone(),
            lhacode: param.lhacode.clone(),
            nature: param.nature.clone(),
            parameter_type: param.parameter_type.clone(),
            value: param.value.map(|value| Complex::new(value.0, value.1)),
            expression: param
                .expression
                .as_ref()
                .map(|expr| utils::parse_python_expression(expr.as_str())),
        }
    }

    pub fn rep_rule(&self) -> Option<[Atom; 2]> {
        let lhs = Atom::parse(&self.name).unwrap();
        let rhs = self.expression.clone();

        //Some([lhs, normalise_complex(&rhs?)])
        Some([lhs, rhs?])
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Order {
    pub name: SmartString<LazyCompact>,
    pub expansion_order: isize,
    pub hierarchy: isize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableModel {
    pub name: SmartString<LazyCompact>,
    pub restriction: Option<SmartString<LazyCompact>>,
    orders: Vec<Order>,
    parameters: Vec<SerializableParameter>,
    particles: Vec<SerializableParticle>,
    propagators: Vec<SerializablePropagator>,
    lorentz_structures: Vec<SerializableLorentzStructure>,
    couplings: Vec<SerializableCoupling>,
    vertex_rules: Vec<SerializableVertexRule>,
}

impl SerializableModel {
    pub fn from_file(file_path: String) -> Result<SerializableModel, Report> {
        let f = File::open(file_path.clone())
            .wrap_err_with(|| format!("Could not open model yaml file {}", file_path))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .map_err(|e| eyre!(format!("Error parsing model yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<SerializableModel, Report> {
        serde_yaml::from_str(yaml_str.as_str())
            .map_err(|e| eyre!(format!("Error parsing model yaml: {}", e)))
            .suggestion("Is it a correct yaml file")
    }

    pub fn from_model(model: &Model) -> SerializableModel {
        SerializableModel {
            name: model.name.clone(),
            restriction: model.restriction.clone(),
            orders: model
                .orders
                .iter()
                .map(|order| order.as_ref().clone())
                .collect(),
            parameters: model
                .parameters
                .iter()
                .map(|parameter| SerializableParameter::from_parameter(parameter.as_ref()))
                .collect(),
            particles: model
                .particles
                .iter()
                .map(|particle| SerializableParticle::from_particle(particle.as_ref()))
                .collect(),
            propagators: model
                .propagators
                .iter()
                .map(|propagator| SerializablePropagator::from_propagator(propagator.as_ref()))
                .collect(),
            lorentz_structures: model
                .lorentz_structures
                .iter()
                .map(|lorentz_structure| {
                    SerializableLorentzStructure::from_lorentz_structure(lorentz_structure.as_ref())
                })
                .collect(),
            couplings: model
                .couplings
                .iter()
                .map(|coupling| SerializableCoupling::from_coupling(coupling.as_ref()))
                .collect(),
            vertex_rules: model
                .vertex_rules
                .iter()
                .map(|vertex_rule| SerializableVertexRule::from_vertex_rule(vertex_rule.as_ref()))
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Model {
    pub name: SmartString<LazyCompact>,
    pub restriction: Option<SmartString<LazyCompact>>,
    pub orders: Vec<Arc<Order>>,
    pub parameters: Vec<Arc<Parameter>>,
    pub particles: Vec<Arc<Particle>>,
    pub propagators: Vec<Arc<Propagator>>,
    pub lorentz_structures: Vec<Arc<LorentzStructure>>,
    pub couplings: Vec<Arc<Coupling>>,
    pub vertex_rules: Vec<Arc<VertexRule>>,
    pub order_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub parameter_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub lorentz_structure_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub particle_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub particle_pdg_to_position: HashMap<isize, usize, RandomState>,
    pub propagator_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub coupling_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
    pub vertex_rule_name_to_position: HashMap<SmartString<LazyCompact>, usize, RandomState>,
}

impl Default for Model {
    fn default() -> Self {
        Model {
            name: SmartString::<LazyCompact>::from("ModelNotLoaded"),
            restriction: None,
            orders: vec![],
            parameters: vec![],
            particles: vec![],
            propagators: vec![],
            lorentz_structures: vec![],
            couplings: vec![],
            vertex_rules: vec![],
            order_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            parameter_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            lorentz_structure_name_to_position: HashMap::<
                SmartString<LazyCompact>,
                usize,
                RandomState,
            >::default(),
            particle_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            particle_pdg_to_position: HashMap::<isize, usize, RandomState>::default(),
            propagator_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            coupling_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
            vertex_rule_name_to_position:
                HashMap::<SmartString<LazyCompact>, usize, RandomState>::default(),
        }
    }
}
impl Model {
    pub fn substitute_model_params(&self, atom: &Atom) -> Atom {
        let mut sub_atom = atom.clone();
        for cpl in self.couplings.iter() {
            let [pattern, rhs] = cpl.rep_rule();

            sub_atom =
                sub_atom.replace_all(&pattern.into_pattern(), &rhs.into_pattern(), None, None);
        }

        for para in self.parameters.iter() {
            if let Some([pattern, rhs]) = para.rep_rule() {
                sub_atom =
                    sub_atom.replace_all(&pattern.into_pattern(), &rhs.into_pattern(), None, None);
            }
        }
        sub_atom
    }

    pub fn dependent_coupling_replacements(&self) -> Vec<(Pattern, Pattern)> {
        let mut reps = vec![];
        for cpl in self.couplings.iter().filter(|c| c.value.is_none()) {
            let [pattern, rhs] = cpl.rep_rule();
            reps.push((pattern.into_pattern(), rhs.into_pattern()));
        }
        reps
    }

    pub fn internal_parameter_replacements(&self) -> Vec<(Pattern, Pattern)> {
        let mut reps = vec![];
        for para in self
            .parameters
            .iter()
            .filter(|p| matches!(p.nature, ParameterNature::Internal))
        {
            if let Some([pattern, rhs]) = para.rep_rule() {
                reps.push((pattern.into_pattern(), rhs.into_pattern()));
            }
        }
        reps
    }

    pub fn valued_coupling_re_im_split(
        &self,
        fn_map: &mut FunctionMap<Rational>,
    ) -> Vec<(Pattern, Pattern)> {
        let mut reps = vec![];
        for cpl in self.couplings.iter().filter(|c| c.value.is_some()) {
            let lhs = Atom::parse(&cpl.name).unwrap().into_pattern();
            if let Some(value) = cpl.value {
                let rhs = if value.im == 0.0 {
                    let name = Atom::new_var(State::get_symbol(format!("{}_re", cpl.name)));
                    fn_map.add_constant(name.clone(), Rational::from(value.re));
                    name.into_pattern()
                } else if value.re == 0.0 {
                    let name = Atom::new_var(State::get_symbol(format!("{}_im", cpl.name)));
                    fn_map.add_constant(name.clone(), Rational::from(value.im));
                    name.into_pattern()
                } else {
                    let name_re = Atom::new_var(State::get_symbol(cpl.name.clone() + "_re"));
                    fn_map.add_constant(name_re.clone(), Rational::from(value.re));
                    let name_im = Atom::new_var(State::get_symbol(cpl.name.clone() + "_im"));
                    fn_map.add_constant(name_im.clone(), Rational::from(value.im));
                    let i = Atom::new_var(State::I);
                    (&name_re + i * &name_im).into_pattern()
                };
                reps.push((lhs, rhs));
            }
        }
        reps
    }

    pub fn substitute_split_model_params(&self, atom: &Atom) -> Atom {
        atom.clone()
    }

    pub fn valued_parameter_re_im_split(
        &self,
        fn_map: &mut FunctionMap<Rational>,
    ) -> Vec<(Pattern, Pattern)> {
        let mut reps = vec![];
        for param in self.parameters.iter().filter(|p| p.value.is_some()) {
            let lhs = Atom::parse(&param.name).unwrap().into_pattern();
            if let Some(value) = param.value {
                let rhs = match param.parameter_type {
                    ParameterType::Imaginary => {
                        if value.re.is_zero() {
                            let name =
                                Atom::new_var(State::get_symbol(format!("{}_im", param.name)));
                            fn_map.add_constant(name.clone(), Rational::from(value.im.0));
                            name.into_pattern()
                        } else {
                            let name_re =
                                Atom::new_var(State::get_symbol(param.name.clone() + "_re"));
                            fn_map.add_constant(name_re.clone(), Rational::from(value.re.0));
                            let name_im =
                                Atom::new_var(State::get_symbol(param.name.clone() + "_im"));
                            fn_map.add_constant(name_im.clone(), Rational::from(value.im.0));
                            let i = Atom::new_var(State::I);
                            (&name_re + i * &name_im).into_pattern()
                        }
                    }
                    ParameterType::Real => {
                        let name = Atom::new_var(State::get_symbol(format!("{}_re", param.name)));
                        fn_map.add_constant(name.clone(), Rational::from(value.re.0));
                        name.into_pattern()
                    }
                };
                reps.push((lhs, rhs));
            }
        }
        reps
    }

    pub fn evaluate_couplings(&self, atom: Atom) -> Atom {
        // let mut atom = atom;
        // for cpl in self.couplings.iter() {
        //     if let Some(value) = cpl.value {
        //         let pat = Atom::parse(&cpl.name).unwrap().into_pattern();

        //         let re = Atom::new_num(value);
        //         atom = atom.replace_all(&pattern.into_pattern(), &rhs.into_pattern(), None, None);
        //     }
        //     let [pattern, rhs] = cpl.rep_rule();
        //     atom = atom.replace_all(&pattern.into_pattern(), &rhs.into_pattern(), None, None);
        // }
        atom
    }

    pub fn append_coupling_eval<'a, T: FloatLike>(
        &'a self,
        const_map: &mut HashMap<AtomView<'a>, Complex<F<T>>>,
    ) {
        // let mut atom = atom;
        for cpl in self.couplings.iter() {
            if let Some(value) = cpl.value {
                let val = Complex::new(F::<T>::from_f64(value.re), F::<T>::from_f64(value.im));
                const_map.insert(cpl.expression.as_view(), val);
            }
        }
    }

    pub fn append_parameter_map(&self, const_map: &mut AHashMap<Atom, Complex<F<f64>>>) {
        // let mut atom = atom;
        for cpl in self.parameters.iter() {
            if let Some(value) = cpl.value {
                let key = Atom::parse(&cpl.name).unwrap();
                const_map.insert(key, value);
            }
        }
    }
    pub fn is_empty(&self) -> bool {
        self.name == "ModelNotLoaded" || self.particles.is_empty()
    }

    pub fn export_coupling_replacement_rules(
        &self,
        export_root: &str,
        print_ops: PrintOptions,
    ) -> Result<(), Report> {
        let path = Path::new(export_root).join("sources").join("model");

        if !path.exists() {
            fs::create_dir_all(&path)?;
        }
        let mut reps = Vec::new();

        for cpl in self.couplings.iter() {
            reps.push(
                cpl.rep_rule()
                    .map(|a| format!("{}", AtomPrinter::new_with_options(a.as_view(), print_ops))),
            );
        }

        for para in self.parameters.iter() {
            if let Some(rule) = para.rep_rule() {
                reps.push(
                    rule.map(|a| {
                        format!("{}", AtomPrinter::new_with_options(a.as_view(), print_ops))
                    }),
                );
            }
        }

        fs::write(
            path.join("model_replacements.json"),
            serde_json::to_string_pretty(&reps)?,
        )?;

        Ok(())
    }

    pub fn from_serializable_model(serializable_model: SerializableModel) -> Model {
        let mut model: Model = Model::default();
        model.name = serializable_model.name;
        model.restriction = serializable_model.restriction;

        // Extract coupling orders
        model.orders = serializable_model
            .orders
            .iter()
            .enumerate()
            .map(|(i_order, serializable_order)| {
                let order = Arc::new(Order {
                    name: serializable_order.name.clone(),
                    expansion_order: serializable_order.expansion_order,
                    hierarchy: serializable_order.hierarchy,
                });
                model
                    .order_name_to_position
                    .insert(order.name.clone(), i_order);
                order
            })
            .collect();

        // Extract parameters
        model.parameters = serializable_model
            .parameters
            .iter()
            .enumerate()
            .map(|(i_param, serializable_param)| {
                let parameter =
                    Arc::new(Parameter::from_serializable_parameter(serializable_param));
                model
                    .parameter_name_to_position
                    .insert(parameter.name.clone(), i_param);
                parameter
            })
            .collect();

        // Extract particles
        model.particles = serializable_model
            .particles
            .iter()
            .enumerate()
            .map(|(i_part, serializable_particle)| {
                let particle = Arc::new(Particle::from_serializable_particle(
                    &model,
                    serializable_particle,
                ));
                model
                    .particle_name_to_position
                    .insert(particle.name.clone(), i_part);
                model
                    .particle_pdg_to_position
                    .insert(particle.pdg_code, i_part);
                particle
            })
            .collect();

        // Extract propagators

        model.propagators = serializable_model
            .propagators
            .iter()
            .enumerate()
            .map(|(i_prop, serializable_propagator)| {
                let propagator = Arc::new(Propagator::from_serializable_propagator(
                    &model,
                    serializable_propagator,
                ));
                model
                    .propagator_name_to_position
                    .insert(propagator.name.clone(), i_prop);
                propagator
            })
            .collect();

        // Extract Lorentz structures
        model.lorentz_structures = serializable_model
            .lorentz_structures
            .iter()
            .enumerate()
            .map(|(i_lor, serializable_lorentz_structure)| {
                let lorentz_structure =
                    Arc::new(LorentzStructure::from_serializable_lorentz_structure(
                        serializable_lorentz_structure,
                    ));
                model
                    .lorentz_structure_name_to_position
                    .insert(lorentz_structure.name.clone(), i_lor);
                lorentz_structure
            })
            .collect();

        // Extract couplings
        model.couplings = serializable_model
            .couplings
            .iter()
            .enumerate()
            .map(|(i_coupl, serializable_coupling)| {
                let coupling =
                    Arc::new(Coupling::from_serializable_coupling(serializable_coupling));
                model
                    .coupling_name_to_position
                    .insert(coupling.name.clone(), i_coupl);
                coupling
            })
            .collect();

        // Extract vertex rules
        model.vertex_rules = serializable_model
            .vertex_rules
            .iter()
            .enumerate()
            .map(|(i_vr, serializable_vertex_rule)| {
                let vertex_rule = Arc::new(VertexRule::from_serializable_vertex_rule(
                    &model,
                    serializable_vertex_rule,
                ));
                model
                    .vertex_rule_name_to_position
                    .insert(vertex_rule.name.clone(), i_vr);
                vertex_rule
            })
            .collect();

        model
    }

    pub fn to_serializable(&self) -> SerializableModel {
        SerializableModel::from_model(self)
    }

    pub fn to_yaml(&self) -> Result<String, Error> {
        serde_yaml::to_string(&self.to_serializable())
    }

    pub fn from_file(file_path: String) -> Result<Model, Report> {
        SerializableModel::from_file(file_path).map(Model::from_serializable_model)
    }

    pub fn from_yaml_str(yaml_str: String) -> Result<Model, Report> {
        SerializableModel::from_yaml_str(yaml_str).map(Model::from_serializable_model)
    }

    #[inline]
    pub fn get_particle(&self, name: &SmartString<LazyCompact>) -> Arc<Particle> {
        if let Some(position) = self.particle_name_to_position.get(name) {
            self.particles[*position].clone()
        } else {
            panic!("Particle '{}' not found in model '{}'.", name, self.name);
        }
    }
    #[inline]
    pub fn get_particle_from_pdg(&self, pdg: isize) -> Arc<Particle> {
        if let Some(position) = self.particle_pdg_to_position.get(&pdg) {
            self.particles[*position].clone()
        } else {
            panic!(
                "Particle with PDG {} not found in model '{}'.",
                pdg, self.name
            );
        }
    }

    #[inline]
    pub fn get_propagator(&self, name: &SmartString<LazyCompact>) -> Arc<Propagator> {
        if let Some(position) = self.propagator_name_to_position.get(name) {
            self.propagators[*position].clone()
        } else {
            panic!("Propagator '{}' not found in model '{}'.", name, self.name);
        }
    }

    #[inline]
    pub fn get_parameter(&self, name: &SmartString<LazyCompact>) -> Arc<Parameter> {
        if let Some(position) = self.parameter_name_to_position.get(name) {
            self.parameters[*position].clone()
        } else {
            panic!("Parameter '{}' not found in model '{}'.", name, self.name);
        }
    }
    #[inline]
    pub fn get_order(&self, name: &SmartString<LazyCompact>) -> Arc<Order> {
        if let Some(position) = self.order_name_to_position.get(name) {
            self.orders[*position].clone()
        } else {
            panic!(
                "Coupling order '{}' not found in model '{}'.",
                name, self.name
            );
        }
    }
    #[inline]
    pub fn get_lorentz_structure(&self, name: &SmartString<LazyCompact>) -> Arc<LorentzStructure> {
        if let Some(position) = self.lorentz_structure_name_to_position.get(name) {
            self.lorentz_structures[*position].clone()
        } else {
            panic!(
                "Lorentz structure '{}' not found in model '{}'.",
                name, self.name
            );
        }
    }
    #[inline]
    pub fn get_coupling(&self, name: &SmartString<LazyCompact>) -> Arc<Coupling> {
        if let Some(position) = self.coupling_name_to_position.get(name) {
            self.couplings[*position].clone()
        } else {
            panic!("Coupling '{}' not found in model '{}'.", name, self.name);
        }
    }
    #[inline]
    pub fn get_vertex_rule(&self, name: &SmartString<LazyCompact>) -> Arc<VertexRule> {
        if let Some(position) = self.vertex_rule_name_to_position.get(name) {
            self.vertex_rules[*position].clone()
        } else {
            panic!("Vertex rule '{}' not found in model '{}'.", name, self.name);
        }
    }
}
