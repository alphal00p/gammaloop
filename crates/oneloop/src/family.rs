use symbolica::atom::Atom;

#[derive(Debug, Clone)]
pub struct Propagator {
    pub momentum: Atom,
    pub mass_sq: Atom,
}

#[derive(Debug, Clone)]
pub struct Isp {
    pub expression: Atom,
}

#[derive(Debug, Clone, Default)]
pub struct Kinematics {
    pub invariants: Vec<Atom>,
    pub masses_sq: Vec<Atom>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Integral {
    pub propagator_exponents: Vec<i32>,
    pub isp_exponents: Vec<i32>,
}

#[derive(Debug, Clone)]
pub struct IntegralFamily {
    pub propagators: Vec<Propagator>,
    pub isps: Vec<Isp>,
    pub kinematics: Kinematics,
    pub targets: Vec<Integral>,
    pub numerator: Atom,
}
