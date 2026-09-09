use symbolica::atom::Atom;

/// A loop propagator `1 / ((l + momentum)^2 - mass_sq)`.
#[derive(Debug, Clone)]
pub struct Propagator {
    pub momentum: Atom,
    /// The propagator mass squared.
    pub mass_sq: Atom,
}

/// An irreducible scalar product.
#[derive(Debug, Clone)]
pub struct Isp {
    pub expression: Atom,
}

/// External kinematics.
#[derive(Debug, Clone, Default)]
pub struct Kinematics {
    /// The `C(N, 2)` pairwise invariants `(r_i - r_j)^2`
    pub invariants: Vec<Atom>,
}

/// A target integral
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Integral {
    /// Power of each propagator
    pub propagator_exponents: Vec<i32>,
    /// Reserved.
    pub isp_exponents: Vec<i32>,
}

#[derive(Debug, Clone)]
pub struct IntegralFamily {
    /// Loop propagators
    pub propagators: Vec<Propagator>,
    /// Reserved.
    pub isps: Vec<Isp>,
    /// External kinematics.
    pub kinematics: Kinematics,
    /// Target integrals
    pub targets: Vec<Integral>,
    /// Numerator: a polynomial in `dot(k,k)` and `dot(k,q_i)`
    pub numerator: Atom,
}
