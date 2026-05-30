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
}

#[cfg(test)]
mod tests {
    use super::{Integral, IntegralFamily, Kinematics, Propagator};
    use symbolica::atom::Atom;

    #[test]
    fn bubble_family_has_two_propagators_and_one_target() {
        let family = IntegralFamily {
            propagators: vec![
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: Atom::Zero,
                },
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: Atom::Zero,
                },
            ],
            isps: vec![],
            kinematics: Kinematics::default(),
            targets: vec![Integral {
                propagator_exponents: vec![1, 1],
                isp_exponents: vec![],
            }],
        };
        assert_eq!(family.propagators.len(), 2);
        assert_eq!(family.targets.len(), 1);
        assert_eq!(family.targets[0].propagator_exponents, vec![1, 1]);
    }
}
