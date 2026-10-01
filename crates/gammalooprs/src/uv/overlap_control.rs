//! Exact points for a negative control after physical soft-coefficient extraction.
//! A nonzero specialization disproves an identity; a zero specialization proves
//! nothing. Momentum substitution preserves numerator products and tensor slots.

use color_eyre::eyre::{Result, ensure, eyre};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    domains::{float::Complex, rational::Rational},
    parse_lit, symbol,
};

use crate::{cff::expression::OrientationID, model::UFOSymbol, utils::GS};

const MAX_BYTES: usize = 8 * 1024 * 1024;
const MAX_DEPTH: usize = 128;

pub(crate) struct ExactControlPoint {
    bubble_x: Rational,
}

impl ExactControlPoint {
    pub(crate) fn new() -> Self {
        Self { bubble_x: 4.into() }
    }

    /// Both bubble energies and their mass-three counterparts are rational.
    /// The second point separates the bubble from the hard energy |r+q|=4.
    pub(crate) fn points() -> [Self; 2] {
        [
            Self::new(),
            Self {
                bubble_x: (9, 4).into(),
            },
        ]
    }

    /// The input has already undergone Q6=λq, Q5=r+(1−λ)q, Q8=b and
    /// coefficient extraction. Do not shift Q5 a second time here.
    pub(crate) fn apply(&self, atom: &Atom) -> Result<Atom> {
        ensure!(
            atom.as_view().get_byte_size() <= MAX_BYTES,
            "control point input budget"
        );
        self.substitute(atom.as_view(), 0)
    }

    /// A nonzero exact Q(i) coefficient times an integer power of pi is a
    /// point witness. Unknown parameters never silently become numbers.
    pub(crate) fn nonzero_scalar_witness(&self, atom: &Atom) -> Result<bool> {
        let value = self.apply(atom)?;
        let (coefficient, _) = Self::rational_pi_monomial(&value)?;
        Ok(!coefficient.is_zero())
    }

    /// Keep the known nonzero pi factor symbolic. This is not evaluation at
    /// pi=1, nor an appeal to transcendence for sums of different pi powers.
    pub(crate) fn rational_pi_monomial(atom: &Atom) -> Result<(Atom, i64)> {
        ensure!(
            atom.as_view().get_byte_size() <= MAX_BYTES,
            "control pi input budget"
        );
        Self::pi_monomial(atom.as_view(), 0)
    }

    fn pi_monomial(atom: AtomView<'_>, depth: usize) -> Result<(Atom, i64)> {
        ensure!(depth <= MAX_DEPTH, "control pi recursion budget");
        let (coefficient, exponent) = match atom {
            AtomView::Num(_) => {
                Complex::<Rational>::try_from(atom)
                    .map_err(|_| eyre!("unproven control point: coefficient is not exact Q(i)"))?;
                (atom.to_owned(), 0)
            }
            AtomView::Var(variable) if variable.get_symbol() == Symbol::PI => (Atom::one(), 1),
            AtomView::Add(sum) => {
                let mut result = (Atom::Zero, 0);
                for term in sum {
                    let (coefficient, exponent) = Self::pi_monomial(term, depth + 1)?;
                    if coefficient.is_zero() {
                        continue;
                    }
                    if result.0.is_zero() {
                        result = (coefficient, exponent);
                    } else {
                        ensure!(
                            result.1 == exponent,
                            "unproven control point: sum has different pi powers"
                        );
                        result.0 += coefficient;
                    }
                }
                result
            }
            AtomView::Mul(product) => {
                let mut result = (Atom::one(), 0i64);
                for factor in product {
                    let (coefficient, exponent) = Self::pi_monomial(factor, depth + 1)?;
                    result.0 *= coefficient;
                    result.1 = result
                        .1
                        .checked_add(exponent)
                        .ok_or_else(|| eyre!("control pi exponent overflow"))?;
                }
                result
            }
            AtomView::Pow(power) => {
                let exponent = i64::try_from(power.get_exp())
                    .map_err(|_| eyre!("unproven control point: noninteger pi power"))?;
                ensure!(exponent.unsigned_abs() <= 64, "control pi exponent budget");
                let (coefficient, base_exponent) = Self::pi_monomial(power.get_base(), depth + 1)?;
                ensure!(
                    exponent >= 0 || !coefficient.is_zero(),
                    "singular control point: inverse base vanishes"
                );
                ensure!(
                    coefficient
                        .as_view()
                        .get_byte_size()
                        .saturating_mul(exponent.unsigned_abs() as usize)
                        <= MAX_BYTES,
                    "control pi power budget"
                );
                (
                    coefficient.pow(exponent),
                    base_exponent
                        .checked_mul(exponent)
                        .ok_or_else(|| eyre!("control pi exponent overflow"))?,
                )
            }
            _ => {
                return Err(eyre!(
                    "unproven control point: scalar is not a rational pi monomial: {atom}"
                ));
            }
        };
        ensure!(
            coefficient.as_view().get_byte_size() <= MAX_BYTES,
            "control pi coefficient budget"
        );
        Ok(if coefficient.is_zero() {
            (Atom::Zero, 0)
        } else {
            (coefficient, exponent)
        })
    }

    fn substitute(&self, atom: AtomView<'_>, depth: usize) -> Result<Atom> {
        ensure!(depth <= MAX_DEPTH, "control point recursion budget");
        let result = match atom {
            AtomView::Num(_) => {
                Complex::<Rational>::try_from(atom).map_err(|_| {
                    eyre!("unproven control point: inexact or undefined coefficient")
                })?;
                atom.to_owned()
            }
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                if [GS.m_uv_vacuum, GS.m_uv_expansion].contains(&symbol) {
                    Atom::num(3)
                } else if symbol == UFOSymbol::from("GC_1").0 {
                    // Exact SM definition in assets/models/json/sm/sm.json.
                    let ee =
                        self.substitute(Atom::from(UFOSymbol::from("ee")).as_view(), depth + 1)?;
                    parse_lit!(-1i / 3) * ee
                } else if symbol == UFOSymbol::from("GC_11").0 {
                    let strong =
                        self.substitute(Atom::from(UFOSymbol::from("G")).as_view(), depth + 1)?;
                    parse_lit!(1i) * strong
                } else if [UFOSymbol::from("ee").0, UFOSymbol::from("G").0].contains(&symbol) {
                    // The SM declares ee=2*sqrt(aEW*pi), G=2*sqrt(aS*pi).
                    // These two independent bare couplings are explicitly set
                    // to one; no other model parameter is assigned a value.
                    Atom::one()
                } else {
                    atom.to_owned()
                }
            }
            AtomView::Fun(call) => {
                ensure!(
                    ![
                        Symbol::IF,
                        OrientationID::symbol(),
                        GS.theta,
                        GS.orientation_delta,
                        symbol!("gammalooprs::uv::numerator_family"),
                    ]
                    .contains(&call.get_symbol()),
                    "unproven control point: guard or unresolved numerator family"
                );
                if [GS.emr_vec, GS.emr_mom].contains(&call.get_symbol()) {
                    ensure!(
                        call.get_nargs() == 2,
                        "unproven control point: momentum arity"
                    );
                    let edge = i64::try_from(call.get(0))
                        .map_err(|_| eyre!("unproven control point: nonliteral momentum owner"))?;
                    self.momentum(edge, call.get(1), call.get_symbol() == GS.emr_vec)?
                } else if call.get_symbol() == GS.energy_surface {
                    ensure!(
                        call.get_nargs() == 2 && call.get(0).is_zero(),
                        "unproven control point: energy has not been physically normalized"
                    );
                    let invariant = self.substitute(call.get(1), depth + 1)?;
                    Self::positive_root(&invariant)?
                } else {
                    let mut builder = FunctionBuilder::new(call.get_symbol());
                    for argument in call {
                        builder = builder.add_arg(self.substitute(argument, depth + 1)?);
                    }
                    builder.finish()
                }
            }
            AtomView::Add(sum) => Atom::add_many(
                sum.iter()
                    .map(|part| self.substitute(part, depth + 1))
                    .collect::<Result<Vec<_>>>()?,
            ),
            AtomView::Mul(product) => Atom::mul_many(
                product
                    .iter()
                    .map(|part| self.substitute(part, depth + 1))
                    .collect::<Result<Vec<_>>>()?,
            ),
            AtomView::Pow(power) => {
                let exponent = i64::try_from(power.get_exp())
                    .map_err(|_| eyre!("unproven control point: noninteger power"))?;
                ensure!(
                    exponent.unsigned_abs() <= 64,
                    "control point exponent budget"
                );
                let base = self.substitute(power.get_base(), depth + 1)?;
                if exponent < 0 {
                    // Check before constructing the reciprocal, including when
                    // a different factor in the same product becomes zero.
                    let (coefficient, _) = Self::rational_pi_monomial(&base)?;
                    ensure!(
                        !coefficient.is_zero(),
                        "singular control point: inverse base vanishes"
                    );
                }
                ensure!(
                    base.as_view()
                        .get_byte_size()
                        .saturating_mul(exponent.unsigned_abs() as usize)
                        <= MAX_BYTES,
                    "control point power budget"
                );
                base.pow(exponent)
            }
        };
        ensure!(
            result.as_view().get_byte_size() <= MAX_BYTES,
            "control point output budget"
        );
        Ok(result)
    }

    fn momentum(&self, edge: i64, index: AtomView<'_>, spatial: bool) -> Result<Atom> {
        let components: [Atom; 4] = match edge {
            0 => [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero],
            6 => [Atom::Zero, Atom::Zero, Atom::one(), Atom::Zero],
            5 => [Atom::Zero, Atom::num(4), Atom::num(-1), Atom::Zero],
            8 => [
                Atom::Zero,
                Atom::num(self.bubble_x.clone()),
                Atom::Zero,
                Atom::Zero,
            ],
            _ => {
                return Err(eyre!(
                    "unproven control point: unrouted momentum owner {edge}"
                ));
            }
        };
        if let Some(component) = (0..=3).find(|component| index == GS.cind(*component).as_view()) {
            if spatial && component == 0 {
                return Ok(Atom::Zero);
            }
            ensure!(
                spatial || edge == 0 || component != 0,
                "unproven control point: internal temporal momentum remains unresolved"
            );
            return Ok(components[component].clone());
        }
        ensure!(
            spatial || edge == 0,
            "unproven control point: internal four-momentum remains unresolved"
        );
        // Preserve the actual tensor index/representation, not its printed name.
        let start = usize::from(spatial);
        Ok(Atom::add_many((start..=3).map(|component| {
            &components[component]
                * GS.delta_vec
                    .call_args([GS.cind(component), index.to_owned()])
        })))
    }

    fn positive_root(invariant: &Atom) -> Result<Atom> {
        let value = Rational::try_from(invariant).map_err(|_| {
            eyre!("unproven control point: nonrational energy invariant: {invariant}")
        })?;
        ensure!(
            value > 0,
            "unproven control point: nonpositive energy invariant"
        );
        let numerator = value.numerator().root(2);
        let denominator = value.denominator().root(2);
        ensure!(
            &numerator * &numerator == value.numerator()
                && &denominator * &denominator == value.denominator(),
            "unproven control point: energy invariant is not a rational square: {invariant}"
        );
        Ok(Atom::num(Rational::new(numerator, denominator)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;
    use symbolica::function;

    #[test]
    fn physical_components_and_positive_energy_branch_are_exact() -> Result<()> {
        test_initialise()?;
        let points = ExactControlPoint::points();
        let index = parse_lit!(spenso::mink(4, 17));
        let q = function!(GS.emr_vec, 6, &index);
        assert_eq!(
            points[0].apply(&q)?,
            GS.delta_vec.call_args([GS.cind(2), index.clone()])
        );
        let q0 = function!(GS.emr_mom, 0, &index);
        assert_eq!(points[0].apply(&q0)?, GS.energy_delta(&index));
        let hard_x = function!(GS.emr_mom, 5, GS.cind(1)) + function!(GS.emr_mom, 6, GS.cind(1));
        let hard_y = function!(GS.emr_mom, 5, GS.cind(2)) + function!(GS.emr_mom, 6, GS.cind(2));
        let hard = function!(
            GS.energy_surface,
            0,
            hard_x.pow(2) + hard_y.pow(2) + Atom::var(GS.m_uv_vacuum).pow(2)
        );
        let bubble = function!(
            GS.energy_surface,
            0,
            function!(GS.emr_mom, 8, GS.cind(1)).pow(2) + Atom::var(GS.m_uv_vacuum).pow(2)
        );
        assert_eq!(points[0].apply(&hard)?, Atom::num(5));
        assert_eq!(points[0].apply(&bubble)?, Atom::num(5));
        assert_eq!(points[1].apply(&bubble)?, Atom::num((15, 4)));
        assert!(points[0].apply(&(&hard - &bubble).pow(-1)).is_err());
        assert_eq!(
            points[1].apply(&(&hard - &bubble).pow(-1))?,
            Atom::num((4, 5))
        );
        Ok(())
    }

    #[test]
    fn control_rejects_singular_or_unsupported_points_without_assigning_unknowns() -> Result<()> {
        test_initialise()?;
        let point = ExactControlPoint::new();
        let unknown = parse_lit!(unassigned_parameter);
        assert_eq!(point.apply(&unknown)?, unknown);
        assert!(point.nonzero_scalar_witness(&unknown).is_err());
        let qx = function!(GS.emr_mom, 6, GS.cind(1));
        // Even the vanishing numerator does not license a singular denominator.
        let denominator = function!(GS.emr_mom, 8, GS.cind(1)) - 4;
        assert!(point.apply(&(qx * denominator.pow(-1))).is_err());
        assert!(point.apply(&function!(GS.energy_surface, 0, 2)).is_err());
        assert!(point.apply(&function!(GS.energy_surface, 0, 0)).is_err());
        assert!(point.apply(&function!(GS.energy_surface, 1, 4)).is_err());
        assert!(point.apply(&function!(GS.emr_mom, 5, GS.cind(0))).is_err());
        assert!(point.apply(&Atom::num(0.25f64)).is_err());
        assert!(point.nonzero_scalar_witness(&parse_lit!(1i / 3))?);
        assert!(!point.nonzero_scalar_witness(&Atom::Zero)?);
        assert_eq!(
            point.apply(&(Atom::from(UFOSymbol::from("ee")) * Atom::from(UFOSymbol::from("G"))))?,
            Atom::one()
        );
        Ok(())
    }

    #[test]
    fn sm_coupling_definitions_and_common_pi_factors_remain_exact() -> Result<()> {
        test_initialise()?;
        let point = ExactControlPoint::new();
        let gc1 = Atom::from(UFOSymbol::from("GC_1"));
        let gc11 = Atom::from(UFOSymbol::from("GC_11"));
        assert_eq!(point.apply(&gc1)?, parse_lit!(-1i / 3));
        assert_eq!(point.apply(&gc11)?, parse_lit!(1i));
        let pi = Atom::var(Symbol::PI);
        let numerator = Atom::num(2) * gc1.pow(2) * gc11.pow(4);
        let value = point.apply(&(numerator * pi.pow(-9)))?;
        assert_eq!(
            ExactControlPoint::rational_pi_monomial(&value)?,
            (parse_lit!(-2 / 9), -9)
        );
        assert!(value.contains_symbol(Symbol::PI));
        let same_power = (&pi * 2 + &pi * parse_lit!(1i)) / pi.pow(4);
        assert_eq!(
            ExactControlPoint::rational_pi_monomial(&same_power)?,
            (parse_lit!(2 + 1i), -3)
        );
        assert!(point.nonzero_scalar_witness(&same_power)?);
        Ok(())
    }

    #[test]
    fn pi_witness_checks_inverse_domain_and_rejects_unknown_couplings() -> Result<()> {
        test_initialise()?;
        let point = ExactControlPoint::new();
        let pi = Atom::var(Symbol::PI);
        let qy = function!(GS.emr_mom, 6, GS.cind(2));
        let zero_base = pi.pow(2) * (qy - 1);
        assert!(point.apply(&zero_base.pow(-1)).is_err());
        let unknown = Atom::from(UFOSymbol::from("GC_unknown"));
        assert_eq!(point.apply(&unknown)?, unknown);
        assert!(point.nonzero_scalar_witness(&unknown).is_err());
        assert!(point.apply(&(unknown * &pi).pow(-1)).is_err());
        assert!(ExactControlPoint::rational_pi_monomial(&(&pi + 1)).is_err());
        assert!(ExactControlPoint::rational_pi_monomial(&pi.pow(Atom::num((1, 2)))).is_err());
        assert_eq!(
            ExactControlPoint::rational_pi_monomial(&Atom::Zero)?,
            (Atom::Zero, 0)
        );
        Ok(())
    }
}
