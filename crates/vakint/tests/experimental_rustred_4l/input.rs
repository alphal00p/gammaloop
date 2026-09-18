//! Test input construction shared by the candidate and independent FMFT oracle.

use rustred::family::IntegralFamily;
use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
use std::sync::Arc;
use symbolica::atom::{Atom, AtomCore};
use symbolica::function;
use vakint::symbols::S;
use vakint::vakint_parse as vk_parse;

pub struct ParentInput {
    pub family: Arc<IntegralFamily>,
    pub physical_momenta: Vec<Atom>,
    momenta: Vec<Atom>,
    edges: Vec<(i64, i64)>,
}

impl ParentInput {
    pub fn from_csv(source: &str) -> Self {
        let rows = source
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|line| {
                line.split(',')
                    .map(|entry| entry.trim().parse::<i64>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            rows.len(),
            10,
            "complete four-loop family needs ten denominators"
        );
        let mut physical_momenta = Vec::new();
        let mut momenta = Vec::new();
        let mut edges = Vec::new();
        let mut propagators = Vec::new();
        let mut auxiliary_seen = false;
        for (slot, row) in rows.iter().enumerate() {
            assert_eq!(
                row.len(),
                6,
                "edge pair plus four momentum coefficients required"
            );
            let auxiliary = row[0] == 0 && row[1] == 0;
            assert!(
                !auxiliary_seen || auxiliary,
                "physical slots must precede auxiliary slots"
            );
            auxiliary_seen |= auxiliary;
            let mut momentum = Atom::Zero;
            let mut components = Vec::new();
            for (axis, coefficient) in row[2..].iter().enumerate() {
                if *coefficient != 0 {
                    momentum += Atom::num(*coefficient) * function!(S.k, Atom::num(axis + 1));
                    components.push(format!("({coefficient})*k{}", axis + 1));
                }
            }
            assert!(!momentum.is_zero(), "zero denominator momentum");
            propagators.push(TextPropagator {
                id: format!("D{}", slot + 1),
                expression: format!("({})^2-1", components.join("+")),
                target_power: i64::from(!auxiliary),
                power_shift: None,
            });
            if !auxiliary {
                physical_momenta.push(momentum.clone());
            }
            momenta.push(momentum);
            edges.push((row[0], row[1]));
        }
        let family = Compiler::new(Limits::default())
            .unwrap()
            .compile_text(TextProject {
                name: None,
                parameters: None,
                loop_momenta: (1..=4).map(|axis| format!("k{axis}")).collect(),
                external_momenta: Vec::new(),
                dimension: "d".into(),
                propagators,
                external_gram: Vec::new(),
                numerator: None,
            })
            .unwrap()
            .into_lowered(LoweringLimits::default())
            .unwrap()
            .into_family();
        Self {
            family: Arc::new(family),
            physical_momenta,
            momenta,
            edges,
        }
    }

    /// Convert the exact candidate key to a common-mass oracle integral.
    ///
    /// Vakint's legacy oracle requires a symbolic mass until after evaluation.
    /// Its caller subsequently sets `muvsq = 1`, matching the RustRed family.
    /// Negative powers are scalar numerator factors; no terminal basis is guessed.
    pub fn integral(&self, powers: &[i64]) -> Atom {
        assert_eq!(powers.len(), self.momenta.len());
        let mut propagators = Atom::num(1);
        let mut numerator = Atom::num(1);
        let mass_squared = vk_parse!("muvsq").unwrap();
        for (slot, (&power, momentum)) in powers.iter().zip(&self.momenta).enumerate() {
            if slot < self.physical_momenta.len() {
                let (left, right) = self.edges[slot];
                propagators *= vk_parse!(format!(
                    "prop({},edge({left},{right}),{},muvsq,{})",
                    slot + 1,
                    momentum.to_canonical_string(),
                    power.max(0),
                ))
                .unwrap();
            } else {
                assert!(power <= 0, "auxiliary denominator cannot become physical");
            }
            if power < 0 {
                let exponent = power.checked_neg().expect("bounded numerator power");
                numerator *=
                    (function!(S.dot, momentum, momentum) - &mass_squared).pow(Atom::num(exponent));
            }
        }
        numerator * function!(S.topo, propagators)
    }
}
