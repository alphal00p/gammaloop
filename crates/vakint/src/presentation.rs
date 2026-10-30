//! One read-only presentation model for terminal and notebook inspection.

use std::fmt;

use colored::Colorize;
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::{
    AvailablePrecision, CoefficientPrecision, ComponentPrecision, DataProvenance,
    MasterCoefficient, MasterSymbolicExpression, NumericalCoefficient, NumericalEvaluationResult,
    PrecisionSource, PrecisionWarning, VakintConstant, VakintMaster,
};

/// Expressions retain their native Symbolica renderer in notebook output.
pub(crate) enum Content {
    Text(String),
    Expression(Atom),
}

impl Content {
    /// Hide a single snapshot's transport tag only in the rendered formula.
    /// Keep mixed snapshots explicit so different data cannot appear to cancel.
    pub(crate) fn display_expression(value: &Atom) -> Atom {
        let identity = |part: AtomView<'_>| {
            let AtomView::Fun(function) = part else {
                return None;
            };
            let symbol = function.get_symbol();
            if symbol.get_namespace() != crate::NAMESPACE
                || symbol.get_stripped_name() != "RustRedConstant"
                || function.get_nargs() != 2
            {
                return None;
            }
            let mut arguments = function.iter();
            Some((
                arguments.next().unwrap().to_owned(),
                arguments.next().unwrap().to_owned(),
            ))
        };
        let mut snapshots = std::collections::BTreeSet::new();
        value.visitor(&mut |part| {
            if let Some((snapshot, _)) = identity(part) {
                snapshots.insert(snapshot);
                false
            } else {
                true
            }
        });
        if snapshots.len() != 1 {
            return value.clone();
        }
        value.replace_map(|part, _context, out| {
            if let Some((_, symbol)) = identity(part) {
                **out = symbol;
            }
        })
    }

    pub(crate) fn plain(&self) -> String {
        match self {
            Self::Text(value) => value.clone(),
            Self::Expression(value) => Self::display_expression(value).to_plain_string(),
        }
    }
}

impl From<String> for Content {
    fn from(value: String) -> Self {
        Self::Text(value)
    }
}

impl From<&str> for Content {
    fn from(value: &str) -> Self {
        Self::Text(value.into())
    }
}

#[derive(Default)]
pub(crate) struct Presentation {
    pub title: String,
    pub fields: Vec<(String, Content)>,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<Content>>,
    pub warnings: Vec<String>,
}

pub(crate) trait Present {
    fn presentation(&self) -> Presentation;
}

impl Presentation {
    pub(crate) fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            ..Self::default()
        }
    }

    pub(crate) fn field(&mut self, name: impl Into<String>, value: impl Into<Content>) {
        self.fields.push((name.into(), value.into()));
    }

    pub(crate) fn shortened(value: &str, limit: usize) -> String {
        let value = value.replace(['\n', '\r'], " ");
        if value.chars().count() <= limit {
            value
        } else {
            value
                .chars()
                .take(limit.saturating_sub(1))
                .collect::<String>()
                + "…"
        }
    }

    fn border(f: &mut fmt::Formatter<'_>, widths: &[usize], edges: [char; 3]) -> fmt::Result {
        write!(f, "{}", edges[0])?;
        for (index, width) in widths.iter().enumerate() {
            write!(f, "{}", "─".repeat(width + 2))?;
            write!(
                f,
                "{}",
                if index + 1 == widths.len() {
                    edges[2]
                } else {
                    edges[1]
                }
            )?;
        }
        writeln!(f)
    }
}

impl fmt::Display for Presentation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let color = std::env::var_os("NO_COLOR").is_none()
            && colored::control::SHOULD_COLORIZE.should_colorize();
        let title = if color {
            self.title.bold().cyan().to_string()
        } else {
            self.title.clone()
        };
        writeln!(f, "{title}")?;
        let label_width = self
            .fields
            .iter()
            .map(|(label, _)| label.chars().count())
            .max()
            .unwrap_or(0);
        let limit = if f.alternate() { usize::MAX } else { 84 };
        for (label, value) in &self.fields {
            let label = format!("{label:label_width$}");
            let label = if color {
                label.dimmed().to_string()
            } else {
                label
            };
            writeln!(f, "  {label}  {}", Self::shortened(&value.plain(), limit))?;
        }
        if !self.columns.is_empty() {
            let rows = std::iter::once(self.columns.clone())
                .chain(self.rows.iter().map(|row| {
                    row.iter()
                        .map(|value| {
                            Self::shortened(
                                &value.plain(),
                                if f.alternate() { usize::MAX } else { 36 },
                            )
                        })
                        .collect()
                }))
                .collect::<Vec<Vec<String>>>();
            let widths = (0..self.columns.len())
                .map(|index| {
                    rows.iter()
                        .map(|row| row.get(index).map_or(0, |cell| cell.chars().count()))
                        .max()
                        .unwrap_or(0)
                })
                .collect::<Vec<_>>();
            Self::border(f, &widths, ['┌', '┬', '┐'])?;
            for (index, row) in rows.iter().enumerate() {
                write!(f, "│")?;
                for (cell, width) in row.iter().zip(&widths) {
                    let padding = width.saturating_sub(cell.chars().count());
                    let cell = if color && index == 0 {
                        cell.bold().to_string()
                    } else {
                        cell.clone()
                    };
                    write!(f, " {cell}{} │", " ".repeat(padding))?;
                }
                writeln!(f)?;
                if index == 0 {
                    Self::border(f, &widths, ['├', '┼', '┤'])?;
                }
            }
            Self::border(f, &widths, ['└', '┴', '┘'])?;
        }
        for warning in &self.warnings {
            let heading = if color {
                "Warning".yellow().bold().to_string()
            } else {
                "Warning".into()
            };
            writeln!(f, "  {heading}: {warning}")?;
        }
        Ok(())
    }
}

pub(crate) fn digits(value: &ComponentPrecision) -> String {
    if value.exact_zero() {
        "exact zero".into()
    } else {
        value
            .significant_digits()
            .map_or_else(|| "unresolved".into(), |digits| format!("{digits} digits"))
    }
}

impl fmt::Display for AvailablePrecision {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ArbitraryPrecision => write!(f, "arbitrary precision"),
            Self::FixedPrecision(digits) => write!(f, "up to {digits} decimal digits"),
        }
    }
}

impl Present for AvailablePrecision {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Available precision");
        card.field("Capability", self.to_string());
        card
    }
}

impl Present for ComponentPrecision {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Component precision");
        card.field("Achieved", digits(self));
        card.field(
            "Absolute error",
            Content::Expression(Atom::num(self.absolute_error().clone())),
        );
        card.field("Error quality", format!("{:?}", self.quality()));
        card.field("Exact zero", self.exact_zero().to_string());
        card
    }
}

impl Present for CoefficientPrecision {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Coefficient precision");
        card.columns = ["Component", "Achieved", "Absolute error", "Error quality"]
            .map(String::from)
            .into();
        for (name, precision) in [("Real", self.real()), ("Imaginary", self.imaginary())] {
            card.rows.push(vec![
                name.into(),
                digits(precision).into(),
                Content::Expression(Atom::num(precision.absolute_error().clone())),
                format!("{:?}", precision.quality()).into(),
            ]);
        }
        card
    }
}

impl Present for PrecisionSource {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Precision source");
        card.field("Identity", self.identity.clone());
        if let Some(master) = &self.master_id {
            card.field("Master", master.clone());
        }
        if let Some(power) = self.epsilon_power {
            card.field("Epsilon power", power.to_string());
        }
        card.field(
            "Available digits",
            self.available_decimal_precision
                .map_or_else(|| "arbitrary".into(), |digits| digits.to_string()),
        );
        card.field("Source SHA-256", self.source_sha256.clone());
        card
    }
}

impl fmt::Display for PrecisionWarning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} [requested {} digits",
            self.message, self.requested_decimal_precision
        )?;
        if let Some(achieved) = self.achieved_decimal_precision {
            write!(f, "; achieved {achieved} digits")?;
        }
        write!(f, "; data {}]", self.data_version)?;
        if let Some(power) = self.epsilon_power {
            write!(f, "; ε^{power}")?;
        }
        for source in &self.sources {
            write!(f, "; {}", source.identity)?;
            if let Some(master) = &source.master_id {
                write!(f, " [master {master}]")?;
            }
            if let Some(power) = source.epsilon_power {
                write!(f, " ε^{power}")?;
            }
            if let Some(digits) = source.available_decimal_precision {
                write!(f, " ({digits} source digits)")?;
            }
        }
        Ok(())
    }
}

impl Present for PrecisionWarning {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Precision warning");
        card.field("Message", self.message.clone());
        card.field(
            "Requested digits",
            self.requested_decimal_precision.to_string(),
        );
        if let Some(achieved) = self.achieved_decimal_precision {
            card.field("Achieved digits", achieved.to_string());
        }
        card.field("Data version", self.data_version.clone());
        if let Some(power) = self.epsilon_power {
            card.field("Epsilon power", power.to_string());
        }
        card.columns = ["Source", "Master", "ε power", "Available digits"]
            .map(String::from)
            .into();
        for source in &self.sources {
            card.rows.push(vec![
                source.identity.clone().into(),
                source.master_id.clone().unwrap_or_default().into(),
                source
                    .epsilon_power
                    .map_or_else(String::new, |power| power.to_string())
                    .into(),
                source
                    .available_decimal_precision
                    .map_or_else(|| "arbitrary".into(), |digits| digits.to_string())
                    .into(),
            ]);
        }
        card
    }
}

impl Present for DataProvenance {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Data provenance");
        card.field("Data version", self.version());
        card.field("Source SHA-256", self.source_sha256());
        if let Some(master) = self.master_id() {
            card.field("Master", master);
        }
        if let Some(power) = self.epsilon_power() {
            card.field("Epsilon power", power.to_string());
        }
        card
    }
}

impl Present for NumericalCoefficient {
    fn presentation(&self) -> Presentation {
        let mut card = self.precision().presentation();
        card.title = "Numerical coefficient".into();
        card.field(
            "Central value",
            Content::Expression(Atom::num(self.value().clone())),
        );
        card
    }
}

impl Present for NumericalEvaluationResult {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Vakint numerical result");
        card.field(
            "Requested precision",
            format!("{} decimal digits", self.requested_decimal_precision()),
        );
        card.columns = [
            "ε power",
            "Central value",
            "Real accuracy",
            "Imaginary accuracy",
        ]
        .map(String::from)
        .into();
        for (power, coefficient) in self.coefficients() {
            card.rows.push(vec![
                power.to_string().into(),
                Content::Expression(Atom::num(coefficient.value().clone())),
                digits(coefficient.precision().real()).into(),
                digits(coefficient.precision().imaginary()).into(),
            ]);
        }
        card.warnings = self.warnings().iter().map(ToString::to_string).collect();
        card
    }
}

impl Present for VakintConstant {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Vakint constant");
        card.field(
            "Requested precision",
            format!("{} decimal digits", self.requested_decimal_precision()),
        );
        card.field(
            "Symbol",
            Content::Expression(self.symbol().as_atom().clone()),
        );
        card.field("Definition", Content::Expression(self.definition().clone()));
        card.field(
            "Available precision",
            self.available_precision().to_string(),
        );
        card.field(
            "Numerical value",
            Content::Expression(Atom::num(self.numerical_value().clone())),
        );
        card.field("Achieved precision", digits(self.precision().real()));
        card.field(
            "Error quality",
            format!("{:?}", self.precision().real().quality()),
        );
        card.field(
            "Absolute error",
            Content::Expression(Atom::num(self.precision().real().absolute_error().clone())),
        );
        card.field("Data version", self.provenance().version());
        card.field("Source SHA-256", self.provenance().source_sha256());
        card.warnings = self.warnings().iter().map(ToString::to_string).collect();
        card
    }
}

impl Present for MasterSymbolicExpression {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Master symbolic expression");
        match self {
            Self::Expression(value) => card.field("Expression", Content::Expression(value.clone())),
            Self::NumericalOnly(digits) => {
                card.field("Numerical only", format!("up to {digits} decimal digits"))
            }
        }
        card
    }
}

impl Present for MasterCoefficient {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new("Master coefficient");
        card.field(
            "Requested precision",
            format!("{} decimal digits", self.requested_decimal_precision()),
        );
        card.field("Symbol", Content::Expression(self.symbol().clone()));
        card.field("Epsilon power", self.epsilon_power().to_string());
        card.field(
            "Available precision",
            self.available_precision().to_string(),
        );
        card.field(
            "Numerical value",
            Content::Expression(Atom::num(self.numerical_value().clone())),
        );
        card.field("Achieved precision", digits(self.precision().real()));
        card.field(
            "Error quality",
            format!("{:?}", self.precision().real().quality()),
        );
        card.field(
            "Absolute error",
            Content::Expression(Atom::num(self.precision().real().absolute_error().clone())),
        );
        card.fields
            .extend(self.symbolic_expression().presentation().fields);
        card.field("Data version", self.provenance().version());
        card.warnings = self.warnings().iter().map(ToString::to_string).collect();
        card
    }
}

impl Present for VakintMaster {
    fn presentation(&self) -> Presentation {
        let mut card = Presentation::new(format!(
            "Vakint master · {} loop{}",
            self.loop_count(),
            if self.loop_count() == 1 { "" } else { "s" }
        ));
        card.field(
            "Symbol",
            Content::Expression(self.symbol().as_atom().clone()),
        );
        card.field("Master identity", self.master_id());
        card.field(
            "Requested precision",
            format!("{} decimal digits", self.requested_decimal_precision()),
        );
        card.field("Family", self.family_sha256());
        card.field(
            "Status",
            if self.is_computed() {
                "computed coefficients available"
            } else {
                "declared terminal · coefficients uncomputed"
            },
        );
        card.field("Normalization", self.normalization());
        card.field("Data version", self.provenance().version());
        let vertices = self.graph().vertices().count();
        let edges = self.graph().edges().count();
        card.field(
            "Graph",
            format!(
                "{vertices} {} · {edges} edge{}",
                if vertices == 1 { "vertex" } else { "vertices" },
                if edges == 1 { "" } else { "s" }
            ),
        );
        card.columns = [
            "ε power",
            "Numerical value",
            "Available precision",
            "Achieved",
            "Symbolic expression",
        ]
        .map(String::from)
        .into();
        for (power, coefficient) in self.coefficients() {
            let formula = match coefficient.symbolic_expression() {
                MasterSymbolicExpression::Expression(value) => Content::Expression(value.clone()),
                MasterSymbolicExpression::NumericalOnly(digits) => {
                    format!("numerical only · {digits} digits").into()
                }
            };
            card.rows.push(vec![
                power.to_string().into(),
                Content::Expression(Atom::num(coefficient.numerical_value().clone())),
                coefficient.available_precision().to_string().into(),
                digits(coefficient.precision().real()).into(),
                formula,
            ]);
            card.warnings
                .extend(coefficient.warnings().iter().map(ToString::to_string));
        }
        card.warnings.sort();
        card.warnings.dedup();
        card
    }
}

macro_rules! display_record {
    ($($ty:ty),* $(,)?) => { $(
        impl fmt::Display for $ty {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt::Display::fmt(&self.presentation(), f)
            }
        }
    )* };
}
display_record!(
    DataProvenance,
    NumericalCoefficient,
    NumericalEvaluationResult,
    VakintConstant,
    MasterCoefficient,
    MasterSymbolicExpression,
    VakintMaster
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_values_are_bounded_without_losing_alternate_detail() {
        let mut card = Presentation::new("Numerical result");
        card.field("Value", "1234567890".repeat(20));
        assert!(format!("{card}").contains('…'));
        assert!(format!("{card:#}").contains(&"1234567890".repeat(20)));
    }

    #[test]
    fn bound_constants_render_mathematically_without_changing_stored_identity() {
        let value = crate::vakint_parse!("2*RustRedConstant(snapshot_a,K5Moment(1))").unwrap();
        let content = Content::Expression(value.clone());
        assert_eq!(
            Content::display_expression(&value),
            crate::vakint_parse!("2*K5Moment(1)").unwrap()
        );
        assert!(!content.plain().contains("RustRedConstant"));
        let Content::Expression(stored) = content else {
            unreachable!()
        };
        assert_eq!(stored, value);
        assert!(stored.contains_symbol(crate::vakint_symbol!("RustRedConstant")));
        let mixed = crate::vakint_parse!(
            "RustRedConstant(snapshot_a,K5Moment(1))-RustRedConstant(snapshot_b,K5Moment(1))"
        )
        .unwrap();
        assert_eq!(Content::display_expression(&mixed), mixed);
    }
}
