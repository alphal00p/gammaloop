use std::{collections::BTreeMap, fs, path::Path};

use serde_json::Value;
use thiserror::Error;

use crate::{
    OrientationID, ParsedGraph, ThreeDExpression,
    surface::{HybridSurfaceID, LinearEnergyExpr, RationalAtomExt},
    tree::{NodeId, Tree},
};

#[derive(Debug, Error)]
pub enum EvaluationError {
    #[error("invalid JSON: {0}")]
    Json(#[from] serde_json::Error),
    #[error("{0}")]
    GraphIo(#[from] crate::graph_io::GraphIoError),
    #[error("failed to read `{path}`: {source}")]
    Read {
        path: String,
        source: std::io::Error,
    },
    #[error("invalid numeric value `{0}`")]
    InvalidNumber(String),
    #[error("invalid vector JSON; expected [[...], ...]")]
    InvalidVectorJson,
    #[error("mass map must be a JSON object")]
    InvalidMassMap,
    #[error("missing mass value for `{0}`")]
    MissingMass(String),
    #[error("uniform sampling scale M is required by this expression")]
    MissingUniformScale,
    #[error("uniform sampling scale M must be nonzero")]
    ZeroUniformScale,
    #[error("numerator parse error: {0}")]
    NumeratorParse(String),
    #[error("numerator evaluation error: {0}")]
    NumeratorEval(String),
    #[error(
        "residual four-dimensional denominator edge {edge_id} is outside the generated edge-energy map"
    )]
    ResidualDenominatorEdgeOutOfRange { edge_id: usize },
    #[error("residual four-dimensional denominator power {power} is outside the supported range")]
    ResidualDenominatorPowerOutOfRange { power: usize },
}

pub type Result<T> = std::result::Result<T, EvaluationError>;

#[derive(Debug, Clone)]
pub struct EvaluationInput {
    pub external_momenta: Vec<[f64; 4]>,
    pub loop_spatial_momenta: Vec<[f64; 3]>,
    pub masses: Vec<f64>,
    pub uniform_scale: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct EvaluationResult {
    pub value: f64,
    pub numerator_calls: usize,
}

impl EvaluationInput {
    pub fn deterministic(
        parsed: &ParsedGraph,
        seed: u64,
        mass_overrides: &BTreeMap<String, f64>,
        uniform_scale: Option<f64>,
    ) -> Result<Self> {
        let mut rng = DeterministicRng::new(seed);
        let external_momenta = (0..parsed.external_names.len())
            .map(|_| {
                [
                    rng.uniform(-0.7, 0.7),
                    rng.uniform(-0.3, 0.3),
                    rng.uniform(-0.3, 0.3),
                    rng.uniform(-0.3, 0.3),
                ]
            })
            .collect::<Vec<_>>();
        let loop_spatial_momenta = (0..parsed.loop_names.len())
            .map(|_| {
                [
                    rng.uniform(-0.25, 0.25),
                    rng.uniform(-0.25, 0.25),
                    rng.uniform(-0.25, 0.25),
                ]
            })
            .collect::<Vec<_>>();
        let mut generated_masses = BTreeMap::<String, f64>::new();
        let masses = parsed
            .internal_edges
            .iter()
            .map(|edge| {
                if let Some(key) = &edge.mass_key {
                    mass_overrides.get(key).copied().unwrap_or_else(|| {
                        *generated_masses
                            .entry(key.clone())
                            .or_insert_with(|| round_to_6(rng.uniform(0.4, 1.3)))
                    })
                } else {
                    0.0
                }
            })
            .collect::<Vec<_>>();
        Ok(Self {
            external_momenta,
            loop_spatial_momenta,
            masses,
            uniform_scale,
        })
    }

    pub fn with_overrides(
        mut self,
        external_momenta: Option<Vec<[f64; 4]>>,
        loop_spatial_momenta: Option<Vec<[f64; 3]>>,
        mass_overrides: &BTreeMap<String, f64>,
        parsed: &ParsedGraph,
    ) -> Result<Self> {
        if let Some(external_momenta) = external_momenta {
            self.external_momenta = external_momenta;
        }
        if let Some(loop_spatial_momenta) = loop_spatial_momenta {
            self.loop_spatial_momenta = loop_spatial_momenta;
        }
        for (edge_id, edge) in parsed.internal_edges.iter().enumerate() {
            if let Some(key) = &edge.mass_key
                && let Some(value) = mass_overrides.get(key)
            {
                self.masses[edge_id] = *value;
            }
        }
        Ok(self)
    }
}
pub fn load_expression_json(path_or_json: &str) -> Result<ThreeDExpression<OrientationID>> {
    let text = if Path::new(path_or_json).is_file() {
        fs::read_to_string(path_or_json).map_err(|source| EvaluationError::Read {
            path: path_or_json.to_string(),
            source,
        })?
    } else {
        path_or_json.to_string()
    };
    let value = serde_json::from_str::<Value>(&text)?;
    if let Some(expression) = value.get("expression") {
        Ok(serde_json::from_value(expression.clone())?)
    } else {
        Ok(serde_json::from_value(value)?)
    }
}

pub fn parse_vectors4(value: Option<&str>) -> Result<Option<Vec<[f64; 4]>>> {
    let Some(value) = value else {
        return Ok(None);
    };
    parse_vectors::<4>(value).map(Some)
}

pub fn parse_vectors3(value: Option<&str>) -> Result<Option<Vec<[f64; 3]>>> {
    let Some(value) = value else {
        return Ok(None);
    };
    parse_vectors::<3>(value).map(Some)
}

pub fn parse_mass_map(value: Option<&str>, file: Option<&Path>) -> Result<BTreeMap<String, f64>> {
    let text = if let Some(path) = file {
        Some(
            fs::read_to_string(path).map_err(|source| EvaluationError::Read {
                path: path.display().to_string(),
                source,
            })?,
        )
    } else {
        value.map(str::to_string)
    };
    let Some(text) = text else {
        return Ok(BTreeMap::new());
    };
    let value = serde_json::from_str::<Value>(&text)?;
    let object = value.as_object().ok_or(EvaluationError::InvalidMassMap)?;
    object
        .iter()
        .map(|(key, value)| Ok((key.clone(), value_to_f64(value)?)))
        .collect()
}

pub fn parse_uniform_scale(value: Option<&str>) -> Result<Option<f64>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let parsed = value
        .parse::<f64>()
        .map_err(|_| EvaluationError::InvalidNumber(value.to_string()))?;
    if parsed == 0.0 {
        return Err(EvaluationError::ZeroUniformScale);
    }
    Ok(Some(parsed))
}

pub fn evaluate_expression(
    parsed: &ParsedGraph,
    expression: &ThreeDExpression<OrientationID>,
    numerator_expr: &str,
    input: &EvaluationInput,
) -> Result<EvaluationResult> {
    let numerator = NumeratorExpr::parse(numerator_expr)?;
    let evaluator = ExpressionEvaluator::new(parsed, expression, input);
    evaluator.evaluate(&numerator)
}

struct ExpressionEvaluator<'a> {
    parsed: &'a ParsedGraph,
    expression: &'a ThreeDExpression<OrientationID>,
    input: &'a EvaluationInput,
    internal_energies: Vec<f64>,
    external_energies: Vec<f64>,
}

impl<'a> ExpressionEvaluator<'a> {
    fn new(
        parsed: &'a ParsedGraph,
        expression: &'a ThreeDExpression<OrientationID>,
        input: &'a EvaluationInput,
    ) -> Self {
        let internal_energies = parsed
            .internal_edges
            .iter()
            .enumerate()
            .map(|(edge_id, edge)| {
                let spatial = edge_spatial_momentum(
                    &edge.signature,
                    &input.loop_spatial_momenta,
                    &input.external_momenta,
                );
                (spatial.iter().map(|x| x * x).sum::<f64>() + input.masses[edge_id].powi(2)).sqrt()
            })
            .collect();
        let external_energies = input
            .external_momenta
            .iter()
            .map(|momentum| momentum[0])
            .collect();
        Self {
            parsed,
            expression,
            input,
            internal_energies,
            external_energies,
        }
    }

    fn evaluate(&self, numerator: &NumeratorExpr) -> Result<EvaluationResult> {
        if self.expression_uses_uniform_scale() && self.input.uniform_scale.is_none() {
            return Err(EvaluationError::MissingUniformScale);
        }
        let mut total = 0.0;
        let mut numerator_cache = BTreeMap::<String, f64>::new();
        let mut numerator_calls = 0usize;
        for orientation in &self.expression.orientations {
            let residual_denominator_factor =
                self.expression.residual_denominators.iter().try_fold(
                    1.0,
                    |factor, denominator| {
                        let edge_id = denominator.edge_id.0;
                        let energy_map = orientation.edge_energy_map.get(edge_id).ok_or(
                            EvaluationError::ResidualDenominatorEdgeOutOfRange { edge_id },
                        )?;
                        let on_shell_energy = self.internal_energies.get(edge_id).ok_or(
                            EvaluationError::ResidualDenominatorEdgeOutOfRange { edge_id },
                        )?;
                        let power = i32::try_from(denominator.power).map_err(|_| {
                            EvaluationError::ResidualDenominatorPowerOutOfRange {
                                power: denominator.power,
                            }
                        })?;
                        let q0 = self.linear_expr_value(energy_map)?;
                        Ok::<_, EvaluationError>(
                            factor / (q0 * q0 - on_shell_energy * on_shell_energy).powi(power),
                        )
                    },
                )?;
            for variant in &orientation.variants {
                let mut denom = residual_denominator_factor * rational_to_f64(&variant.prefactor);
                for edge in &variant.half_edges {
                    denom /= 2.0 * self.internal_energies[edge.0];
                }
                if variant.uniform_scale_power > 0 {
                    let scale = self
                        .input
                        .uniform_scale
                        .ok_or(EvaluationError::MissingUniformScale)?;
                    denom /= scale.powi(variant.uniform_scale_power as i32);
                }
                denom *= self.tree_sum(&variant.denominator)?;

                let mut numerator_surface_factor = 1.0;
                for surface_id in &variant.numerator_surfaces {
                    numerator_surface_factor *= self.surface_value(*surface_id)?;
                }

                let num_key = format!(
                    "{:?}|{:?}",
                    orientation.loop_energy_map, orientation.edge_energy_map
                );
                let numerator_value = if let Some(value) = numerator_cache.get(&num_key) {
                    *value
                } else {
                    numerator_calls += 1;
                    let loop_four = self.loop_four_vectors(
                        &orientation.loop_energy_map,
                        &orientation.edge_energy_map,
                    )?;
                    let edge_four = self.edge_four_vectors(&orientation.edge_energy_map)?;
                    let value = numerator.eval(&EvalContext {
                        loops: &loop_four,
                        edges: &edge_four,
                        external: &self.input.external_momenta,
                    })?;
                    numerator_cache.insert(num_key, value);
                    value
                };
                let contribution = denom * numerator_surface_factor * numerator_value;
                total += contribution;
            }
        }
        Ok(EvaluationResult {
            value: total,
            numerator_calls,
        })
    }

    fn expression_uses_uniform_scale(&self) -> bool {
        self.expression
            .surfaces
            .linear_surface_cache
            .iter()
            .any(|surface| surface.expression.uses_uniform_scale())
            || self.expression.orientations.iter().any(|orientation| {
                orientation
                    .loop_energy_map
                    .iter()
                    .chain(&orientation.edge_energy_map)
                    .any(LinearEnergyExpr::uses_uniform_scale)
                    || orientation
                        .variants
                        .iter()
                        .any(|variant| variant.uniform_scale_power > 0)
            })
    }

    fn tree_sum(&self, tree: &Tree<HybridSurfaceID>) -> Result<f64> {
        self.tree_node_value(tree, NodeId::root())
    }

    fn tree_node_value(&self, tree: &Tree<HybridSurfaceID>, node_id: NodeId) -> Result<f64> {
        let node = tree.get_node(node_id);
        let mut factor = match node.data {
            HybridSurfaceID::Unit => 1.0,
            surface_id => {
                let value = self.surface_value(surface_id)?;
                1.0 / if value == 0.0 { 1.0e-80 } else { value }
            }
        };
        if !node.children.is_empty() {
            factor *= node
                .children
                .iter()
                .map(|child| self.tree_node_value(tree, *child))
                .sum::<Result<f64>>()?;
        }
        Ok(factor)
    }

    fn surface_value(&self, surface_id: HybridSurfaceID) -> Result<f64> {
        match surface_id {
            HybridSurfaceID::Unit => Ok(1.0),
            HybridSurfaceID::Infinite => Ok(0.0),
            HybridSurfaceID::Linear(id) => self
                .linear_expr_value(&self.expression.surfaces.linear_surface_cache[id].expression),
            HybridSurfaceID::Esurface(_) | HybridSurfaceID::Hsurface(_) => {
                Err(EvaluationError::NumeratorEval(
                    "legacy non-linear surfaces are not supported by three-dimensional-reps evaluator"
                        .to_string(),
                ))
            }
        }
    }

    fn linear_expr_value(&self, expr: &LinearEnergyExpr) -> Result<f64> {
        let mut total = rational_to_f64(&expr.constant);
        for (edge_id, coeff) in &expr.internal_terms {
            total += rational_to_f64(coeff) * self.internal_energies[edge_id.0];
        }
        for (edge_id, coeff) in &expr.external_terms {
            total += rational_to_f64(coeff) * self.external_energies[edge_id.0];
        }
        if !expr.uniform_scale_coeff.is_zero_coeff() {
            total += rational_to_f64(&expr.uniform_scale_coeff)
                * self
                    .input
                    .uniform_scale
                    .ok_or(EvaluationError::MissingUniformScale)?;
        }
        Ok(total)
    }

    fn loop_carrier_edge(&self, loop_id: usize) -> Option<usize> {
        let loop_name = self.parsed.loop_names.get(loop_id)?;
        self.parsed
            .internal_edges
            .iter()
            .position(|edge| &edge.label == loop_name)
    }

    fn loop_four_vectors(
        &self,
        loop_energy_map: &[LinearEnergyExpr],
        edge_energy_map: &[LinearEnergyExpr],
    ) -> Result<Vec<[f64; 4]>> {
        loop_energy_map
            .iter()
            .enumerate()
            .map(|(loop_id, expr)| {
                let carrier_edge = self.loop_carrier_edge(loop_id);
                let spatial = if let Some(edge_id) = carrier_edge {
                    edge_spatial_momentum(
                        &self.parsed.internal_edges[edge_id].signature,
                        &self.input.loop_spatial_momenta,
                        &self.input.external_momenta,
                    )
                } else {
                    self.input.loop_spatial_momenta[loop_id]
                };
                let energy = if let Some(edge_id) = carrier_edge {
                    if let Some(edge_expr) = edge_energy_map.get(edge_id) {
                        self.linear_expr_value(edge_expr)?
                    } else {
                        self.linear_expr_value(expr)?
                    }
                } else {
                    self.linear_expr_value(expr)?
                };
                Ok([energy, spatial[0], spatial[1], spatial[2]])
            })
            .collect()
    }

    fn edge_four_vectors(&self, edge_energy_map: &[LinearEnergyExpr]) -> Result<Vec<[f64; 4]>> {
        edge_energy_map
            .iter()
            .enumerate()
            .map(|(edge_id, expr)| {
                let spatial = edge_spatial_momentum(
                    &self.parsed.internal_edges[edge_id].signature,
                    &self.input.loop_spatial_momenta,
                    &self.input.external_momenta,
                );
                Ok([
                    self.linear_expr_value(expr)?,
                    spatial[0],
                    spatial[1],
                    spatial[2],
                ])
            })
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq)]
enum EvalValue {
    Scalar(f64),
    Vector([f64; 4]),
}

impl EvalValue {
    fn scalar(self) -> Result<f64> {
        match self {
            Self::Scalar(value) => Ok(value),
            Self::Vector(_) => Err(EvaluationError::NumeratorEval(
                "expected scalar expression".to_string(),
            )),
        }
    }
}

struct EvalContext<'a> {
    loops: &'a [[f64; 4]],
    edges: &'a [[f64; 4]],
    external: &'a [[f64; 4]],
}

#[derive(Debug, Clone)]
struct NumeratorExpr {
    root: ExprNode,
}

impl NumeratorExpr {
    fn parse(input: &str) -> Result<Self> {
        let tokens = Lexer::new(input).tokens()?;
        let mut parser = Parser::new(tokens);
        let root = parser.parse_expression(0)?;
        parser.expect_end()?;
        Ok(Self { root })
    }

    fn eval(&self, context: &EvalContext<'_>) -> Result<f64> {
        self.root.eval(context)?.scalar()
    }
}

#[derive(Debug, Clone)]
enum ExprNode {
    Number(f64),
    VectorRef(VectorKind, isize),
    ComponentRef(VectorKind, isize, usize),
    Dot(Box<ExprNode>, Box<ExprNode>),
    Neg(Box<ExprNode>),
    Add(Box<ExprNode>, Box<ExprNode>),
    Sub(Box<ExprNode>, Box<ExprNode>),
    Mul(Box<ExprNode>, Box<ExprNode>),
    Div(Box<ExprNode>, Box<ExprNode>),
    Pow(Box<ExprNode>, i32),
}

impl ExprNode {
    fn eval(&self, context: &EvalContext<'_>) -> Result<EvalValue> {
        match self {
            Self::Number(value) => Ok(EvalValue::Scalar(*value)),
            Self::VectorRef(kind, index) => Ok(EvalValue::Vector(kind.vector(context, *index)?)),
            Self::ComponentRef(kind, index, component) => {
                Ok(EvalValue::Scalar(kind.vector(context, *index)?[*component]))
            }
            Self::Dot(lhs, rhs) => {
                let EvalValue::Vector(lhs) = lhs.eval(context)? else {
                    return Err(EvaluationError::NumeratorEval(
                        "dot expects vector arguments".to_string(),
                    ));
                };
                let EvalValue::Vector(rhs) = rhs.eval(context)? else {
                    return Err(EvaluationError::NumeratorEval(
                        "dot expects vector arguments".to_string(),
                    ));
                };
                Ok(EvalValue::Scalar(
                    lhs[0] * rhs[0] - lhs[1] * rhs[1] - lhs[2] * rhs[2] - lhs[3] * rhs[3],
                ))
            }
            Self::Neg(expr) => Ok(EvalValue::Scalar(-expr.eval(context)?.scalar()?)),
            Self::Add(lhs, rhs) => Ok(EvalValue::Scalar(
                lhs.eval(context)?.scalar()? + rhs.eval(context)?.scalar()?,
            )),
            Self::Sub(lhs, rhs) => Ok(EvalValue::Scalar(
                lhs.eval(context)?.scalar()? - rhs.eval(context)?.scalar()?,
            )),
            Self::Mul(lhs, rhs) => Ok(EvalValue::Scalar(
                lhs.eval(context)?.scalar()? * rhs.eval(context)?.scalar()?,
            )),
            Self::Div(lhs, rhs) => Ok(EvalValue::Scalar(
                lhs.eval(context)?.scalar()? / rhs.eval(context)?.scalar()?,
            )),
            Self::Pow(base, exponent) => Ok(EvalValue::Scalar(
                base.eval(context)?.scalar()?.powi(*exponent),
            )),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum VectorKind {
    Loops,
    Edges,
    External,
}

impl VectorKind {
    fn from_ident(value: &str) -> Option<Self> {
        match value {
            "loops" => Some(Self::Loops),
            "edges" => Some(Self::Edges),
            "ext" => Some(Self::External),
            _ => None,
        }
    }

    fn vector(self, context: &EvalContext<'_>, index: isize) -> Result<[f64; 4]> {
        let slice = match self {
            Self::Loops => context.loops,
            Self::Edges => context.edges,
            Self::External => context.external,
        };
        let index = normalize_index(index, slice.len())?;
        Ok(slice[index])
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Number(f64),
    Ident(String),
    Integer(isize),
    Plus,
    Minus,
    Star,
    Slash,
    Pow,
    LParen,
    RParen,
    LBracket,
    RBracket,
    Comma,
    End,
}

struct Lexer<'a> {
    input: &'a str,
    chars: std::iter::Peekable<std::str::CharIndices<'a>>,
}

impl<'a> Lexer<'a> {
    fn new(input: &'a str) -> Self {
        Self {
            input,
            chars: input.char_indices().peekable(),
        }
    }

    fn tokens(mut self) -> Result<Vec<Token>> {
        let mut tokens = Vec::new();
        while let Some((_, ch)) = self.chars.peek().copied() {
            match ch {
                ch if ch.is_whitespace() => {
                    self.chars.next();
                }
                '0'..='9' | '.' => tokens.push(self.number()?),
                'a'..='z' | 'A'..='Z' | '_' => tokens.push(self.ident()),
                '+' => {
                    self.chars.next();
                    tokens.push(Token::Plus);
                }
                '-' => {
                    self.chars.next();
                    tokens.push(Token::Minus);
                }
                '*' => {
                    self.chars.next();
                    if self.chars.peek().is_some_and(|(_, next)| *next == '*') {
                        self.chars.next();
                        tokens.push(Token::Pow);
                    } else {
                        tokens.push(Token::Star);
                    }
                }
                '/' => {
                    self.chars.next();
                    tokens.push(Token::Slash);
                }
                '(' => {
                    self.chars.next();
                    tokens.push(Token::LParen);
                }
                ')' => {
                    self.chars.next();
                    tokens.push(Token::RParen);
                }
                '[' => {
                    self.chars.next();
                    tokens.push(Token::LBracket);
                }
                ']' => {
                    self.chars.next();
                    tokens.push(Token::RBracket);
                }
                ',' => {
                    self.chars.next();
                    tokens.push(Token::Comma);
                }
                _ => {
                    return Err(EvaluationError::NumeratorParse(format!(
                        "unexpected character `{ch}`"
                    )));
                }
            }
        }
        tokens.push(Token::End);
        Ok(tokens)
    }

    fn number(&mut self) -> Result<Token> {
        let start = self.chars.peek().map(|(idx, _)| *idx).unwrap_or_default();
        let mut end = start;
        let mut has_dot = false;
        while let Some((idx, ch)) = self.chars.peek().copied() {
            if ch.is_ascii_digit() {
                end = idx + ch.len_utf8();
                self.chars.next();
            } else if ch == '.' && !has_dot {
                has_dot = true;
                end = idx + 1;
                self.chars.next();
            } else {
                break;
            }
        }
        let text = &self.input[start..end];
        let value = text
            .parse::<f64>()
            .map_err(|_| EvaluationError::InvalidNumber(text.to_string()))?;
        if text.chars().all(|ch| ch.is_ascii_digit())
            && let Ok(integer) = text.parse::<isize>()
        {
            Ok(Token::Integer(integer))
        } else {
            Ok(Token::Number(value))
        }
    }

    fn ident(&mut self) -> Token {
        let start = self.chars.peek().map(|(idx, _)| *idx).unwrap_or_default();
        let mut end = start;
        while let Some((idx, ch)) = self.chars.peek().copied() {
            if ch.is_ascii_alphanumeric() || ch == '_' {
                end = idx + ch.len_utf8();
                self.chars.next();
            } else {
                break;
            }
        }
        Token::Ident(self.input[start..end].to_string())
    }
}

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn new(tokens: Vec<Token>) -> Self {
        Self { tokens, pos: 0 }
    }

    fn parse_expression(&mut self, min_bp: u8) -> Result<ExprNode> {
        let mut lhs = self.parse_prefix()?;
        loop {
            let op = self.peek().clone();
            let Some((left_bp, right_bp)) = infix_binding_power(&op) else {
                break;
            };
            if left_bp < min_bp {
                break;
            }
            self.next();
            if matches!(op, Token::Pow) {
                let exponent = self.parse_signed_integer()?;
                lhs = ExprNode::Pow(Box::new(lhs), exponent);
                continue;
            }
            let rhs = self.parse_expression(right_bp)?;
            lhs = match op {
                Token::Plus => ExprNode::Add(Box::new(lhs), Box::new(rhs)),
                Token::Minus => ExprNode::Sub(Box::new(lhs), Box::new(rhs)),
                Token::Star => ExprNode::Mul(Box::new(lhs), Box::new(rhs)),
                Token::Slash => ExprNode::Div(Box::new(lhs), Box::new(rhs)),
                _ => unreachable!("binding power only returned arithmetic operators"),
            };
        }
        Ok(lhs)
    }

    fn parse_prefix(&mut self) -> Result<ExprNode> {
        match self.next().clone() {
            Token::Number(value) => Ok(ExprNode::Number(value)),
            Token::Integer(value) => Ok(ExprNode::Number(value as f64)),
            Token::Minus => Ok(ExprNode::Neg(Box::new(self.parse_expression(9)?))),
            Token::Plus => self.parse_expression(9),
            Token::LParen => {
                let expr = self.parse_expression(0)?;
                self.expect(Token::RParen)?;
                Ok(expr)
            }
            Token::Ident(ident) if ident == "dot" => {
                self.expect(Token::LParen)?;
                let lhs = self.parse_expression(0)?;
                self.expect(Token::Comma)?;
                let rhs = self.parse_expression(0)?;
                self.expect(Token::RParen)?;
                Ok(ExprNode::Dot(Box::new(lhs), Box::new(rhs)))
            }
            Token::Ident(ident) => self.parse_reference(&ident),
            other => Err(EvaluationError::NumeratorParse(format!(
                "unexpected token `{other:?}`"
            ))),
        }
    }

    fn parse_reference(&mut self, ident: &str) -> Result<ExprNode> {
        let kind = VectorKind::from_ident(ident).ok_or_else(|| {
            EvaluationError::NumeratorParse(format!("unknown identifier `{ident}`"))
        })?;
        self.expect(Token::LBracket)?;
        let vector_index = self.parse_signed_isize()?;
        self.expect(Token::RBracket)?;
        if !matches!(self.peek(), Token::LBracket) {
            return Ok(ExprNode::VectorRef(kind, vector_index));
        }
        self.expect(Token::LBracket)?;
        let component = self.parse_signed_isize()?;
        self.expect(Token::RBracket)?;
        if !(0..=3).contains(&component) {
            return Err(EvaluationError::NumeratorParse(format!(
                "four-vector component index {component} is out of range"
            )));
        }
        Ok(ExprNode::ComponentRef(
            kind,
            vector_index,
            component as usize,
        ))
    }

    fn parse_signed_integer(&mut self) -> Result<i32> {
        let sign = if matches!(self.peek(), Token::Minus) {
            self.next();
            -1
        } else {
            if matches!(self.peek(), Token::Plus) {
                self.next();
            }
            1
        };
        let token = self.next().clone();
        match token {
            Token::Integer(value) => Ok(sign * value as i32),
            other => Err(EvaluationError::NumeratorParse(format!(
                "expected integer exponent, got `{other:?}`"
            ))),
        }
    }

    fn parse_signed_isize(&mut self) -> Result<isize> {
        let sign = if matches!(self.peek(), Token::Minus) {
            self.next();
            -1
        } else {
            if matches!(self.peek(), Token::Plus) {
                self.next();
            }
            1
        };
        let token = self.next().clone();
        match token {
            Token::Integer(value) => Ok(sign * value),
            other => Err(EvaluationError::NumeratorParse(format!(
                "expected integer index, got `{other:?}`"
            ))),
        }
    }

    fn expect(&mut self, expected: Token) -> Result<()> {
        let token = self.next().clone();
        if std::mem::discriminant(&token) == std::mem::discriminant(&expected) {
            Ok(())
        } else {
            Err(EvaluationError::NumeratorParse(format!(
                "expected `{expected:?}`, got `{token:?}`"
            )))
        }
    }

    fn expect_end(&self) -> Result<()> {
        if matches!(self.peek(), Token::End) {
            Ok(())
        } else {
            Err(EvaluationError::NumeratorParse(format!(
                "unexpected trailing token `{:?}`",
                self.peek()
            )))
        }
    }

    fn peek(&self) -> &Token {
        &self.tokens[self.pos]
    }

    fn next(&mut self) -> &Token {
        let pos = self.pos;
        self.pos += 1;
        &self.tokens[pos]
    }
}

fn infix_binding_power(token: &Token) -> Option<(u8, u8)> {
    match token {
        Token::Plus | Token::Minus => Some((1, 2)),
        Token::Star | Token::Slash => Some((3, 4)),
        Token::Pow => Some((7, 6)),
        _ => None,
    }
}

fn parse_vectors<const N: usize>(value: &str) -> Result<Vec<[f64; N]>> {
    let value = if Path::new(value).is_file() {
        fs::read_to_string(value).map_err(|source| EvaluationError::Read {
            path: value.to_string(),
            source,
        })?
    } else {
        value.to_string()
    };
    let value = serde_json::from_str::<Value>(&value)?;
    let rows = value.as_array().ok_or(EvaluationError::InvalidVectorJson)?;
    rows.iter()
        .map(|row| {
            let values = row.as_array().ok_or(EvaluationError::InvalidVectorJson)?;
            if values.len() != N {
                return Err(EvaluationError::InvalidVectorJson);
            }
            let mut out = [0.0; N];
            for (idx, value) in values.iter().enumerate() {
                out[idx] = value_to_f64(value)?;
            }
            Ok(out)
        })
        .collect()
}

fn value_to_f64(value: &Value) -> Result<f64> {
    if let Some(value) = value.as_f64() {
        return Ok(value);
    }
    if let Some(value) = value.as_str() {
        return value
            .parse::<f64>()
            .map_err(|_| EvaluationError::InvalidNumber(value.to_string()));
    }
    Err(EvaluationError::InvalidNumber(value.to_string()))
}

fn edge_spatial_momentum(
    signature: &crate::MomentumSignature,
    loop_spatial_momenta: &[[f64; 3]],
    external_momenta: &[[f64; 4]],
) -> [f64; 3] {
    let mut out = [0.0; 3];
    for (loop_id, coeff) in signature.loop_signature.iter().enumerate() {
        for (dim, item) in out.iter_mut().enumerate() {
            *item += *coeff as f64 * loop_spatial_momenta[loop_id][dim];
        }
    }
    for (external_id, coeff) in signature.external_signature.iter().enumerate() {
        for (dim, item) in out.iter_mut().enumerate() {
            *item += *coeff as f64 * external_momenta[external_id][dim + 1];
        }
    }
    out
}

fn rational_to_f64(value: &symbolica::atom::Atom) -> f64 {
    value.rational_coeff().to_f64()
}

fn normalize_index(index: isize, len: usize) -> Result<usize> {
    let resolved = if index < 0 {
        len as isize + index
    } else {
        index
    };
    if resolved < 0 || resolved as usize >= len {
        return Err(EvaluationError::NumeratorEval(format!(
            "index {index} is out of range for length {len}"
        )));
    }
    Ok(resolved as usize)
}

fn round_to_6(value: f64) -> f64 {
    (value * 1_000_000.0).round() / 1_000_000.0
}

#[derive(Debug, Clone)]
struct DeterministicRng {
    state: u64,
}

impl DeterministicRng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    fn uniform(&mut self, low: f64, high: f64) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let unit = ((self.state >> 11) as f64) / ((1_u64 << 53) as f64);
        low + (high - low) * unit
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::{
        Generate3DExpressionOptions, MomentumSignature, NumeratorSamplingScaleMode,
        generate_3d_expression,
        generation::generate_3d_expression_from_parsed,
        graph_io::{ParsedGraphInitialStateCutEdge, ParsedGraphInternalEdge},
    };

    use super::*;

    #[test]
    fn evaluator_includes_projected_residual_tree_denominators() {
        let parsed = crate::graph_io::test_graphs::pure_tree_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                preserve_internal_edges_as_four_d_denominators: vec![0, 1],
                ..Default::default()
            },
        )
        .unwrap();
        let input = EvaluationInput {
            external_momenta: vec![[3.0, 0.1, 0.2, 0.3], [2.0, -0.2, 0.1, 0.05]],
            loop_spatial_momenta: Vec::new(),
            masses: vec![0.7, 1.1],
            uniform_scale: None,
        };

        let value = evaluate_expression(&parsed, &expression, "1", &input)
            .unwrap()
            .value;
        let edge_0_energy_squared =
            0.1_f64.powi(2) + 0.2_f64.powi(2) + 0.3_f64.powi(2) + 0.7_f64.powi(2);
        let edge_1_energy_squared =
            (-0.1_f64).powi(2) + 0.3_f64.powi(2) + 0.35_f64.powi(2) + 1.1_f64.powi(2);
        let expected = 1.0
            / ((3.0_f64.powi(2) - edge_0_energy_squared)
                * (5.0_f64.powi(2) - edge_1_energy_squared));

        assert!((value - expected).abs() < 1.0e-13);
    }

    #[test]
    fn numerator_parser_evaluates_edge_energy_and_dot_products() {
        let expr = NumeratorExpr::parse("edges[0][0]**2 + dot(edges[-1], ext[0])").unwrap();
        let value = expr
            .eval(&EvalContext {
                loops: &[],
                edges: &[[2.0, 1.0, 0.0, 0.0], [3.0, 0.0, 1.0, 0.0]],
                external: &[[5.0, 1.0, 1.0, 1.0]],
            })
            .unwrap();
        assert_eq!(value, 18.0);
    }

    #[test]
    fn loop_energy_uses_named_emr_carrier_before_diagnostic_fallback() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        let expression = ThreeDExpression::new_empty();
        let input = EvaluationInput {
            external_momenta: vec![[0.0; 4]; 3],
            loop_spatial_momenta: vec![[0.0; 3]],
            masses: vec![2.0, 3.0, 4.0, 5.0],
            uniform_scale: None,
        };
        let loop_energy_map = vec![LinearEnergyExpr::ose(
            linnet::half_edge::involution::EdgeIndex(0),
            1,
        )];
        let edge_energy_map = vec![
            LinearEnergyExpr::ose(linnet::half_edge::involution::EdgeIndex(1), 1),
            LinearEnergyExpr::zero(),
            LinearEnergyExpr::zero(),
            LinearEnergyExpr::zero(),
        ];

        let without_carrier = ExpressionEvaluator::new(&parsed, &expression, &input)
            .loop_four_vectors(&loop_energy_map, &edge_energy_map)
            .unwrap();
        assert_eq!(without_carrier[0][0], 2.0);

        parsed.internal_edges[0].label = parsed.loop_names[0].clone();
        let with_carrier = ExpressionEvaluator::new(&parsed, &expression, &input)
            .loop_four_vectors(&loop_energy_map, &edge_energy_map)
            .unwrap();
        assert_eq!(with_carrier[0][0], 3.0);
    }

    #[test]
    fn evaluator_requires_nonzero_uniform_scale_when_expression_uses_m() {
        let parsed = crate::graph_io::test_graphs::box_pow3_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(3, 4)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                ..Default::default()
            },
        )
        .unwrap();
        let input = EvaluationInput::deterministic(&parsed, 1337, &BTreeMap::new(), None).unwrap();

        assert!(matches!(
            evaluate_expression(&parsed, &expression, "edges[3][0]**4", &input),
            Err(EvaluationError::MissingUniformScale)
        ));
        assert!(matches!(
            parse_uniform_scale(Some("0")),
            Err(EvaluationError::ZeroUniformScale)
        ));
    }

    #[test]
    fn initial_cut_repeated_spectator_quadratic_matches_finite_pole_identity() {
        let internal_edges = [
            (0, 1, 0, vec![0, 0, 0], vec![1], "m0"),
            (1, 0, 1, vec![1, 0, 0], vec![0], "m0"),
            (2, 0, 5, vec![-1, 0, 0], vec![1], "m0"),
            (3, 1, 3, vec![1, 0, 0], vec![-1], "m0"),
            (4, 2, 3, vec![0, 1, 0], vec![0], "m0"),
            (5, 2, 3, vec![-1, -1, 0], vec![1], "m0"),
            (6, 2, 4, vec![1, 0, 0], vec![-1], "m0"),
            (7, 4, 5, vec![0, 0, 1], vec![0], "m0"),
            (8, 4, 5, vec![1, 0, -1], vec![-1], "m1"),
        ]
        .into_iter()
        .map(
            |(edge_id, tail, head, loop_signature, external_signature, mass_key)| {
                ParsedGraphInternalEdge {
                    edge_id,
                    tail,
                    head,
                    label: format!("q{edge_id}"),
                    mass_key: Some(mass_key.to_string()),
                    signature: MomentumSignature {
                        loop_signature,
                        external_signature,
                    },
                    had_pow: false,
                }
            },
        )
        .collect();
        let parsed = ParsedGraph {
            internal_edges,
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 0,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["e1".to_string(), "e4".to_string(), "e7".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..6).map(|node| (format!("v{node}"), node)).collect(),
        };
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(7, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let input = EvaluationInput::deterministic(
            &parsed,
            17,
            &BTreeMap::from([("m0".to_string(), 0.0), ("m1".to_string(), 1.0)]),
            None,
        )
        .unwrap();
        let numerator = "edges[7][0]**2";
        let value = evaluate_expression(&parsed, &expression, numerator, &input)
            .unwrap()
            .value;

        assert!(
            (value - -17.21131191349002).abs() < 1.0e-10,
            "initial-cut repeated-spectator finite-pole identity mismatch: {value}"
        );
    }

    #[test]
    fn cut_aware_vector_matroid_factorization_preserves_cut_surfaces() {
        let edge = |edge_id: usize,
                    tail: usize,
                    head: usize,
                    loop_signature: [i32; 2],
                    external_signature: i32,
                    mass_key: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass_key.to_string()),
                signature: MomentumSignature {
                    loop_signature: loop_signature.to_vec(),
                    external_signature: vec![external_signature],
                },
                had_pow: false,
            }
        };
        // The two bubbles share vertex 0, so the denominator graph is
        // connected, while their loop-energy rows form two independent vector-
        // matroid components. The initial-state cut is internal to the first
        // bubble and supplies its external energy shift.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, [1, 0], 0, "m0"),
                edge(1, 1, 0, [1, 0], 1, "m1"),
                edge(2, 0, 2, [0, 1], 0, "m2"),
                edge(3, 2, 0, [0, 1], 0, "m3"),
                edge(4, 0, 1, [0, 0], 1, "m_cut"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 4,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..3).map(|node| (format!("n{node}"), node)).collect(),
        };
        let cases = [
            generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(Vec::new()),
                    ..Default::default()
                },
            )
            .unwrap(),
            generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(0, 2)]),
                    ..Default::default()
                },
            )
            .unwrap(),
        ];
        assert!(
            cases[1]
                .expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| variant
                    .origin
                    .as_deref()
                    .is_some_and(|origin| origin.contains("contact")))
        );
        let cut_component = ParsedGraph {
            internal_edges: [0, 1, 4]
                .into_iter()
                .enumerate()
                .map(|(local_id, full_id)| {
                    let mut edge = parsed.internal_edges[full_id].clone();
                    edge.edge_id = local_id;
                    edge.signature.loop_signature.truncate(1);
                    edge
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 2,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["k0".to_string()],
            external_names: parsed.external_names.clone(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let component = generate_3d_expression_from_parsed(
            &cut_component,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
                ..Default::default()
            },
        )
        .unwrap();
        let expected_cut_surfaces = component
            .surfaces
            .linear_surface_cache
            .iter()
            .filter(|surface| !surface.expression.external_terms.is_empty())
            .map(|surface| surface.expression.clone())
            .collect::<Vec<_>>();
        let factorized_cut_surfaces = cases[0]
            .expression
            .surfaces
            .linear_surface_cache
            .iter()
            .filter(|surface| !surface.expression.external_terms.is_empty())
            .map(|surface| surface.expression.clone())
            .collect::<Vec<_>>();
        assert!(!expected_cut_surfaces.is_empty());
        assert_eq!(factorized_cut_surfaces.len(), expected_cut_surfaces.len());
        assert!(
            expected_cut_surfaces
                .iter()
                .all(|surface| factorized_cut_surfaces.contains(surface))
        );
    }

    #[test]
    fn cut_aware_attached_uv_bubble_matches_disconnected_cff_product() {
        let edge = |edge_id: usize,
                    tail: usize,
                    head: usize,
                    loop_signature: [i32; 3],
                    external_signature: i32,
                    mass_key: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass_key.to_string()),
                signature: MomentumSignature {
                    loop_signature: loop_signature.to_vec(),
                    external_signature: vec![external_signature],
                },
                had_pow: false,
            }
        };
        // This is the GL24 quartic cograph with a repeated UV bubble attached
        // at vertex 0. Incidence therefore says one connected denominator
        // graph, while the rational loop-energy rows split into the physical
        // cograph and an independent UV vector-matroid component.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 3, [1, 0, 0], 0, "m0"),
                edge(1, 4, 0, [1, 0, 0], 1, "m0"),
                edge(2, 1, 2, [0, 1, 0], 0, "m0"),
                edge(3, 3, 1, [0, 1, 0], -1, "m0"),
                edge(4, 3, 4, [1, -1, 0], 1, "m0"),
                edge(5, 2, 4, [0, 1, 0], 0, "m1"),
                edge(6, 0, 1, [0, 0, 0], 1, "m_cut"),
                edge(7, 0, 5, [0, 0, 1], 0, "m_uv"),
                edge(8, 5, 0, [0, 0, 1], 0, "m_uv"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: vec![ParsedGraphInitialStateCutEdge {
                edge_id: 6,
                external_id: 0,
                external_sign: 1,
            }],
            loop_names: vec!["q1".to_string(), "q3".to_string(), "q_uv".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..6).map(|node| (format!("n{node}"), node)).collect(),
        };
        let generated = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 4)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::None,
                ..Default::default()
            },
        )
        .unwrap();
        // The lower contact contracts only the active cograph edge. The UV
        // bubble remains on its two physical occurrences and own OSE maps;
        // this structural support is independent of internal origin labels.
        let lower_contact_orientations = generated
            .expression
            .orientations
            .iter()
            .filter(|orientation| {
                orientation.variants.iter().any(|variant| {
                    variant.denominator_edges.len() == 7
                        && !variant
                            .denominator_edges
                            .contains(&linnet::half_edge::involution::EdgeIndex(0))
                        && variant
                            .denominator_edges
                            .contains(&linnet::half_edge::involution::EdgeIndex(7))
                        && variant
                            .denominator_edges
                            .contains(&linnet::half_edge::involution::EdgeIndex(8))
                })
            })
            .collect::<Vec<_>>();
        assert!(!lower_contact_orientations.is_empty());
        assert!(lower_contact_orientations.iter().all(|orientation| {
            (orientation.edge_energy_map[7]
                == LinearEnergyExpr::ose(linnet::half_edge::involution::EdgeIndex(7), -1)
                && orientation.edge_energy_map[8]
                    == LinearEnergyExpr::ose(linnet::half_edge::involution::EdgeIndex(8), 1))
                || (orientation.edge_energy_map[7]
                    == LinearEnergyExpr::ose(linnet::half_edge::involution::EdgeIndex(7), 1)
                    && orientation.edge_energy_map[8]
                        == LinearEnergyExpr::ose(linnet::half_edge::involution::EdgeIndex(8), -1))
        }));
        let quadratic_generated = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::None,
                ..Default::default()
            },
        )
        .unwrap();
        let double_quadratic_generated = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2), (2, 2)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::None,
                ..Default::default()
            },
        )
        .unwrap();
        let mut disconnected = parsed.clone();
        disconnected.internal_edges[7].tail = 5;
        disconnected.internal_edges[7].head = 6;
        disconnected.internal_edges[8].tail = 6;
        disconnected.internal_edges[8].head = 5;
        disconnected.node_name_to_internal =
            (0..7).map(|node| (format!("n{node}"), node)).collect();
        let disconnected_generated = generate_3d_expression(
            &disconnected,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 4)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::None,
                ..Default::default()
            },
        )
        .unwrap();
        let disconnected_quadratic_generated = generate_3d_expression(
            &disconnected,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::None,
                ..Default::default()
            },
        )
        .unwrap();
        let cograph = ParsedGraph {
            internal_edges: parsed.internal_edges[..7]
                .iter()
                .cloned()
                .map(|mut edge| {
                    edge.signature.loop_signature.truncate(2);
                    edge
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: parsed.initial_state_cut_edges.clone(),
            loop_names: parsed.loop_names[..2].to_vec(),
            external_names: parsed.external_names.clone(),
            node_name_to_internal: (0..5).map(|node| (format!("n{node}"), node)).collect(),
        };
        let cograph_double_quadratic_generated = generate_3d_expression(
            &cograph,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2), (2, 2)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::None,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            generated
                .expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .any(|variant| !variant
                    .denominator_edges
                    .contains(&linnet::half_edge::involution::EdgeIndex(0)))
        );
        for seed in [17, 1337, 9100] {
            let input = EvaluationInput::deterministic(
                &parsed,
                seed,
                &BTreeMap::from([
                    ("m0".to_string(), 0.0),
                    ("m1".to_string(), 0.1),
                    ("m_cut".to_string(), 1.0),
                    ("m_uv".to_string(), 20.0),
                ]),
                None,
            )
            .unwrap();
            let numerator = "edges[0][0]**4";
            let cograph_input = EvaluationInput {
                external_momenta: input.external_momenta.clone(),
                loop_spatial_momenta: input.loop_spatial_momenta[..2].to_vec(),
                masses: input.masses[..7].to_vec(),
                uniform_scale: None,
            };
            let uv_energy = (input.loop_spatial_momenta[2]
                .iter()
                .map(|component| component * component)
                .sum::<f64>()
                + input.masses[7].powi(2))
            .sqrt();
            let repeated_bubble_factor = 1.0 / (4.0 * uv_energy.powi(3));
            let disconnected_value = evaluate_expression(
                &disconnected,
                &disconnected_generated.expression,
                numerator,
                &input,
            )
            .unwrap()
            .value
                * disconnected_generated.core_global_prefactor_sign.factor() as f64;
            let generalized =
                evaluate_expression(&parsed, &generated.expression, numerator, &input)
                    .unwrap()
                    .value
                    * generated.core_global_prefactor_sign.factor() as f64;
            let scale = generalized
                .abs()
                .max(disconnected_value.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                (generalized - disconnected_value).abs() <= 1.0e-10 * scale,
                "cut-aware attached-UV-bubble quartic product differs from its disconnected CFF product at seed {seed}: connected={generalized:e}, disconnected={disconnected_value:e}"
            );
            let quadratic_numerator = "edges[0][0]**2";
            let quadratic = evaluate_expression(
                &parsed,
                &quadratic_generated.expression,
                quadratic_numerator,
                &input,
            )
            .unwrap()
            .value
                * quadratic_generated.core_global_prefactor_sign.factor() as f64;
            let disconnected_quadratic = evaluate_expression(
                &disconnected,
                &disconnected_quadratic_generated.expression,
                quadratic_numerator,
                &input,
            )
            .unwrap()
            .value
                * disconnected_quadratic_generated
                    .core_global_prefactor_sign
                    .factor() as f64;
            let quadratic_scale = quadratic
                .abs()
                .max(disconnected_quadratic.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                (quadratic - disconnected_quadratic).abs() <= 1.0e-10 * quadratic_scale,
                "cut-aware attached-UV-bubble quadratic contact differs from its disconnected CFF product at seed {seed}: connected={quadratic:e}, disconnected={disconnected_quadratic:e}"
            );

            let double_quadratic_numerator = "edges[0][0]**2*edges[2][0]**2";
            let double_quadratic = evaluate_expression(
                &parsed,
                &double_quadratic_generated.expression,
                double_quadratic_numerator,
                &input,
            )
            .unwrap()
            .value
                * double_quadratic_generated
                    .core_global_prefactor_sign
                    .factor() as f64;
            let cograph_double_quadratic = evaluate_expression(
                &cograph,
                &cograph_double_quadratic_generated.expression,
                double_quadratic_numerator,
                &cograph_input,
            )
            .unwrap()
            .value
                * cograph_double_quadratic_generated
                    .core_global_prefactor_sign
                    .factor() as f64
                * repeated_bubble_factor;
            let double_quadratic_scale = double_quadratic
                .abs()
                .max(cograph_double_quadratic.abs())
                .max(f64::MIN_POSITIVE);
            assert!(
                (double_quadratic - cograph_double_quadratic).abs()
                    <= 1.0e-10 * double_quadratic_scale,
                "cut-aware attached-UV-bubble double-quadratic recursion differs from its cograph CFF times the repeated bubble at seed {seed}: generated={double_quadratic:e}, cograph generated times bubble={cograph_double_quadratic:e}, full core sign={} cograph core sign={}",
                double_quadratic_generated
                    .core_global_prefactor_sign
                    .factor(),
                cograph_double_quadratic_generated
                    .core_global_prefactor_sign
                    .factor(),
            );
        }
    }

    #[test]
    fn higher_power_cff_is_invariant_under_nonzero_uniform_scale() {
        let parsed = crate::graph_io::test_graphs::box_pow3_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(3, 4)]),
                numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                ..Default::default()
            },
        )
        .unwrap();
        let mut input =
            EvaluationInput::deterministic(&parsed, 1337, &BTreeMap::new(), Some(0.75)).unwrap();
        let first = evaluate_expression(&parsed, &expression, "edges[3][0]**4", &input)
            .unwrap()
            .value;
        input.uniform_scale = Some(2.25);
        let second = evaluate_expression(&parsed, &expression, "edges[3][0]**4", &input)
            .unwrap()
            .value;
        let scale = first.abs().max(second.abs()).max(1.0e-14);

        assert!(
            (first - second).abs() <= 1.0e-9 * scale,
            "higher-power CFF changed with auxiliary M: {first:e} vs {second:e}"
        );
    }

    #[test]
    fn repeated_cubic_pole_with_squared_denominator_numerator_pinches_to_simple_pole() {
        let repeated = crate::graph_io::test_graphs::box_pow3_graph();
        let simple = crate::graph_io::test_graphs::box_graph();
        let repeated_expression = generate_3d_expression_from_parsed(
            &repeated,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(3, 4)]),
                ..Default::default()
            },
        )
        .unwrap();
        let simple_expression = generate_3d_expression_from_parsed(
            &simple,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
                ..Default::default()
            },
        )
        .unwrap();
        let repeated_input = EvaluationInput::deterministic(
            &repeated,
            1337,
            &BTreeMap::from([
                ("m1".to_string(), 0.41),
                ("m2".to_string(), 0.53),
                ("m3".to_string(), 0.67),
                ("m4".to_string(), 0.7),
            ]),
            None,
        )
        .unwrap();
        let simple_input = EvaluationInput {
            masses: repeated_input.masses[..4].to_vec(),
            ..repeated_input.clone()
        };
        let repeated_value = evaluate_expression(
            &repeated,
            &repeated_expression,
            "(dot(edges[3],edges[3])-0.49)**2",
            &repeated_input,
        )
        .unwrap()
        .value;
        let simple_value = evaluate_expression(&simple, &simple_expression, "1", &simple_input)
            .unwrap()
            .value;
        let scale = repeated_value
            .abs()
            .max(simple_value.abs())
            .max(f64::MIN_POSITIVE);

        assert!(
            (repeated_value - simple_value).abs() <= 1.0e-11 * scale,
            "D(q)^2 / D(q)^3 did not pinch to 1 / D(q): repeated={repeated_value:e}, simple={simple_value:e}"
        );
    }

    #[test]
    fn repeated_quintic_channel_denominator_factors_lower_to_quartic_and_cubic_channels() {
        let repeated_channel = |power: usize, reverse_last: bool| ParsedGraph {
            internal_edges: (0..power)
                .map(|edge_id| {
                    let mut tail = edge_id;
                    let mut head = (edge_id + 1) % power;
                    let mut loop_sign = 1;
                    if reverse_last && edge_id + 1 == power {
                        std::mem::swap(&mut tail, &mut head);
                        loop_sign = -1;
                    }
                    ParsedGraphInternalEdge {
                        edge_id,
                        tail,
                        head,
                        label: format!("q{edge_id}"),
                        mass_key: Some("m".to_string()),
                        signature: MomentumSignature {
                            loop_signature: vec![loop_sign],
                            external_signature: Vec::new(),
                        },
                        had_pow: false,
                    }
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q0".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power).map(|node| (format!("v{node}"), node)).collect(),
        };
        let generate = |parsed: &ParsedGraph, degree, numerator_sampling_scale| {
            generate_3d_expression(
                parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(0, degree)]),
                    numerator_sampling_scale,
                    ..Default::default()
                },
            )
            .unwrap()
        };

        let mut diagnostics = Vec::new();
        let mut failures = Vec::new();
        for (sampling_mode, uniform_scale) in [
            (NumeratorSamplingScaleMode::None, None),
            (NumeratorSamplingScaleMode::BeyondQuadratic, Some(0.91)),
            (NumeratorSamplingScaleMode::BeyondQuadratic, Some(1.37)),
        ] {
            for reverse_last in [false, true] {
                let graphs = [5, 4, 3].map(|power| repeated_channel(power, reverse_last));
                let generated = [
                    generate(&graphs[0], 6, sampling_mode),
                    generate(&graphs[1], 4, sampling_mode),
                    generate(&graphs[2], 2, sampling_mode),
                ];
                let distributed = [
                    generate_3d_expression(
                        &graphs[0],
                        &Generate3DExpressionOptions {
                            energy_degree_bounds: Some(vec![
                                (0, 2),
                                (1, 1),
                                (2, 1),
                                (3, 1),
                                (4, 1),
                            ]),
                            numerator_sampling_scale: sampling_mode,
                            ..Default::default()
                        },
                    )
                    .unwrap(),
                    generate_3d_expression(
                        &graphs[1],
                        &Generate3DExpressionOptions {
                            energy_degree_bounds: Some(vec![(0, 1), (1, 1), (2, 1), (3, 1)]),
                            numerator_sampling_scale: sampling_mode,
                            ..Default::default()
                        },
                    )
                    .unwrap(),
                ];
                let mut contact = generated[0].expression.clone();
                let mut remainder = generated[0].expression.clone();
                for orientation in &mut contact.orientations {
                    orientation.variants.retain(|variant| {
                        variant
                            .origin
                            .as_deref()
                            .is_some_and(|origin| origin.contains("contact"))
                    });
                }
                for orientation in &mut remainder.orientations {
                    orientation.variants.retain(|variant| {
                        !variant
                            .origin
                            .as_deref()
                            .is_some_and(|origin| origin.contains("contact"))
                    });
                }
                contact
                    .orientations
                    .retain(|orientation| !orientation.variants.is_empty());
                remainder
                    .orientations
                    .retain(|orientation| !orientation.variants.is_empty());

                for (point, loop_spatial_momentum) in [
                    [0.31, -0.47, 0.83],
                    [-0.59, 0.23, 1.07],
                    [1.19, -0.71, 0.43],
                ]
                .into_iter()
                .enumerate()
                {
                    let inputs = [5, 4, 3].map(|power| EvaluationInput {
                        external_momenta: Vec::new(),
                        loop_spatial_momenta: vec![loop_spatial_momentum],
                        masses: vec![0.73; power],
                        uniform_scale,
                    });
                    let denominator = "dot(edges[0],edges[0])-0.5329";
                    for (label, quintic_numerator, lower_numerator, lower_index) in [
                        (
                            "D(q)*q0^4 / D(q)^5 = q0^4 / D(q)^4",
                            format!("({denominator})*edges[0][0]**4"),
                            "edges[0][0]**4",
                            1,
                        ),
                        (
                            "D(q)*q0^3 / D(q)^5 = q0^3 / D(q)^4",
                            format!("({denominator})*edges[0][0]**3"),
                            "edges[0][0]**3",
                            1,
                        ),
                        (
                            "D(q)^2*q0^2 / D(q)^5 = q0^2 / D(q)^3",
                            format!("({denominator})**2*edges[0][0]**2"),
                            "edges[0][0]**2",
                            2,
                        ),
                    ] {
                        let raw_quintic = evaluate_expression(
                            &graphs[0],
                            &generated[0].expression,
                            &quintic_numerator,
                            &inputs[0],
                        )
                        .unwrap()
                        .value;
                        let raw_lower = evaluate_expression(
                            &graphs[lower_index],
                            &generated[lower_index].expression,
                            lower_numerator,
                            &inputs[lower_index],
                        )
                        .unwrap()
                        .value;
                        let source_frame = |index: usize| {
                            generated[index].core_global_prefactor_sign.factor() as f64
                                * if graphs[index].internal_edges.len().is_multiple_of(2) {
                                    1.0
                                } else {
                                    -1.0
                                }
                        };
                        let quintic = raw_quintic;
                        let lower = raw_lower;
                        let contact_value = evaluate_expression(
                            &graphs[0],
                            &contact,
                            &quintic_numerator,
                            &inputs[0],
                        )
                        .unwrap()
                        .value;
                        let remainder_value = evaluate_expression(
                            &graphs[0],
                            &remainder,
                            &quintic_numerator,
                            &inputs[0],
                        )
                        .unwrap()
                        .value;
                        let denominator_edge_bins = generated[0]
                            .expression
                            .orientations
                            .iter()
                            .flat_map(|orientation| &orientation.variants)
                            .map(|variant| variant.denominator_edges.len())
                            .collect::<std::collections::BTreeSet<_>>()
                            .into_iter()
                            .map(|denominator_edges| {
                                let mut sector = generated[0].expression.clone();
                                for orientation in &mut sector.orientations {
                                    orientation.variants.retain(|variant| {
                                        variant.denominator_edges.len() == denominator_edges
                                    });
                                }
                                sector
                                    .orientations
                                    .retain(|orientation| !orientation.variants.is_empty());
                                let value = evaluate_expression(
                                    &graphs[0],
                                    &sector,
                                    &quintic_numerator,
                                    &inputs[0],
                                )
                                .unwrap()
                                .value;
                                (denominator_edges, value)
                            })
                            .collect::<Vec<_>>();
                        let lower_denominator_edge_bins = generated[lower_index]
                            .expression
                            .orientations
                            .iter()
                            .flat_map(|orientation| &orientation.variants)
                            .map(|variant| variant.denominator_edges.len())
                            .collect::<std::collections::BTreeSet<_>>()
                            .into_iter()
                            .map(|denominator_edges| {
                                let mut sector = generated[lower_index].expression.clone();
                                for orientation in &mut sector.orientations {
                                    orientation.variants.retain(|variant| {
                                        variant.denominator_edges.len() == denominator_edges
                                    });
                                }
                                sector
                                    .orientations
                                    .retain(|orientation| !orientation.variants.is_empty());
                                let value = evaluate_expression(
                                    &graphs[lower_index],
                                    &sector,
                                    lower_numerator,
                                    &inputs[lower_index],
                                )
                                .unwrap()
                                .value;
                                (denominator_edges, value)
                            })
                            .collect::<Vec<_>>();
                        let scale = quintic.abs().max(lower.abs()).max(f64::MIN_POSITIVE);
                        let diagnostic = format!(
                            "{label} with {sampling_mode:?}, M={uniform_scale:?}, point {point}, reverse_last={reverse_last}: quintic={quintic:e}, lower={lower:e}, GammaLoop source frames=({}, {}), contact={contact_value:e}, remainder={remainder_value:e}, quintic bins={denominator_edge_bins:?}, lower bins={lower_denominator_edge_bins:?}",
                            source_frame(0),
                            source_frame(lower_index),
                        );
                        if (quintic - lower).abs() > 1.0e-10 * scale.max(1.0) {
                            failures.push(diagnostic.clone());
                        }
                        diagnostics.push(diagnostic);
                    }

                    // Mirror the exact-source minimax plan: the retained
                    // quartic uses four distinct occurrences, while the two
                    // temporal factors of D(q) use the fifth occurrence and
                    // one already carrying a retained factor.
                    let physical_energy = |occurrence: usize, power: usize| {
                        if reverse_last && occurrence + 1 == power {
                            format!("(-edges[{occurrence}][0])")
                        } else {
                            format!("edges[{occurrence}][0]")
                        }
                    };
                    let retained = |power: usize| {
                        (0..4)
                            .map(|occurrence| {
                                format!("({}+0.37)", physical_energy(occurrence, power))
                            })
                            .collect::<Vec<_>>()
                            .join("*")
                    };
                    let common_numerator = format!(
                        "({}*{}+dot(edges[0],edges[0])-edges[0][0]**2-0.5329)*{}",
                        physical_energy(4, 5),
                        physical_energy(0, 5),
                        retained(5),
                    );
                    let lower_numerator = retained(4);
                    let common = evaluate_expression(
                        &graphs[0],
                        &distributed[0].expression,
                        &common_numerator,
                        &inputs[0],
                    )
                    .unwrap()
                    .value;
                    let lower = evaluate_expression(
                        &graphs[1],
                        &distributed[1].expression,
                        &lower_numerator,
                        &inputs[1],
                    )
                    .unwrap()
                    .value;
                    let scale = common.abs().max(lower.abs()).max(f64::MIN_POSITIVE);
                    let diagnostic = format!(
                        "distributed exact plan with {sampling_mode:?}, M={uniform_scale:?}, point {point}, reverse_last={reverse_last}: common={common:e}, lower={lower:e}, common numerator={common_numerator}, lower numerator={lower_numerator}"
                    );
                    if (common - lower).abs() > 1.0e-10 * scale.max(1.0) {
                        failures.push(diagnostic.clone());
                    }
                    diagnostics.push(diagnostic);
                }
            }
        }
        assert!(
            failures.is_empty(),
            "failures:\n{}\nall diagnostics:\n{}",
            failures.join("\n"),
            diagnostics.join("\n")
        );
    }

    #[test]
    fn repeated_quintic_rank_eight_denominator_pinches_match_lower_channels() {
        let repeated_channel = |power: usize| ParsedGraph {
            internal_edges: (0..power)
                .map(|edge_id| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: (edge_id + 1) % power,
                    label: format!("q{edge_id}"),
                    mass_key: Some("m".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q0".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power).map(|node| (format!("v{node}"), node)).collect(),
        };
        let graphs = [5, 4, 3, 2].map(repeated_channel);

        for (sampling_mode, uniform_scale) in [
            (NumeratorSamplingScaleMode::None, None),
            (NumeratorSamplingScaleMode::BeyondQuadratic, Some(0.91)),
            (NumeratorSamplingScaleMode::BeyondQuadratic, Some(1.37)),
        ] {
            let common = generate_3d_expression(
                &graphs[0],
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(0, 8)]),
                    numerator_sampling_scale: sampling_mode,
                    ..Default::default()
                },
            )
            .unwrap();
            let lowers = [6, 4, 2].map(|degree| {
                generate_3d_expression(
                    &graphs[4 - degree / 2],
                    &Generate3DExpressionOptions {
                        energy_degree_bounds: Some(vec![(0, degree)]),
                        numerator_sampling_scale: sampling_mode,
                        ..Default::default()
                    },
                )
                .unwrap()
            });

            for (point, loop_spatial_momentum) in [
                [0.31, -0.47, 0.83],
                [-0.59, 0.23, 1.07],
                [1.19, -0.71, 0.43],
            ]
            .into_iter()
            .enumerate()
            {
                let common_input = EvaluationInput {
                    external_momenta: Vec::new(),
                    loop_spatial_momenta: vec![loop_spatial_momentum],
                    masses: vec![0.73; 5],
                    uniform_scale,
                };
                let denominator = "dot(edges[0],edges[0])-0.5329";

                for (pinched_factors, lower) in (1..=3).zip(&lowers) {
                    let energy_degree = 8 - 2 * pinched_factors;
                    let common_numerator =
                        format!("({denominator})**{pinched_factors}*edges[0][0]**{energy_degree}");
                    let lower_numerator = format!("edges[0][0]**{energy_degree}");
                    let common_value = evaluate_expression(
                        &graphs[0],
                        &common.expression,
                        &common_numerator,
                        &common_input,
                    )
                    .unwrap()
                    .value;
                    let lower_input = EvaluationInput {
                        masses: vec![0.73; 5 - pinched_factors],
                        ..common_input.clone()
                    };
                    let lower_value = evaluate_expression(
                        &graphs[pinched_factors],
                        &lower.expression,
                        &lower_numerator,
                        &lower_input,
                    )
                    .unwrap()
                    .value;
                    let scale = common_value
                        .abs()
                        .max(lower_value.abs())
                        .max(f64::MIN_POSITIVE);
                    assert!(
                        (common_value - lower_value).abs() <= 1.0e-10 * scale,
                        "rank-eight quintic lowering failed with {sampling_mode:?}, M={uniform_scale:?}, point={point}, pinches={pinched_factors}: common={common_value:e}, lower={lower_value:e}"
                    );
                }
            }
        }
    }

    #[test]
    fn repeated_sextic_cycle_exact_mapper_dispatch_lowers_to_quintic_cycle() {
        let repeated_cycle = |power: usize| {
            let (endpoints, node_name_to_internal) = match power {
                6 => (
                    vec![(0, 4), (1, 5), (2, 3), (3, 1), (4, 2), (5, 0)],
                    [
                        ("__gammaloop_exact_power_0_0", 5),
                        ("__gammaloop_exact_power_0_1", 1),
                        ("__gammaloop_exact_power_0_2", 3),
                        ("__gammaloop_exact_power_0_3", 2),
                        ("n0", 0),
                        ("n1", 4),
                    ]
                    .into_iter()
                    .map(|(name, node)| (name.to_string(), node))
                    .collect(),
                ),
                5 => (
                    vec![(0, 3), (1, 4), (2, 1), (3, 2), (4, 0)],
                    [
                        ("__gammaloop_exact_power_0_0", 4),
                        ("__gammaloop_exact_power_0_1", 1),
                        ("__gammaloop_exact_power_0_2", 2),
                        ("n0", 0),
                        ("n1", 3),
                    ]
                    .into_iter()
                    .map(|(name, node)| (name.to_string(), node))
                    .collect(),
                ),
                _ => unreachable!("the exact-source oracle only uses P6 and P5"),
            };
            ParsedGraph {
                internal_edges: endpoints
                    .into_iter()
                    .enumerate()
                    .map(|(edge_id, (tail, head))| ParsedGraphInternalEdge {
                        edge_id,
                        tail,
                        head,
                        label: format!("__gammaloop_exact_edge_{edge_id}"),
                        mass_key: Some("1".to_string()),
                        signature: MomentumSignature {
                            loop_signature: vec![-1],
                            external_signature: Vec::new(),
                        },
                        had_pow: false,
                    })
                    .collect(),
                external_edges: Vec::new(),
                initial_state_cut_edges: Vec::new(),
                loop_names: vec!["ell0".to_string()],
                external_names: Vec::new(),
                node_name_to_internal,
            }
        };
        let sextic = repeated_cycle(6);
        let quintic = repeated_cycle(5);
        let mut failures = Vec::new();
        for (sampling_mode, uniform_scale) in [
            (NumeratorSamplingScaleMode::None, None),
            (NumeratorSamplingScaleMode::BeyondQuadratic, Some(0.91)),
            (NumeratorSamplingScaleMode::BeyondQuadratic, Some(1.37)),
        ] {
            let split = generate_3d_expression(
                &sextic,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(1, 2), (2, 1), (3, 1), (4, 1), (5, 1)]),
                    numerator_sampling_scale: sampling_mode,
                    ..Default::default()
                },
            )
            .unwrap();
            let coherent = generate_3d_expression(
                &sextic,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(1, 1), (2, 1), (3, 1), (4, 1), (5, 2)]),
                    numerator_sampling_scale: sampling_mode,
                    ..Default::default()
                },
            )
            .unwrap();
            let lower = generate_3d_expression(
                &quintic,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(1, 1), (2, 1), (3, 1), (4, 1)]),
                    numerator_sampling_scale: sampling_mode,
                    ..Default::default()
                },
            )
            .unwrap();
            let retained = (1..=4)
                .map(|occurrence| format!("(-edges[{occurrence}][0]+0.37)"))
                .collect::<Vec<_>>()
                .join("*");
            let spatial_mass = "dot(edges[5],edges[5])-edges[5][0]**2-1.0";
            let split_numerator =
                format!("((-edges[5][0])*(-edges[1][0])+{spatial_mass})*{retained}");
            let coherent_numerator = format!("((-edges[5][0])**2+{spatial_mass})*{retained}");

            for (point, loop_spatial_momentum) in [
                [0.31, -0.47, 0.83],
                [-0.59, 0.23, 1.07],
                [1.19, -0.71, 0.43],
            ]
            .into_iter()
            .enumerate()
            {
                let sextic_input = EvaluationInput {
                    external_momenta: Vec::new(),
                    loop_spatial_momenta: vec![loop_spatial_momentum],
                    masses: vec![1.0; 6],
                    uniform_scale,
                };
                let quintic_input = EvaluationInput {
                    masses: vec![1.0; 5],
                    ..sextic_input.clone()
                };
                let split_value = evaluate_expression(
                    &sextic,
                    &split.expression,
                    &split_numerator,
                    &sextic_input,
                )
                .unwrap()
                .value;
                let coherent_value = evaluate_expression(
                    &sextic,
                    &coherent.expression,
                    &coherent_numerator,
                    &sextic_input,
                )
                .unwrap()
                .value;
                let lower_value =
                    evaluate_expression(&quintic, &lower.expression, &retained, &quintic_input)
                        .unwrap()
                        .value;
                for (label, candidate) in [
                    ("split denominator", split_value),
                    ("coherent denominator", coherent_value),
                ] {
                    let scale = candidate
                        .abs()
                        .max(lower_value.abs())
                        .max(f64::MIN_POSITIVE);
                    if (candidate - lower_value).abs() > 1.0e-10 * scale.max(1.0) {
                        failures.push(format!(
                            "{label}, {sampling_mode:?}, M={uniform_scale:?}, point={point}: common={candidate:e}, lower={lower_value:e}, common numerator={split_numerator}, coherent numerator={coherent_numerator}"
                        ));
                    }
                }
            }
        }
        assert!(failures.is_empty(), "failures:\n{}", failures.join("\n"));
    }

    #[test]
    fn repeated_quintic_channel_lowers_inside_a_mixed_cograph_product() {
        let mixed_channel = |power: usize| ParsedGraph {
            internal_edges: [
                ParsedGraphInternalEdge {
                    edge_id: 0,
                    tail: 0,
                    head: 1,
                    label: "c0".to_string(),
                    mass_key: Some("mc0".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![1, 0],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
                ParsedGraphInternalEdge {
                    edge_id: 1,
                    tail: 1,
                    head: 0,
                    label: "c1".to_string(),
                    mass_key: Some("mc1".to_string()),
                    signature: MomentumSignature {
                        loop_signature: vec![-1, 0],
                        external_signature: Vec::new(),
                    },
                    had_pow: false,
                },
            ]
            .into_iter()
            .chain((0..power).map(|occurrence| ParsedGraphInternalEdge {
                edge_id: occurrence + 2,
                tail: occurrence + 2,
                head: (occurrence + 1) % power + 2,
                label: format!("u{occurrence}"),
                mass_key: Some("muv".to_string()),
                signature: MomentumSignature {
                    loop_signature: vec![0, 1],
                    external_signature: Vec::new(),
                },
                had_pow: false,
            }))
            .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["kc".to_string(), "kuv".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..power + 2)
                .map(|node| (format!("v{node}"), node))
                .collect(),
        };
        let graphs = [mixed_channel(5), mixed_channel(4)];
        let generated = [
            generate_3d_expression(
                &graphs[0],
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(0, 1), (2, 6)]),
                    numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                    ..Default::default()
                },
            )
            .unwrap(),
            generate_3d_expression(
                &graphs[1],
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(0, 1), (2, 4)]),
                    numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                    ..Default::default()
                },
            )
            .unwrap(),
        ];
        let origins = generated
            .iter()
            .map(|generated| {
                generated
                    .expression
                    .orientations
                    .iter()
                    .flat_map(|orientation| &orientation.variants)
                    .filter_map(|variant| variant.origin.clone())
                    .collect::<std::collections::BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        let mut failures = Vec::new();
        for (point, (cograph_momentum, uv_momentum)) in [
            ([0.31, -0.47, 0.83], [-0.59, 0.23, 1.07]),
            ([-0.41, 0.73, 0.29], [1.19, -0.71, 0.43]),
            ([0.67, 0.37, -0.53], [-0.89, 1.13, 0.17]),
        ]
        .into_iter()
        .enumerate()
        {
            let inputs = [5, 4].map(|power| EvaluationInput {
                external_momenta: Vec::new(),
                loop_spatial_momenta: vec![cograph_momentum, uv_momentum],
                masses: [0.41, 0.67]
                    .into_iter()
                    .chain(std::iter::repeat_n(0.73, power))
                    .collect(),
                uniform_scale: Some(0.91),
            });
            let common = evaluate_expression(
                &graphs[0],
                &generated[0].expression,
                "(dot(edges[2],edges[2])-0.5329)*edges[2][0]**4",
                &inputs[0],
            )
            .unwrap()
            .value;
            let lower = evaluate_expression(
                &graphs[1],
                &generated[1].expression,
                "edges[2][0]**4",
                &inputs[1],
            )
            .unwrap()
            .value;
            let scale = common.abs().max(lower.abs()).max(f64::MIN_POSITIVE);
            if (common - lower).abs() > 1.0e-10 * scale {
                failures.push(format!(
                    "factorized UV cancellation, point {point}: common={common:e}, lower={lower:e}, ratio={}, common core={}, lower core={}, origins={origins:?}",
                    common / lower,
                    generated[0].core_global_prefactor_sign.factor(),
                    generated[1].core_global_prefactor_sign.factor(),
                ));
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }

    #[test]
    fn lower_sector_powered_pole_contact_reconstructs_numerator_derivatives() {
        let edge = |edge_id: usize,
                    tail: usize,
                    head: usize,
                    label: &str,
                    loop_signature: Vec<i32>,
                    mass_key: &str| {
            crate::graph_io::ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: label.to_string(),
                mass_key: Some(mass_key.to_string()),
                signature: crate::MomentumSignature {
                    loop_signature,
                    external_signature: Vec::new(),
                },
                had_pow: false,
            }
        };
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 0, "z", vec![0, 0, 1], "mz"),
                edge(1, 0, 3, "r", vec![1, -1, 0], "mr"),
                edge(2, 0, 1, "q2", vec![1, 1, 0], "mq"),
                edge(3, 1, 2, "q3", vec![1, 1, 0], "mq"),
                edge(4, 2, 0, "q4", vec![1, 1, 0], "mq"),
                edge(5, 3, 4, "r1", vec![1, -1, 0], "mr"),
                edge(6, 4, 0, "r2", vec![1, -1, 0], "mr"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k0".to_string(), "k1".to_string(), "k2".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..5).map(|node| (format!("n{node}"), node)).collect(),
        };
        let contact_terms = |mut expression: ThreeDExpression<OrientationID>| {
            for orientation in &mut expression.orientations {
                orientation.variants.retain(|variant| {
                    variant.origin.as_deref().is_some_and(|origin| {
                        origin.starts_with("bounded_degree_quadratic_recursive_contact")
                    })
                });
            }
            expression
                .orientations
                .retain(|orientation| !orientation.variants.is_empty());
            expression
        };
        let contact = contact_terms(
            generate_3d_expression_from_parsed(
                &parsed,
                &Generate3DExpressionOptions {
                    // z^2 is carried by edge 0, while k0*k1 is bounded by the
                    // independent r and q edge-energy coordinates.
                    energy_degree_bounds: Some(vec![(0, 2), (1, 2), (2, 2)]),
                    ..Default::default()
                },
            )
            .unwrap(),
        );
        let scalar_contact = contact_terms(
            generate_3d_expression_from_parsed(
                &parsed,
                &Generate3DExpressionOptions {
                    energy_degree_bounds: Some(vec![(0, 2)]),
                    ..Default::default()
                },
            )
            .unwrap(),
        );

        assert!(!contact.orientations.is_empty());
        assert!(
            contact
                .orientations
                .iter()
                .all(|orientation| orientation.loop_energy_map.len() == 3),
            "lower-sector samples must be lifted back to the original loop basis"
        );
        let map_coeff = |expr: &LinearEnergyExpr, edge_id: usize| {
            expr.internal_terms
                .iter()
                .find_map(|(candidate, coefficient)| {
                    (candidate.0 == edge_id).then(|| rational_to_f64(coefficient))
                })
                .unwrap_or(0.0)
        };
        assert!(contact.orientations.iter().any(|orientation| {
            let [k0, k1, _] = orientation.loop_energy_map.as_slice() else {
                return false;
            };
            map_coeff(k0, 1) == 0.5
                && map_coeff(k0, 2) == 0.5
                && map_coeff(k1, 1) == -0.5
                && map_coeff(k1, 2) == 0.5
        }));

        let input = EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![[0.0; 3]; 3],
            masses: vec![0.73, 1.11, 0.89, 0.89, 0.89, 1.11, 1.11],
            uniform_scale: None,
        };
        let actual = evaluate_expression(
            &parsed,
            &contact,
            "edges[0][0]**2*(edges[2][0]**2-edges[1][0]**2)/4",
            &input,
        )
        .unwrap()
        .value;
        let actual_q_squared =
            evaluate_expression(&parsed, &contact, "edges[0][0]**2*edges[2][0]**2", &input)
                .unwrap()
                .value;
        let actual_r_squared =
            evaluate_expression(&parsed, &contact, "edges[0][0]**2*edges[1][0]**2", &input)
                .unwrap()
                .value;
        let parent_scalar = evaluate_expression(&parsed, &scalar_contact, "edges[0][0]**2", &input)
            .unwrap()
            .value;

        let residual = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 3, "r", vec![1, 0], "mr"),
                edge(1, 3, 4, "r1", vec![1, 0], "mr"),
                edge(2, 4, 0, "r2", vec![1, 0], "mr"),
                edge(3, 0, 1, "q1", vec![0, 1], "mq"),
                edge(4, 1, 2, "q2", vec![0, 1], "mq"),
                edge(5, 2, 0, "q3", vec![0, 1], "mq"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["r0".to_string(), "q0".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..5).map(|node| (format!("n{node}"), node)).collect(),
        };
        let residual_scalar_expression = generate_3d_expression_from_parsed(
            &residual,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(Vec::new()),
                ..Default::default()
            },
        )
        .unwrap();
        let residual_scalar = evaluate_expression(
            &residual,
            &residual_scalar_expression,
            "1",
            &EvaluationInput {
                external_momenta: Vec::new(),
                loop_spatial_momenta: vec![[0.0; 3]; 2],
                masses: vec![1.11, 1.11, 1.11, 0.89, 0.89, 0.89],
                uniform_scale: None,
            },
        )
        .unwrap()
        .value;

        // At z contact, k0*k1=(q^2-r^2)/4. Three ordinary denominators in
        // each of the r and q channels keep every energy subcycle convergent,
        // while both quadratic source factors retain their sampled value and
        // beta=2 numerator derivative.
        let er = 1.11_f64;
        let eq = 0.89_f64;
        let i2_r = -1.0 / (4.0 * er.powi(3));
        let i3_r = 3.0 / (16.0 * er.powi(5));
        let i2_q = -1.0 / (4.0 * eq.powi(3));
        let i3_q = 3.0 / (16.0 * eq.powi(5));
        let sampled_value = i3_r * (eq.powi(2) - er.powi(2)) * i3_q / 4.0;
        let numerator_derivative = (i3_r * i2_q - i2_r * i3_q) / 4.0;
        let expected = sampled_value + numerator_derivative;
        let expected_q_squared = i3_r * (i2_q + eq.powi(2) * i3_q);
        let expected_r_squared = (i2_r + er.powi(2) * i3_r) * i3_q;

        assert!(parent_scalar.is_sign_positive());
        assert!(residual_scalar.is_sign_negative());
        assert!((parent_scalar + residual_scalar).abs() < 1.0e-13);
        assert!((residual_scalar + i3_r * i3_q).abs() < 1.0e-13);
        assert!(sampled_value.is_sign_negative());
        assert!(numerator_derivative.is_sign_positive());
        assert!(
            (actual_q_squared - expected_q_squared).abs() < 1.0e-13,
            "q-squared contact {actual_q_squared} differs from analytic {expected_q_squared}"
        );
        assert!(
            (actual_r_squared - expected_r_squared).abs() < 1.0e-13,
            "r-squared contact {actual_r_squared} differs from analytic {expected_r_squared}"
        );
        assert!(
            (actual - expected).abs() < 1.0e-13,
            "rank-projected contact {actual} differs from analytic {expected}"
        );
    }

    #[test]
    fn rank_projected_constant_contact_preserves_even_duplicate_parity() {
        let edge = |edge_id: usize,
                    tail: usize,
                    head: usize,
                    label: &str,
                    loop_signature: Vec<i32>,
                    mass_key: &str| {
            crate::graph_io::ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: label.to_string(),
                mass_key: Some(mass_key.to_string()),
                signature: crate::MomentumSignature {
                    loop_signature,
                    external_signature: Vec::new(),
                },
                had_pow: false,
            }
        };
        let parent = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 0, "z", vec![0, 1], "mz"),
                edge(1, 0, 1, "q1", vec![1, 0], "mq"),
                edge(2, 1, 0, "q2", vec![1, 0], "mq"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["kq".to_string(), "kz".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let options = Generate3DExpressionOptions {
            energy_degree_bounds: Some(vec![(0, 2)]),
            ..Default::default()
        };
        let generated = generate_3d_expression(&parent, &options).unwrap();
        let bounded_component = generated
            .energy_factor_components
            .iter()
            .find(|component| component.internal_edge_ids.contains(&0))
            .expect("the rank-projected edge must retain a bounded energy-factor component");
        assert_eq!(
            bounded_component.ownership,
            crate::CffEnergyFactorOwnership::VariantLocal,
            "the rank-projected contact must keep its bounded energy factors variant-local"
        );
        let contact_origins = generated
            .expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .filter_map(|variant| variant.origin.as_deref())
            .filter(|origin| origin.contains("contact"))
            .collect::<Vec<_>>();
        for branch in ["e0=-", "e0=0", "e0=+"] {
            assert!(
                contact_origins
                    .iter()
                    .any(|origin| origin.ends_with(branch)),
                "the component-wrapped quadratic contact lacks its deterministic {branch} branch: {contact_origins:?}"
            );
        }
        let contact_terms = |mut expression: ThreeDExpression<OrientationID>| {
            for orientation in &mut expression.orientations {
                orientation.variants.retain(|variant| {
                    variant
                        .origin
                        .as_deref()
                        .is_some_and(|origin| origin.contains("contact"))
                });
            }
            expression
                .orientations
                .retain(|orientation| !orientation.variants.is_empty());
            expression
        };
        let projected = contact_terms(generated.expression);
        let parent_input = EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![[0.0; 3]; 2],
            masses: vec![0.73, 0.89, 0.89],
            uniform_scale: None,
        };
        let projected_value =
            evaluate_expression(&parent, &projected, "edges[0][0]**2", &parent_input)
                .unwrap()
                .value;
        let residual = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, "q1", vec![1], "mq"),
                edge(1, 1, 0, "q2", vec![1], "mq"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["kq".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: BTreeMap::from([("n0".to_string(), 0), ("n1".to_string(), 1)]),
        };
        let residual_expression =
            generate_3d_expression_from_parsed(&residual, &Generate3DExpressionOptions::default())
                .unwrap();
        let residual_value = evaluate_expression(
            &residual,
            &residual_expression,
            "1",
            &EvaluationInput {
                external_momenta: Vec::new(),
                loop_spatial_momenta: vec![[0.0; 3]],
                masses: vec![0.89, 0.89],
                uniform_scale: None,
            },
        )
        .unwrap()
        .value;
        assert!(residual_value.abs() > 1.0e-12);
        assert!((projected_value - residual_value).abs() < 1.0e-13);
    }

    #[test]
    fn gl24_rank_four_contact_preserves_denominator_cancellation() {
        let edge = |edge_id: usize,
                    tail: usize,
                    head: usize,
                    loop_signature: [i32; 2],
                    external_signature: i32,
                    mass_key: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass_key.to_string()),
                signature: MomentumSignature {
                    loop_signature: loop_signature.to_vec(),
                    external_signature: vec![external_signature],
                },
                had_pow: false,
            }
        };
        let parent = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 2, [-1, 1], -1, "m0"),
                edge(1, 0, 4, [0, -1], 0, "m1"),
                edge(2, 1, 0, [-1, 0], -1, "m0"),
                edge(3, 2, 1, [-1, 0], 0, "m0"),
                edge(4, 3, 2, [0, -1], 1, "m0"),
                edge(5, 4, 3, [0, -1], 0, "m0"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["k0".to_string(), "k1".to_string()],
            external_names: vec!["p0".to_string()],
            node_name_to_internal: (0..5).map(|node| (format!("v{node}"), node)).collect(),
        };
        let parent_generated = generate_3d_expression(
            &parent,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 4)]),
                ..Default::default()
            },
        )
        .unwrap();

        let residual = ParsedGraph {
            // Contracting edge 0 identifies parent nodes 0 and 2.
            internal_edges: vec![
                edge(0, 0, 3, [0, -1], 0, "m1"),
                edge(1, 1, 0, [-1, 0], -1, "m0"),
                edge(2, 0, 1, [-1, 0], 0, "m0"),
                edge(3, 2, 0, [0, -1], 1, "m0"),
                edge(4, 3, 2, [0, -1], 0, "m0"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: parent.loop_names.clone(),
            external_names: parent.external_names.clone(),
            node_name_to_internal: (0..4).map(|node| (format!("r{node}"), node)).collect(),
        };
        let residual_generated =
            generate_3d_expression(&residual, &Generate3DExpressionOptions::default()).unwrap();
        let bounded_residual_generated = generate_3d_expression(
            &residual,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2), (1, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let parent_input = EvaluationInput {
            external_momenta: vec![[1.3, 0.11, -0.17, 0.07]],
            loop_spatial_momenta: vec![[0.13, -0.09, 0.21], [-0.16, 0.19, 0.08]],
            masses: vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            uniform_scale: None,
        };
        let production_frame = |generated: &crate::GeneratedThreeDExpression| {
            generated.core_global_prefactor_sign.factor() as f64
        };
        let parent_raw_value = evaluate_expression(
            &parent,
            &parent_generated.expression,
            "dot(edges[0],edges[0])**2",
            &parent_input,
        )
        .unwrap()
        .value;
        let parent_value = parent_raw_value * production_frame(&parent_generated);
        let quadratic_parent_generated = generate_3d_expression(
            &parent,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let quadratic_parent_raw_value = evaluate_expression(
            &parent,
            &quadratic_parent_generated.expression,
            "dot(edges[0],edges[0])",
            &parent_input,
        )
        .unwrap()
        .value;
        let quadratic_parent_value =
            quadratic_parent_raw_value * production_frame(&quadratic_parent_generated);
        let residual_input = EvaluationInput {
            masses: parent_input.masses[1..].to_vec(),
            ..parent_input.clone()
        };
        let scalar_residual_raw_value = evaluate_expression(
            &residual,
            &residual_generated.expression,
            "1",
            &residual_input,
        )
        .unwrap()
        .value;
        let scalar_residual_value =
            scalar_residual_raw_value * production_frame(&residual_generated);
        let residual_value = evaluate_expression(
            &residual,
            &residual_generated.expression,
            "(-loops[0][0]+loops[1][0]-ext[0][0])**2\
             -(-loops[0][1]+loops[1][1]-ext[0][1])**2\
             -(-loops[0][2]+loops[1][2]-ext[0][2])**2\
             -(-loops[0][3]+loops[1][3]-ext[0][3])**2",
            &residual_input,
        )
        .unwrap()
        .value
            * production_frame(&residual_generated);
        let bounded_residual_raw_value = evaluate_expression(
            &residual,
            &bounded_residual_generated.expression,
            "(-loops[0][0]+loops[1][0]-ext[0][0])**2\
             -(-loops[0][1]+loops[1][1]-ext[0][1])**2\
             -(-loops[0][2]+loops[1][2]-ext[0][2])**2\
             -(-loops[0][3]+loops[1][3]-ext[0][3])**2",
            &residual_input,
        )
        .unwrap()
        .value;
        let bounded_residual_value =
            bounded_residual_raw_value * production_frame(&bounded_residual_generated);
        let quadratic_scale = quadratic_parent_value
            .abs()
            .max(scalar_residual_value.abs())
            .max(f64::MIN_POSITIVE);
        assert!(
            (quadratic_parent_value - scalar_residual_value).abs() <= 1.0e-11 * quadratic_scale,
            "the quadratic D0/D0 control must reproduce the scalar lower sector: parent={quadratic_parent_value:e} (raw={quadratic_parent_raw_value:e}, frame={}), residual={scalar_residual_value:e} (raw={scalar_residual_raw_value:e}, frame={})",
            production_frame(&quadratic_parent_generated),
            production_frame(&residual_generated),
        );
        let residual_scale = bounded_residual_value
            .abs()
            .max(residual_value.abs())
            .max(f64::MIN_POSITIVE);
        assert!(
            (bounded_residual_value - residual_value).abs() <= 1.0e-11 * residual_scale,
            "the bounded quadratic residual must agree with its convergent pure CFF: bounded={bounded_residual_value:e} (raw={bounded_residual_raw_value:e}, frame={}), pure={residual_value:e}",
            production_frame(&bounded_residual_generated),
        );
        let scale = parent_value
            .abs()
            .max(residual_value.abs())
            .max(f64::MIN_POSITIVE);

        assert!(
            (parent_value - residual_value).abs() <= 1.0e-11 * scale,
            "D0^2/D0 rank-four contact failed to reproduce the residual D0 numerator: parent={parent_value:e} (raw={parent_raw_value:e}, frame={}), residual={residual_value:e}",
            production_frame(&parent_generated),
        );
    }

    #[test]
    fn independent_bubbles_bounded_maps_match_component_cff_product() {
        let edge = |edge_id: usize,
                    tail: usize,
                    head: usize,
                    loop_signature: [i32; 2],
                    external_signature: i32,
                    mass_key: &str| {
            ParsedGraphInternalEdge {
                edge_id,
                tail,
                head,
                label: format!("q{edge_id}"),
                mass_key: Some(mass_key.to_string()),
                signature: MomentumSignature {
                    loop_signature: loop_signature.to_vec(),
                    external_signature: vec![external_signature],
                },
                had_pow: false,
            }
        };
        // Incidence joins the two bubbles at node zero, while their rational
        // loop-energy rows form independent x and y components.
        let parsed = ParsedGraph {
            internal_edges: vec![
                edge(0, 0, 1, [1, 0], 0, "m0"),
                edge(1, 1, 0, [1, 0], 1, "m1"),
                edge(2, 0, 2, [0, 1], 0, "m2"),
                edge(3, 2, 0, [0, 1], -1, "m3"),
            ],
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["x".to_string(), "y".to_string()],
            external_names: vec!["p".to_string()],
            node_name_to_internal: (0..3).map(|node| (format!("v{node}"), node)).collect(),
        };
        let component_local = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let cross_component = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2), (2, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let isolated_bubble =
            |loop_name: &str, external_sign: i32, masses: [&str; 2]| ParsedGraph {
                internal_edges: vec![
                    ParsedGraphInternalEdge {
                        edge_id: 0,
                        tail: 0,
                        head: 1,
                        label: "q0".to_string(),
                        mass_key: Some(masses[0].to_string()),
                        signature: MomentumSignature {
                            loop_signature: vec![1],
                            external_signature: vec![0],
                        },
                        had_pow: false,
                    },
                    ParsedGraphInternalEdge {
                        edge_id: 1,
                        tail: 1,
                        head: 0,
                        label: "q1".to_string(),
                        mass_key: Some(masses[1].to_string()),
                        signature: MomentumSignature {
                            loop_signature: vec![1],
                            external_signature: vec![external_sign],
                        },
                        had_pow: false,
                    },
                ],
                external_edges: Vec::new(),
                initial_state_cut_edges: Vec::new(),
                loop_names: vec![loop_name.to_string()],
                external_names: vec!["p".to_string()],
                node_name_to_internal: (0..2).map(|node| (format!("v{node}"), node)).collect(),
            };
        let x_bubble = isolated_bubble("x", 1, ["m0", "m1"]);
        let y_bubble = isolated_bubble("y", -1, ["m2", "m3"]);
        let x_bounded = generate_3d_expression(
            &x_bubble,
            &Generate3DExpressionOptions {
                energy_degree_bounds: Some(vec![(0, 2)]),
                ..Default::default()
            },
        )
        .unwrap();
        let y_scalar =
            generate_3d_expression(&y_bubble, &Generate3DExpressionOptions::default()).unwrap();
        let component_local_origins = component_local
            .expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .filter_map(|variant| variant.origin.as_deref())
            .collect::<std::collections::BTreeSet<_>>();
        for sample in ["e0=-", "e0=0", "e0=+"] {
            assert!(
                component_local_origins
                    .iter()
                    .any(|origin| { origin.contains("contact") && origin.ends_with(sample) }),
                "the component-local quadratic envelope lacks its deterministic {sample} bounded pinch branch: {component_local_origins:?}"
            );
        }
        assert!(
            component_local
                .energy_factor_components
                .iter()
                .all(|component| component.ownership
                    == crate::CffEnergyFactorOwnership::VariantLocal),
            "the component-local quadratic envelope must remain variant-local"
        );
        assert!(
            cross_component
                .energy_factor_components
                .iter()
                .all(|component| component.ownership
                    == crate::CffEnergyFactorOwnership::VariantLocal),
            "the cross-component envelope must retain variant-local bounded maps"
        );
        let x_edges = std::collections::BTreeSet::from([0, 1]);
        let y_edges = std::collections::BTreeSet::from([2, 3]);
        assert!(
            cross_component
                .expression
                .orientations
                .iter()
                .all(|orientation| {
                    orientation.variants.iter().all(|variant| {
                        let support = variant
                            .denominator_edges
                            .iter()
                            .map(|edge| edge.0)
                            .collect::<std::collections::BTreeSet<_>>();
                        !support.is_disjoint(&x_edges) && !support.is_disjoint(&y_edges)
                    })
                }),
            "a cross-component variant lost one factor of the exact component product"
        );

        let mut failures = Vec::new();
        let frame = |generated: &crate::GeneratedThreeDExpression| {
            generated.core_global_prefactor_sign.factor() as f64
        };
        for seed in [17, 1337, 9100] {
            let input = EvaluationInput::deterministic(
                &parsed,
                seed,
                &BTreeMap::from([
                    ("m0".to_string(), 0.43),
                    ("m1".to_string(), 0.71),
                    ("m2".to_string(), 0.59),
                    ("m3".to_string(), 0.83),
                ]),
                None,
            )
            .unwrap();
            let x_input = EvaluationInput {
                external_momenta: input.external_momenta.clone(),
                loop_spatial_momenta: vec![input.loop_spatial_momenta[0]],
                masses: input.masses[..2].to_vec(),
                uniform_scale: input.uniform_scale,
            };
            let y_input = EvaluationInput {
                external_momenta: input.external_momenta.clone(),
                loop_spatial_momenta: vec![input.loop_spatial_momenta[1]],
                masses: input.masses[2..].to_vec(),
                uniform_scale: input.uniform_scale,
            };
            let x_raw = evaluate_expression(
                &x_bubble,
                &x_bounded.expression,
                "edges[0][0]*edges[0][0]",
                &x_input,
            )
            .unwrap()
            .value;
            let x_framed = x_raw * frame(&x_bounded);
            let y_raw = evaluate_expression(&y_bubble, &y_scalar.expression, "1", &y_input)
                .unwrap()
                .value;
            let y_framed = y_raw * frame(&y_scalar);
            let full_component_raw = evaluate_expression(
                &parsed,
                &component_local.expression,
                "edges[0][0]*edges[0][0]",
                &input,
            )
            .unwrap()
            .value;
            let full_component_framed = full_component_raw * frame(&component_local);
            let isolated_product = x_framed * y_framed;
            let scale = isolated_product
                .abs()
                .max(full_component_framed.abs())
                .max(f64::MIN_POSITIVE);
            if (isolated_product - full_component_framed).abs() > 1.0e-11 * scale {
                failures.push(format!(
                    "isolated generated product differs from the full generated CFF at seed {seed}: isolated={isolated_product:.17e}, full={full_component_framed:.17e}, relative delta={:.17e}",
                    (isolated_product - full_component_framed).abs() / scale,
                ));
            }
            let mut expanded_generated = 0.0;
            for numerator in [
                "edges[0][0]*edges[0][0]",
                "2*edges[0][0]*edges[2][0]",
                "edges[2][0]*edges[2][0]",
            ] {
                let raw =
                    evaluate_expression(&parsed, &cross_component.expression, numerator, &input)
                        .unwrap()
                        .value;
                expanded_generated += raw * frame(&cross_component);
            }
            let factorized_generated = evaluate_expression(
                &parsed,
                &cross_component.expression,
                "(edges[0][0]+edges[2][0])*(edges[0][0]+edges[2][0])",
                &input,
            )
            .unwrap()
            .value
                * frame(&cross_component);
            let scale = expanded_generated
                .abs()
                .max(factorized_generated.abs())
                .max(f64::MIN_POSITIVE);
            if (expanded_generated - factorized_generated).abs() > 1.0e-11 * scale {
                failures.push(format!(
                    "cross-component generated CFF changed under numerator linearity at seed {seed}: expanded={expanded_generated:.17e}, factorized={factorized_generated:.17e}, relative delta={:.17e}",
                    (expanded_generated - factorized_generated).abs() / scale,
                ));
            }
        }
        assert!(
            failures.is_empty(),
            "independent-bubble bounded maps violate their CFF product invariants:\n{}",
            failures.join("\n"),
        );
    }

    #[test]
    fn reversed_vacuum_triangle_incidence_preserves_generalized_cff() {
        let triangle = |reversed: bool| ParsedGraph {
            internal_edges: (0..3)
                .map(|edge_id| {
                    let tail = edge_id;
                    let head = (edge_id + 1) % 3;
                    ParsedGraphInternalEdge {
                        edge_id,
                        tail: if reversed { head } else { tail },
                        head: if reversed { tail } else { head },
                        label: format!("q{edge_id}"),
                        mass_key: Some("m".to_string()),
                        signature: MomentumSignature {
                            loop_signature: vec![1],
                            external_signature: Vec::new(),
                        },
                        had_pow: false,
                    }
                })
                .collect(),
            external_edges: Vec::new(),
            initial_state_cut_edges: Vec::new(),
            loop_names: vec!["q".to_string()],
            external_names: Vec::new(),
            node_name_to_internal: (0..3).map(|node| (format!("n{node}"), node)).collect(),
        };
        let options = Generate3DExpressionOptions {
            energy_degree_bounds: Some(vec![(0, 2)]),
            ..Default::default()
        };
        let input = EvaluationInput {
            external_momenta: Vec::new(),
            loop_spatial_momenta: vec![[0.23, -0.31, 0.47]],
            masses: vec![1.0; 3],
            uniform_scale: None,
        };
        let evaluate = |parsed: &ParsedGraph| {
            let generated = generate_3d_expression(parsed, &options).unwrap();
            evaluate_expression(parsed, &generated.expression, "edges[0][0]**2", &input)
                .unwrap()
                .value
                * generated.core_global_prefactor_sign.factor() as f64
        };
        let forward = evaluate(&triangle(false));
        let reversed = evaluate(&triangle(true));
        let scale = forward.abs().max(reversed.abs()).max(f64::MIN_POSITIVE);
        assert!(
            (forward - reversed).abs() <= 1.0e-12 * scale,
            "reversing only the source incidence of an even vacuum triangle changed its generalized CFF: forward={forward:e}, reversed={reversed:e}"
        );
    }
}
