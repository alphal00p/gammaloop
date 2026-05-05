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

#[derive(Debug, Clone)]
pub struct ComparisonResult {
    pub cff: f64,
    pub ltd: f64,
    pub abs_cff_minus_ltd: f64,
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

pub struct ComparisonRequest<'a> {
    pub parsed: &'a ParsedGraph,
    pub cff_options: &'a crate::Generate3DExpressionOptions,
    pub ltd_options: &'a crate::Generate3DExpressionOptions,
    pub numerator_expr: &'a str,
    pub input: Option<EvaluationInput>,
    pub seed: u64,
    pub external_override: Option<Vec<[f64; 4]>>,
    pub loop_override: Option<Vec<[f64; 3]>>,
    pub mass_overrides: &'a BTreeMap<String, f64>,
}

pub fn compare_cff_ltd(request: ComparisonRequest<'_>) -> Result<ComparisonResult> {
    let input = match request.input {
        Some(input) => input,
        None => EvaluationInput::deterministic(
            request.parsed,
            request.seed,
            request.mass_overrides,
            None,
        )?
        .with_overrides(
            request.external_override,
            request.loop_override,
            request.mass_overrides,
            request.parsed,
        )?,
    };
    let cff = crate::generate_3d_expression_from_parsed(request.parsed, request.cff_options)
        .map_err(|error| EvaluationError::NumeratorEval(error.to_string()))?;
    let ltd = crate::generate_3d_expression_from_parsed(request.parsed, request.ltd_options)
        .map_err(|error| EvaluationError::NumeratorEval(error.to_string()))?;
    let cff_value =
        evaluate_expression(request.parsed, &cff, request.numerator_expr, &input)?.value;
    let ltd_value =
        evaluate_expression(request.parsed, &ltd, request.numerator_expr, &input)?.value;
    Ok(ComparisonResult {
        cff: cff_value,
        ltd: ltd_value,
        abs_cff_minus_ltd: (cff_value - ltd_value).abs(),
    })
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
            for variant in &orientation.variants {
                let mut denom = rational_to_f64(&variant.prefactor);
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
                total += denom * numerator_surface_factor * numerator_value;
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
        Generate3DExpressionOptions, NumeratorSamplingScaleMode, RepresentationMode,
        generate_3d_expression_from_parsed,
    };

    use super::*;

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
    fn evaluator_requires_nonzero_uniform_scale_when_expression_uses_m() {
        let parsed = crate::graph_io::test_graphs::box_pow3_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Ltd,
                energy_degree_bounds: vec![(3, 4)],
                numerator_sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
                include_cff_duplicate_signature_excess_sign: true,
                preserve_internal_edges_as_four_d_denominators: Vec::new(),
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
}
