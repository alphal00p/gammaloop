//! Joint movement in independent native coordinates, certified against every
//! original separator. Grouped coordinates receive the mean *raw* force; pins
//! and affine shifts are eliminated before solving the convex projection.

use std::collections::{BTreeMap, BTreeSet};

use clarabel::algebra::CscMatrix;
use clarabel::solver::{DefaultSettings, DefaultSolver, IPSolver, NonnegativeConeT};
use serde::{Deserialize, Serialize};
use twofloat::TwoFloat as Wide;

use super::super::spring::ShiftDirection;

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum AxisConstraint {
    Free,
    Fixed,
    Grouped {
        point: usize,
        direction: ShiftDirection,
    },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PointRecord {
    pub constraints: [AxisConstraint; 2],
    pub shift: [f64; 2],
    pub reference: [f64; 2],
    pub fixed_axes: [bool; 2],
    pub is_node: bool,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Halfplane {
    pub point: usize,
    pub normal: [f64; 2],
    pub bound: f64,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Certificate {
    pub status: String,
    pub primal_residual: f64,
    pub stationarity_residual: f64,
    pub complementarity_residual: f64,
    pub iterations: u32,
    pub rows: usize,
    pub components: usize,
}

#[derive(Clone, Debug)]
pub struct MovementDofs {
    index: Vec<[Option<usize>; 2]>,
    weights: Vec<f64>,
    fixed: Vec<[bool; 2]>,
    reference: Vec<[f64; 2]>,
    relations: Vec<(usize, usize, usize, f64)>,
    signs: Vec<(usize, usize, f64, f64, f64)>,
}

struct Components(Vec<usize>);

impl Components {
    fn new(count: usize) -> Self {
        Self((0..count).collect())
    }
    fn root(&mut self, mut index: usize) -> usize {
        while self.0[index] != index {
            self.0[index] = self.0[self.0[index]];
            index = self.0[index];
        }
        index
    }
    fn join(&mut self, first: usize, second: usize) {
        let first = self.root(first);
        let second = self.root(second);
        self.0[first] = second;
    }
}

impl MovementDofs {
    pub fn new(records: &[PointRecord], nodes_fixed: bool) -> Result<Self, String> {
        let mut components = Components::new(2 * records.len());
        let mut result = Self {
            index: vec![[None; 2]; records.len()],
            weights: Vec::new(),
            fixed: records.iter().map(|r| r.fixed_axes).collect(),
            reference: records.iter().map(|r| r.reference).collect(),
            relations: Vec::new(),
            signs: Vec::new(),
        };
        for (index, record) in records.iter().enumerate() {
            if !record
                .reference
                .iter()
                .chain(&record.shift)
                .all(|x| x.is_finite())
            {
                return Err("Native reference positions and shifts must be finite pairs".into());
            }
            for axis in 0..2 {
                match record.constraints[axis] {
                    AxisConstraint::Free => {}
                    AxisConstraint::Fixed => result.fixed[index][axis] = true,
                    AxisConstraint::Grouped { point, direction } => {
                        let other = records.get(point).ok_or("Unknown native group reference")?;
                        if nodes_fixed && record.is_node {
                            continue;
                        }
                        components.join(2 * index + axis, 2 * point + axis);
                        result.relations.push((
                            index,
                            point,
                            axis,
                            record.shift[axis] - other.shift[axis],
                        ));
                        if direction != ShiftDirection::Any && !(nodes_fixed && other.is_node) {
                            let sign = if direction == ShiftDirection::PositiveOnly {
                                1.0
                            } else {
                                -1.0
                            };
                            let initial = sign * (record.reference[axis] - record.shift[axis]);
                            if initial <= 0.0 {
                                return Err("Native directed reference has the wrong sign".into());
                            }
                            let floor = (0.5 * initial).min(
                                1e-9 * record.reference[axis]
                                    .abs()
                                    .max(record.shift[axis].abs())
                                    .max(1.0),
                            );
                            result
                                .signs
                                .push((index, axis, sign, record.shift[axis], floor));
                        }
                    }
                }
            }
        }
        let mut frozen = BTreeSet::new();
        for (index, axes) in result.fixed.iter().enumerate() {
            for (axis, fixed) in axes.iter().enumerate() {
                if *fixed {
                    frozen.insert(components.root(2 * index + axis));
                }
            }
        }
        let mut identifiers = BTreeMap::new();
        for (index, axes) in result.index.iter_mut().enumerate() {
            for (axis, value) in axes.iter_mut().enumerate() {
                let root = components.root(2 * index + axis);
                if frozen.contains(&root) {
                    continue;
                }
                let next = identifiers.len();
                let dof = *identifiers.entry(root).or_insert(next);
                *value = Some(dof);
                if result.weights.len() == dof {
                    result.weights.push(0.0);
                }
                result.weights[dof] += 1.0;
            }
        }
        result.validate(&result.reference)?;
        Ok(result)
    }

    pub fn validate(&self, positions: &[[f64; 2]]) -> Result<(), String> {
        if positions.len() != self.index.len() || !positions.iter().flatten().all(|x| x.is_finite())
        {
            return Err("Movement coordinates must be finite and match native keys".into());
        }
        for (index, fixed) in self.fixed.iter().enumerate() {
            for axis in 0..2 {
                if fixed[axis]
                    && (positions[index][axis] - self.reference[index][axis]).abs() > 1e-9
                {
                    return Err("Native fixed reference changed".into());
                }
            }
        }
        for &(first, second, axis, offset) in &self.relations {
            if (positions[first][axis] - positions[second][axis] - offset).abs() > 1e-9 {
                return Err("Native grouped coordinates disagree with their shifts".into());
            }
        }
        for &(index, axis, sign, shift, _) in &self.signs {
            if sign * (positions[index][axis] - shift) <= 0.0 {
                return Err("Native directed coordinate crossed its sign boundary".into());
            }
        }
        Ok(())
    }

    pub fn sign_rows(&self, positions: &[[f64; 2]]) -> Vec<Halfplane> {
        self.signs
            .iter()
            .map(|&(point, axis, sign, shift, floor)| {
                let mut normal = [0.0; 2];
                normal[axis] = -sign;
                Halfplane {
                    point,
                    normal,
                    bound: sign * (positions[point][axis] - shift) - floor,
                }
            })
            .collect()
    }

    /// A halfplane touches at most two independent native coordinates. Preserve
    /// the axis order and zero additions when constructing its sparse row.
    #[inline]
    fn coefficients(&self, row: &Halfplane) -> [(usize, f64); 2] {
        let mut coefficients = [(self.weights.len(), 0.0); 2];
        for (axis, &normal) in row.normal.iter().enumerate() {
            if let Some(column) = self.index[row.point][axis] {
                let slot = if coefficients[0].0 == column { 0 } else { axis };
                coefficients[slot].0 = column;
                coefficients[slot].1 += normal;
            }
        }
        coefficients
    }

    pub fn project(
        &self,
        proposal: &[[f64; 2]],
        rows: &[Halfplane],
        scale: f64,
    ) -> Result<(Vec<[f64; 2]>, Certificate), String> {
        if !scale.is_finite()
            || scale <= 0.0
            || proposal.len() != self.index.len()
            || !proposal.iter().flatten().all(|x| x.is_finite())
        {
            return Err("Nonfinite force or invalid movement scale".into());
        }
        let count = self.weights.len();
        let mut desired = vec![0.0; count];
        for (axes, movement) in self.index.iter().zip(proposal) {
            for axis in 0..2 {
                if let Some(column) = axes[axis] {
                    desired[column] += movement[axis];
                }
            }
        }
        for (value, weight) in desired.iter_mut().zip(&self.weights) {
            *value /= weight * scale;
        }
        if !desired.iter().all(|value| value.is_finite()) {
            return Err("Nonfinite normalized movement problem".into());
        }
        // Validate and test the grouped proposal before allocating sparse matrix
        // storage. Most late-epoch proposals already satisfy every halfplane.
        let mut components = Components::new(count);
        let mut joined = vec![false; self.index.len()];
        let mut feasible = true;
        let mut evidence = Certificate {
            status: "independent convex movement components".into(),
            ..Certificate::default()
        };
        for row in rows {
            if row.point >= self.index.len()
                || !row.bound.is_finite()
                || !row.normal.iter().all(|x| x.is_finite())
            {
                return Err("Nonfinite or invalid movement halfplane".into());
            }
            if row.bound < -1e-12 * scale {
                return Err("Old drawing is outside a movement halfplane".into());
            }
            let coefficients = self.coefficients(row);
            // Every row for this point couples the same two DOFs. Repeated joins
            // cannot change their component; axis-aligned rows do not join them.
            if !joined[row.point] && coefficients.iter().all(|(_, value)| *value != 0.0) {
                components.join(coefficients[0].0, coefficients[1].0);
                joined[row.point] = true;
            }
            if coefficients.iter().any(|(_, value)| *value != 0.0) {
                let bound = row.bound.max(0.0) / scale;
                if !bound.is_finite() {
                    return Err("Nonfinite normalized movement problem".into());
                }
                evidence.rows += 1;
                if feasible {
                    let terms = coefficients.map(|(column, coefficient)| {
                        coefficient * desired.get(column).copied().unwrap_or(0.0)
                    });
                    let values = if coefficients[0].0 <= coefficients[1].0 {
                        terms
                    } else {
                        [terms[1], terms[0]]
                    };
                    feasible = values.into_iter().sum::<f64>() <= bound;
                }
            }
        }
        if count == 0 {
            evidence.status = "all coordinates fixed".into();
            return Ok((vec![[0.0; 2]; proposal.len()], evidence));
        }
        if feasible {
            let mut roots = vec![false; count];
            for column in 0..count {
                roots[components.root(column)] = true;
            }
            evidence.components = roots.into_iter().filter(|root| *root).count();
            return Ok((
                self.index
                    .iter()
                    .map(|axes| {
                        axes.map(|column| column.map_or(0.0, |column| desired[column] * scale))
                    })
                    .collect(),
                evidence,
            ));
        }
        // Only constrained proposals need an explicit matrix. These rows have
        // already been validated in the feasibility pass above.
        let mut matrix = Vec::with_capacity(evidence.rows);
        let mut bounds = Vec::with_capacity(evidence.rows);
        for row in rows {
            let coefficients = self.coefficients(row);
            if coefficients.iter().any(|(_, value)| *value != 0.0) {
                matrix.push(coefficients);
                bounds.push(row.bound.max(0.0) / scale);
            }
        }
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for index in 0..count {
            groups
                .entry(components.root(index))
                .or_default()
                .push(index);
        }
        evidence.components = groups.len();
        let groups: Vec<_> = groups.into_values().collect();
        let mut locations = vec![(0, 0); count];
        let mut blocks = Vec::with_capacity(groups.len());
        for (group, columns) in groups.iter().enumerate() {
            for (local, &column) in columns.iter().enumerate() {
                locations[column] = (group, local);
            }
            blocks.push((Vec::new(), Vec::new()));
        }
        // Most late-epoch proposals are already feasible. Test the sparse
        // rows in the same column order as the dense dot product, without
        // allocating a dense block for those unconstrained optima.
        let mut unconstrained = vec![true; groups.len()];
        for (row, &bound) in matrix.iter().zip(&bounds) {
            let &(column, _) = row.iter().find(|(_, value)| *value != 0.0).unwrap();
            let group = locations[column].0;
            if unconstrained[group] {
                unconstrained[group] = groups[group]
                    .iter()
                    .map(|&column| {
                        let coefficient = row
                            .iter()
                            .find(|(index, _)| *index == column)
                            .map_or(0.0, |(_, value)| *value);
                        coefficient * desired[column]
                    })
                    .sum::<f64>()
                    <= bound;
            }
        }
        for (row, bound) in matrix.iter().zip(bounds) {
            let &(column, _) = row.iter().find(|(_, value)| *value != 0.0).unwrap();
            let group = locations[column].0;
            if unconstrained[group] {
                continue;
            }
            let mut local = vec![0.0; groups[group].len()];
            for &(column, value) in row {
                if value != 0.0 {
                    local[locations[column].1] = value;
                }
            }
            blocks[group].0.push(local);
            blocks[group].1.push(bound);
        }
        let mut values = vec![0.0; count];
        for ((columns, (block, block_bounds)), unconstrained) in
            groups.iter().zip(blocks).zip(unconstrained)
        {
            if unconstrained {
                for &column in columns {
                    values[column] = desired[column] * scale;
                }
                continue;
            }
            let weights: Vec<_> = columns.iter().map(|&column| self.weights[column]).collect();
            let desired: Vec<_> = columns.iter().map(|&column| desired[column]).collect();
            let (candidate, local) =
                ProjectedComponent::solve(&block, &block_bounds, &weights, &desired)?;
            for (&column, value) in columns.iter().zip(candidate) {
                values[column] = value * scale;
            }
            evidence.primal_residual = evidence.primal_residual.max(local.primal_residual);
            evidence.stationarity_residual = evidence
                .stationarity_residual
                .max(local.stationarity_residual);
            evidence.complementarity_residual = evidence
                .complementarity_residual
                .max(local.complementarity_residual);
            evidence.iterations += local.iterations;
        }
        evidence.primal_residual *= scale;
        Ok((
            self.index
                .iter()
                .map(|axes| axes.map(|column| column.map_or(0.0, |column| values[column])))
                .collect(),
            evidence,
        ))
    }
}

struct ProjectedComponent<'a> {
    matrix: &'a [Vec<f64>],
    bounds: &'a [f64],
    quadratic: Vec<f64>,
    linear: Vec<f64>,
}

impl ProjectedComponent<'_> {
    fn solve(
        matrix: &[Vec<f64>],
        bounds: &[f64],
        quadratic: &[f64],
        desired: &[f64],
    ) -> Result<(Vec<f64>, Certificate), String> {
        if matrix
            .iter()
            .zip(bounds)
            .all(|(row, bound)| row.iter().zip(desired).map(|(a, x)| a * x).sum::<f64>() <= *bound)
        {
            return Ok((
                desired.to_vec(),
                Certificate {
                    status: "unconstrained optimum".into(),
                    rows: bounds.len(),
                    ..Certificate::default()
                },
            ));
        }
        let objective_scale = quadratic
            .iter()
            .zip(desired)
            .map(|(q, v)| (q * v).abs())
            .fold(1.0, f64::max);
        let problem = ProjectedComponent {
            matrix,
            bounds,
            quadratic: quadratic.iter().map(|q| q / objective_scale).collect(),
            linear: quadratic
                .iter()
                .zip(desired)
                .map(|(q, v)| -q * v / objective_scale)
                .collect(),
        };
        let mut candidate_residual = None;
        let (mut values, dual, mut status, iterations) = if desired.len() <= 2 {
            let candidate = (desired.len() == 2)
                .then(|| problem.polygon_candidate())
                .flatten()
                .filter(|(values, dual)| {
                    let residual = problem.certify(values, dual);
                    if Self::certified(residual) {
                        candidate_residual = Some(residual);
                        true
                    } else {
                        false
                    }
                });
            let (values, dual, status) = if let Some((values, dual)) = candidate {
                (values, dual, "certified binary64 movement polygon")
            } else {
                let (values, dual) = problem.small_polygon()?;
                (values, dual, "enumerated movement polygon")
            };
            (values, dual, status.to_owned(), 0)
        } else {
            let p = CscMatrix::new(
                desired.len(),
                desired.len(),
                (0..=desired.len()).collect(),
                (0..desired.len()).collect(),
                problem.quadratic.clone(),
            );
            let mut column_ptr = vec![0];
            let mut row_index = Vec::new();
            let mut nonzeros = Vec::new();
            for column in 0..desired.len() {
                for (index, row) in matrix.iter().enumerate() {
                    if row[column] != 0.0 {
                        row_index.push(index);
                        nonzeros.push(row[column]);
                    }
                }
                column_ptr.push(nonzeros.len());
            }
            let a = CscMatrix::new(matrix.len(), desired.len(), column_ptr, row_index, nonzeros);
            let settings = DefaultSettings {
                verbose: false,
                max_threads: 1,
                max_iter: 200,
                direct_solve_method: "qdldl".into(),
                tol_gap_abs: 1e-12,
                tol_gap_rel: 1e-12,
                tol_feas: 1e-12,
                static_regularization_constant: 0.0,
                dynamic_regularization_eps: 1e-15,
                dynamic_regularization_delta: 1e-10,
                iterative_refinement_abstol: 1e-15,
                iterative_refinement_reltol: 1e-15,
                iterative_refinement_max_iter: 30,
                ..DefaultSettings::default()
            };
            let mut solver = DefaultSolver::new(
                &p,
                &problem.linear,
                &a,
                bounds,
                &[NonnegativeConeT(bounds.len())],
                settings,
            )
            .map_err(|e| e.to_string())?;
            solver.solve();
            (
                solver.solution.x,
                solver.solution.z.into_iter().map(Wide::from).collect(),
                format!("Clarabel {:?}", solver.solution.status),
                solver.solution.iterations,
            )
        };
        let mut residual = candidate_residual.unwrap_or_else(|| problem.certify(&values, &dual));
        if desired.len() > 2 && !Self::certified(residual) {
            if let Some((candidate, multipliers)) = problem.polish(&values, &dual) {
                let refined = problem.certify(&candidate, &multipliers);
                if Self::certified(refined) {
                    values = candidate;
                    residual = refined;
                    status.push_str(" (active-face polished)");
                }
            }
        }
        if !Self::certified(residual) {
            return Err(format!(
                "Uncertified joint movement ({status}): primal={}, stationarity={}, complementarity={}, dual_sign={}",
                residual[0], residual[1], residual[2], residual[3]
            ));
        }
        Ok((
            values,
            Certificate {
                status,
                primal_residual: residual[0],
                stationarity_residual: residual[1],
                complementarity_residual: residual[2],
                iterations,
                rows: bounds.len(),
                components: 1,
            },
        ))
    }

    fn certified(values: [f64; 4]) -> bool {
        values
            .into_iter()
            .zip([2e-9, 2e-8, 2e-8, 2e-9])
            .all(|(value, limit)| value.is_finite() && value <= limit)
    }

    fn certify(&self, values: &[f64], dual: &[Wide]) -> [f64; 4] {
        if values.len() != self.quadratic.len()
            || dual.len() != self.bounds.len()
            || !values.iter().all(|v| v.is_finite())
            || !dual.iter().all(Wide::is_valid)
        {
            return [f64::INFINITY; 4];
        }
        let values: Vec<_> = values.iter().copied().map(Wide::from).collect();
        let mut stationarity: Vec<_> = self
            .quadratic
            .iter()
            .zip(&values)
            .zip(&self.linear)
            .map(|((&q, &v), &l)| Wide::from(q) * v + Wide::from(l))
            .collect();
        let mut residual = [0.0_f64; 4];
        let terms = (self.quadratic.len() + 1) as f64;
        let roundoff = 32.0 * f64::EPSILON * terms;
        for ((row, &bound), &multiplier) in self.matrix.iter().zip(self.bounds).zip(dual) {
            if multiplier == Wide::from(0.0) {
                // A zero dual contributes no stationarity or complementarity.
                // Prove its primal slack negative with a conservative sequential
                // dot-product error bound before using wide arithmetic. The
                // absolute floor covers gradual underflow; the magnitude guard
                // excludes overflow in the wide sum. Boundary or extreme rows
                // retain the full-precision calculation below.
                let mut approximate = -bound;
                let mut magnitude = bound.abs();
                for (&a, &x) in row.iter().zip(&values) {
                    let product = a * x.hi();
                    approximate += product;
                    magnitude += product.abs();
                }
                let allowance = roundoff * magnitude + terms * f64::MIN_POSITIVE;
                if roundoff < 0.125
                    && magnitude < f64::MAX / 4.0
                    && allowance.is_finite()
                    && approximate < -allowance
                {
                    continue;
                }
            }
            let excess = row
                .iter()
                .zip(&values)
                .fold(Wide::from(-bound), |sum, (&a, &x)| sum + Wide::from(a) * x);
            if !excess.is_valid() {
                return [f64::INFINITY; 4];
            }
            residual[0] = residual[0].max(f64::from(excess));
            if multiplier != Wide::from(0.0) {
                residual[2] = residual[2].max(f64::from((multiplier * excess).abs()));
                residual[3] = residual[3].max(f64::from(-multiplier));
                for (stationarity, &a) in stationarity.iter_mut().zip(row) {
                    *stationarity += Wide::from(a) * multiplier;
                }
            }
        }
        residual[1] = stationarity
            .into_iter()
            .map(|x| f64::from(x.abs()))
            .fold(0.0, f64::max);
        residual
    }

    // TwoFloat 0.8's reciprocal implementation rounds its high-word product
    // before subtracting one. Recover the quotient from successive extended
    // residuals instead, so nearly parallel boundary intersections retain both
    // words rather than silently reverting to binary64 division accuracy.
    fn divide(numerator: Wide, denominator: Wide) -> Wide {
        let first = Wide::from(numerator.hi() / denominator.hi());
        let remainder = numerator - denominator * first;
        let second = Wide::from(remainder.hi() / denominator.hi());
        let remainder = remainder - denominator * second;
        first + second + Wide::from(remainder.hi() / denominator.hi())
    }

    // Binary64 identifies a likely active face cheaply. This is only a
    // proposal: the extended-precision KKT certificate checks every original
    // halfplane and multiplier before it can replace the double-double solve.
    fn polygon_candidate(&self) -> Option<(Vec<f64>, Vec<Wide>)> {
        debug_assert_eq!(self.quadratic.len(), 2);
        let desired = [
            -self.linear[0] / self.quadratic[0],
            -self.linear[1] / self.quadratic[1],
        ];
        struct Candidate {
            cost: f64,
            point: [f64; 2],
            multipliers: [(usize, f64); 2],
        }
        let mut best: Option<Candidate> = None;
        let mut consider = |point: [f64; 2], multipliers: [(usize, f64); 2]| {
            if !point.iter().all(|x| x.is_finite())
                || multipliers
                    .iter()
                    .any(|&(_, value)| !value.is_finite() || value < -2e-9)
            {
                return;
            }
            // This solver is always two-dimensional. Keep the sequential
            // floating-point sum (including its -0.0 seed) without constructing
            // iterators for every candidate and every feasibility row.
            let cost = -0.0
                + (0.5 * self.quadratic[0] * point[0] * point[0] + self.linear[0] * point[0])
                + (0.5 * self.quadratic[1] * point[1] * point[1] + self.linear[1] * point[1]);
            if !cost.is_finite()
                || best
                    .as_ref()
                    .is_some_and(|candidate| cost >= candidate.cost)
            {
                return;
            }

            // A roundoff allowance finds boundary candidates; it does not
            // relax acceptance, which uses the original wide certificate.
            if self.matrix.iter().zip(self.bounds).any(|(row, &bound)| {
                let x = row[0] * point[0];
                let y = row[1] * point[1];
                let dot = -0.0 + x + y;
                let magnitude = -0.0 + x.abs() + y.abs();
                dot - bound > 32.0 * f64::EPSILON * magnitude.max(bound.abs()).max(1.0)
            }) {
                return;
            }
            best = Some(Candidate {
                cost,
                point,
                multipliers,
            });
        };
        for (index, (row, &bound)) in self.matrix.iter().zip(self.bounds).enumerate() {
            let distance = (-0.0 + row[0] * desired[0] + row[1] * desired[1]) - bound;
            let denominator =
                -0.0 + row[0] * row[0] / self.quadratic[0] + row[1] * row[1] / self.quadratic[1];
            let multiplier = distance / denominator;
            let point = [
                desired[0] - multiplier * row[0] / self.quadratic[0],
                desired[1] - multiplier * row[1] / self.quadratic[1],
            ];
            consider(point, [(index, multiplier), (index, 0.0)]);
        }
        for first in 0..self.matrix.len() {
            for second in first + 1..self.matrix.len() {
                let a = &self.matrix[first];
                let b = &self.matrix[second];
                let determinant = a[0] * b[1] - a[1] * b[0];
                if determinant == 0.0 {
                    continue;
                }
                let point = [
                    (self.bounds[first] * b[1] - a[1] * self.bounds[second]) / determinant,
                    (a[0] * self.bounds[second] - self.bounds[first] * b[0]) / determinant,
                ];
                let gradient = [
                    -(self.quadratic[0] * point[0] + self.linear[0]),
                    -(self.quadratic[1] * point[1] + self.linear[1]),
                ];
                consider(
                    point,
                    [
                        (
                            first,
                            (gradient[0] * b[1] - b[0] * gradient[1]) / determinant,
                        ),
                        (
                            second,
                            (a[0] * gradient[1] - gradient[0] * a[1]) / determinant,
                        ),
                    ],
                );
            }
        }
        let Candidate {
            point, multipliers, ..
        } = best?;
        let mut dual = vec![Wide::from(0.0); self.bounds.len()];
        for (row, value) in multipliers {
            dual[row] += Wide::from(value);
        }
        Some((point.to_vec(), dual))
    }

    fn small_polygon(&self) -> Result<(Vec<f64>, Vec<Wide>), String> {
        let zero = Wide::from(0.0);
        let mut dual = vec![zero; self.bounds.len()];
        let desired: Vec<_> = self
            .linear
            .iter()
            .zip(&self.quadratic)
            .map(|(&l, &q)| Self::divide(-Wide::from(l), Wide::from(q)))
            .collect();
        if desired.len() == 1 {
            let mut lower: Option<(Wide, usize)> = None;
            let mut upper: Option<(Wide, usize)> = None;
            for (index, (row, &bound)) in self.matrix.iter().zip(self.bounds).enumerate() {
                let value = Self::divide(Wide::from(bound), Wide::from(row[0]));
                if row[0] > 0.0 && upper.is_none_or(|(old, _)| value < old) {
                    upper = Some((value, index));
                }
                if row[0] < 0.0 && lower.is_none_or(|(old, _)| value > old) {
                    lower = Some((value, index));
                }
            }
            let mut value = desired[0];
            let mut active = None;
            if let Some((bound, row)) = lower {
                if value < bound {
                    value = bound;
                    active = Some(row);
                }
            }
            if let Some((bound, row)) = upper {
                if value > bound {
                    value = bound;
                    active = Some(row);
                }
            }
            if let Some(row) = active {
                dual[row] = Self::divide(
                    -(Wide::from(self.quadratic[0]) * value + Wide::from(self.linear[0])),
                    Wide::from(self.matrix[row][0]),
                );
            }
            return Ok((vec![f64::from(value)], dual));
        }
        let roots = [
            Wide::from(self.quadratic[0]).sqrt(),
            Wide::from(self.quadratic[1]).sqrt(),
        ];
        let target = [desired[0] * roots[0], desired[1] * roots[1]];
        let mut normals = Vec::with_capacity(self.matrix.len());
        let mut lengths = Vec::with_capacity(self.matrix.len());
        let mut offsets = Vec::with_capacity(self.matrix.len());
        for (row, &bound) in self.matrix.iter().zip(self.bounds) {
            let transformed = [
                Self::divide(Wide::from(row[0]), roots[0]),
                Self::divide(Wide::from(row[1]), roots[1]),
            ];
            let length = (transformed[0] * transformed[0] + transformed[1] * transformed[1]).sqrt();
            normals.push([
                Self::divide(transformed[0], length),
                Self::divide(transformed[1], length),
            ]);
            lengths.push(length);
            offsets.push(Self::divide(Wide::from(bound), length));
        }
        // Retain the accepted long-double feasibility allowance; double-double
        // arithmetic evaluates the predicate and final KKT residual more accurately.
        let roundoff = |point: [Wide; 2]| {
            Wide::from(128.0 * 2_f64.powi(-63))
                * point[0].abs().max(point[1].abs()).max(Wide::from(1.0))
        };
        let distance: Vec<_> = normals
            .iter()
            .zip(&offsets)
            .map(|(n, &b)| n[0] * target[0] + n[1] * target[1] - b)
            .collect();
        if distance.iter().all(|&d| d <= roundoff(target)) {
            return Ok((desired.into_iter().map(f64::from).collect(), dual));
        }
        struct Candidate {
            cost: Wide,
            point: [Wide; 2],
            multipliers: [Option<(usize, Wide)>; 2],
        }
        let cost = |point: [Wide; 2]| {
            Wide::from(0.5) * (point[0] * point[0] + point[1] * point[1])
                - point[0] * target[0]
                - point[1] * target[1]
        };
        let mut best: Option<Candidate> = None;
        let consider = |best: &mut Option<Candidate>,
                        point: [Wide; 2],
                        multipliers: [Option<(usize, Wide)>; 2]| {
            if !point.iter().all(Wide::is_valid)
                || multipliers
                    .iter()
                    .flatten()
                    .any(|(_, value)| !value.is_valid() || *value < Wide::from(-2e-9))
            {
                return;
            }
            let tolerance = roundoff(point);
            if normals.iter().zip(&offsets).any(|(normal, &bound)| {
                normal[0] * point[0] + normal[1] * point[1] - bound > tolerance
            }) {
                return;
            }
            let cost = cost(point);
            if best.as_ref().is_none_or(|candidate| cost < candidate.cost) {
                *best = Some(Candidate {
                    cost,
                    point,
                    multipliers,
                });
            }
        };
        for (index, normal) in normals.iter().enumerate() {
            consider(
                &mut best,
                [
                    target[0] - distance[index] * normal[0],
                    target[1] - distance[index] * normal[1],
                ],
                [
                    Some((index, Self::divide(distance[index], lengths[index]))),
                    None,
                ],
            );
        }
        for first in 0..normals.len() {
            for second in first + 1..normals.len() {
                let a = normals[first];
                let b = normals[second];
                let determinant = a[0] * b[1] - a[1] * b[0];
                if determinant == zero {
                    continue;
                }
                let point = [
                    Self::divide(offsets[first] * b[1] - a[1] * offsets[second], determinant),
                    Self::divide(a[0] * offsets[second] - offsets[first] * b[0], determinant),
                ];
                // The original enumeration only replaces a feasible candidate
                // with a strictly lower objective. A dominated intersection
                // cannot win, so it needs neither dual divisions nor a scan of
                // every halfplane. Keep the same extended-precision objective
                // and enumeration order, including ties.
                if best
                    .as_ref()
                    .is_some_and(|candidate| cost(point) >= candidate.cost)
                {
                    continue;
                }
                let delta = [target[0] - point[0], target[1] - point[1]];
                consider(
                    &mut best,
                    point,
                    [
                        Some((
                            first,
                            Self::divide(
                                Self::divide(delta[0] * b[1] - b[0] * delta[1], determinant),
                                lengths[first],
                            ),
                        )),
                        Some((
                            second,
                            Self::divide(
                                Self::divide(a[0] * delta[1] - delta[0] * a[1], determinant),
                                lengths[second],
                            ),
                        )),
                    ],
                );
            }
        }
        let Candidate {
            point, multipliers, ..
        } = best.ok_or("No primal/dual certified candidate in a movement polygon")?;
        for (index, value) in multipliers.into_iter().flatten() {
            dual[index] = value;
        }
        Ok((
            vec![
                f64::from(Self::divide(point[0], roots[0])),
                f64::from(Self::divide(point[1], roots[1])),
            ],
            dual,
        ))
    }

    fn polish(&self, values: &[f64], dual: &[Wide]) -> Option<(Vec<f64>, Vec<Wide>)> {
        if values.len() != self.quadratic.len()
            || dual.len() != self.bounds.len()
            || !values.iter().all(|v| v.is_finite())
            || !dual.iter().all(Wide::is_valid)
        {
            return None;
        }
        let mut active: Vec<_> = self
            .matrix
            .iter()
            .zip(self.bounds)
            .zip(dual)
            .enumerate()
            .filter_map(|(index, ((row, bound), multiplier))| {
                (*multiplier > Wide::from(1e-7)
                    && (row.iter().zip(values).map(|(a, b)| a * b).sum::<f64>() - bound).abs()
                        < 1e-7)
                    .then_some(index)
            })
            .collect();
        for _attempt in 0..=active.len() {
            let basis = self.independent_rows(&active);
            let count = self.quadratic.len();
            let size = count + basis.len();
            let mut system = vec![vec![0.0; size]; size];
            let mut right = vec![0.0; size];
            for index in 0..count {
                system[index][index] = self.quadratic[index];
                right[index] = -self.linear[index];
            }
            for (index, &row) in basis.iter().enumerate() {
                right[count + index] = self.bounds[row];
                let (coordinate_rows, multiplier_rows) = system.split_at_mut(count);
                for (column, coordinate_row) in coordinate_rows.iter_mut().enumerate() {
                    multiplier_rows[index][column] = self.matrix[row][column];
                    coordinate_row[count + index] = self.matrix[row][column];
                }
            }
            let factor = LuFactor::new(&system)?;
            let mut solution: Vec<_> = factor.solve(&right)?.into_iter().map(Wide::from).collect();
            for _ in 0..6 {
                let residual: Vec<_> = system
                    .iter()
                    .zip(&right)
                    .map(|(row, &b)| {
                        f64::from(
                            row.iter()
                                .zip(&solution)
                                .fold(Wide::from(b), |sum, (&a, &x)| sum - Wide::from(a) * x),
                        )
                    })
                    .collect();
                let correction = factor.solve(&residual)?;
                for (value, correction) in solution.iter_mut().zip(correction) {
                    *value += Wide::from(correction);
                }
            }
            if !solution.iter().all(Wide::is_valid) {
                return None;
            }
            let mut multipliers = vec![Wide::from(0.0); self.bounds.len()];
            for (&row, &value) in basis.iter().zip(&solution[count..]) {
                multipliers[row] = value;
            }
            let negative = multipliers
                .iter()
                .enumerate()
                .filter(|(_, x)| **x < Wide::from(-2e-9))
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
            if let Some((remove, _)) = negative {
                active.retain(|&row| row != remove);
                continue;
            }
            return Some((
                solution[..count].iter().copied().map(f64::from).collect(),
                multipliers,
            ));
        }
        None
    }

    /// Column-pivoted Householder QR of the active transpose selects the same
    /// independent-face criterion as the reference SciPy pivoted QR.
    fn independent_rows(&self, active: &[usize]) -> Vec<usize> {
        if active.is_empty() {
            return Vec::new();
        }
        let count = self.quadratic.len();
        let mut matrix: Vec<Vec<_>> = (0..count)
            .map(|column| active.iter().map(|&row| self.matrix[row][column]).collect())
            .collect();
        let mut pivot: Vec<_> = (0..active.len()).collect();
        let mut diagonals = Vec::new();
        for step in 0..count.min(active.len()) {
            let selected = (step..active.len())
                .max_by(|&a, &b| {
                    let norm = |column| {
                        matrix[step..]
                            .iter()
                            .map(|row| row[column] * row[column])
                            .sum::<f64>()
                    };
                    norm(a).total_cmp(&norm(b)).then_with(|| b.cmp(&a))
                })
                .unwrap();
            for row in &mut matrix {
                row.swap(step, selected);
            }
            pivot.swap(step, selected);
            let norm = matrix[step..]
                .iter()
                .fold(0.0_f64, |sum, row| sum.hypot(row[step]));
            let alpha = -norm.copysign(matrix[step][step]);
            diagonals.push(alpha.abs());
            if norm == 0.0 {
                continue;
            }
            let mut reflector: Vec<_> = matrix[step..].iter().map(|row| row[step]).collect();
            reflector[0] -= alpha;
            let denominator = reflector.iter().map(|value| value * value).sum::<f64>();
            for column in step..active.len() {
                let factor = 2.0
                    * matrix[step..]
                        .iter()
                        .zip(&reflector)
                        .map(|(row, &v)| row[column] * v)
                        .sum::<f64>()
                    / denominator;
                for (row, &value) in matrix[step..].iter_mut().zip(&reflector) {
                    row[column] -= factor * value;
                }
            }
            matrix[step][step] = alpha;
        }
        let maximum = matrix
            .iter()
            .enumerate()
            .flat_map(|(row, values)| values.iter().skip(row).map(|v| v.abs()))
            .fold(1.0, f64::max);
        let rank = diagonals
            .iter()
            .filter(|&&value| value > 1e-12 * maximum)
            .count();
        pivot[..rank].iter().map(|&index| active[index]).collect()
    }
}

/// Partial-pivot LU is reused for six extended-residual refinements. A singular
/// active face is rejected; no zero-movement or relaxed-constraint fallback exists.
struct LuFactor {
    matrix: Vec<Vec<f64>>,
    swaps: Vec<usize>,
}

impl LuFactor {
    fn new(input: &[Vec<f64>]) -> Option<Self> {
        let mut matrix = input.to_vec();
        let count = matrix.len();
        let mut swaps = Vec::new();
        for column in 0..count {
            let pivot = (column..count).max_by(|&a, &b| {
                matrix[a][column]
                    .abs()
                    .total_cmp(&matrix[b][column].abs())
                    .then_with(|| b.cmp(&a))
            })?;
            if matrix[pivot][column] == 0.0 || !matrix[pivot][column].is_finite() {
                return None;
            }
            swaps.push(pivot);
            matrix.swap(column, pivot);
            for row in column + 1..count {
                matrix[row][column] /= matrix[column][column];
                for entry in column + 1..count {
                    matrix[row][entry] -= matrix[row][column] * matrix[column][entry];
                }
            }
        }
        Some(Self { matrix, swaps })
    }
    fn solve(&self, right: &[f64]) -> Option<Vec<f64>> {
        let mut values = right.to_vec();
        let count = values.len();
        for (row, &swap) in self.swaps.iter().enumerate() {
            values.swap(row, swap);
        }
        for row in 0..count {
            for column in 0..row {
                values[row] -= self.matrix[row][column] * values[column];
            }
        }
        for row in (0..count).rev() {
            for column in row + 1..count {
                values[row] -= self.matrix[row][column] * values[column];
            }
            values[row] /= self.matrix[row][row];
        }
        values
            .iter()
            .all(|value| value.is_finite())
            .then_some(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn records(points: &[[f64; 2]]) -> Vec<PointRecord> {
        points
            .iter()
            .map(|&reference| PointRecord {
                constraints: [AxisConstraint::Free; 2],
                shift: [0.0; 2],
                reference,
                fixed_axes: [false; 2],
                is_node: false,
            })
            .collect()
    }

    fn box_rows(count: usize, cap: f64) -> Vec<Halfplane> {
        (0..count)
            .flat_map(|point| {
                [[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]].map(|normal| Halfplane {
                    point,
                    normal,
                    bound: cap,
                })
            })
            .collect()
    }

    #[test]
    fn blocked_axis_preserves_group_motion_and_other_free_axes() {
        let points = [[0.0, 0.0], [0.0, 1.0]];
        let mut metadata = records(&points);
        metadata[1].constraints[0] = AxisConstraint::Grouped {
            point: 0,
            direction: ShiftDirection::Any,
        };
        let dofs = MovementDofs::new(&metadata, false).unwrap();
        let mut rows = box_rows(2, 0.2);
        rows.push(Halfplane {
            point: 0,
            normal: [0.0, -1.0],
            bound: 0.0,
        });
        let (step, certificate) = dofs
            .project(&[[10.0, -0.2], [-1.0, 0.2]], &rows, 0.2)
            .unwrap();
        assert!(step[0][0] > 0.19);
        assert_eq!(step[0][0], step[1][0]);
        assert!(step[0][1] >= -1e-10);
        assert!(step[1][1] > 0.19);
        assert!(certificate.components > 1);
    }

    #[test]
    fn projection_keeps_sparse_axes_and_fixed_rows_separate() {
        let mut metadata = records(&[[0.0; 2]; 4]);
        metadata[1].constraints[0] = AxisConstraint::Fixed;
        metadata[2].constraints = [
            AxisConstraint::Grouped {
                point: 0,
                direction: ShiftDirection::Any,
            },
            AxisConstraint::Fixed,
        ];
        metadata[3].constraints = [AxisConstraint::Fixed; 2];
        let dofs = MovementDofs::new(&metadata, false).unwrap();
        let mut rows = box_rows(4, 0.125);
        rows.extend([
            Halfplane {
                point: 0,
                normal: [-0.0, 0.0],
                bound: 0.0,
            },
            Halfplane {
                point: 1,
                normal: [1.0, -0.0],
                bound: 0.05,
            },
            Halfplane {
                point: 3,
                normal: [1.0, 1.0],
                bound: 0.05,
            },
        ]);
        let (step, certificate) = dofs
            .project(
                &[[2.0, -3.0], [99.0, 2.0], [-1.0, 99.0], [99.0; 2]],
                &rows,
                1.0,
            )
            .unwrap();
        assert_eq!(
            step,
            [[0.125, -0.125], [0.0, 0.125], [0.125, 0.0], [0.0; 2]]
        );
        assert_eq!(certificate.rows, 8);
        assert_eq!(certificate.components, 3);
        assert_eq!(certificate.primal_residual, 0.0);
    }

    #[test]
    fn shifted_cycles_eliminate_pins_and_preserve_directed_signs() {
        let points = [[2.0, 0.0], [3.0, 2.0], [4.0, 4.0]];
        let mut metadata = records(&points);
        for (index, row) in metadata.iter_mut().enumerate() {
            row.shift[0] = index as f64 + 1.0;
            row.constraints[0] = AxisConstraint::Grouped {
                point: (index + 1) % 3,
                direction: ShiftDirection::PositiveOnly,
            };
        }
        metadata[1].constraints[1] = AxisConstraint::Fixed;
        let dofs = MovementDofs::new(&metadata, false).unwrap();
        let mut rows = box_rows(3, 2.0);
        rows.extend(dofs.sign_rows(&points));
        let (step, _) = dofs
            .project(&[[-100.0, 2.0], [-1.0, 4.0], [-1.0, -2.0]], &rows, 2.0)
            .unwrap();
        assert_eq!(step[0][0], step[1][0]);
        assert_eq!(step[1][0], step[2][0]);
        assert_eq!(step[1][1], 0.0);
        let mut after = points;
        for (point, step) in after.iter_mut().zip(step) {
            for axis in 0..2 {
                point[axis] += step[axis];
            }
        }
        dofs.validate(&after).unwrap();
        after[1][1] += 0.1;
        assert!(
            dofs.validate(&after)
                .unwrap_err()
                .contains("fixed reference")
        );
        metadata[0].shift[0] += 0.1;
        assert!(
            MovementDofs::new(&metadata, false)
                .unwrap_err()
                .contains("grouped coordinates")
        );
    }

    #[test]
    fn invalid_forces_bounds_and_group_targets_are_rejected() {
        let mut metadata = records(&[[0.0, 0.0]]);
        metadata[0].constraints[0] = AxisConstraint::Grouped {
            point: 5,
            direction: ShiftDirection::Any,
        };
        assert!(MovementDofs::new(&metadata, false).is_err());
        let dofs = MovementDofs::new(&records(&[[0.0, 0.0]]), false).unwrap();
        assert!(dofs.project(&[[f64::NAN, 0.0]], &[], 1.0).is_err());
        assert!(
            dofs.project(
                &[[0.0, 0.0]],
                &[Halfplane {
                    point: 0,
                    normal: [1.0, 0.0],
                    bound: -1.0
                }],
                1.0
            )
            .is_err()
        );
        assert!(dofs.project(&[[0.0, 0.0]], &[], 0.0).is_err());
    }

    #[test]
    fn inactive_certificate_rows_keep_full_precision_primal_residuals() {
        // A nonzero dual disables the inactive-row filter without changing the
        // primal calculation. Compare that wide reference across cancellation,
        // subnormal products, signed zero, and overflow-risk scales.
        let cases = [
            ([1.0, 0.0], [1.0, 1.0], 2.0),
            ([1.0, 0.0], [1.0, 1.0], 1.0 - 4e-9),
            ([1e300, -1e300], [1.0, 1.0], 1e285),
            ([1e300, -1e300], [1.0, 1.0], -1e285),
            (
                [f64::MIN_POSITIVE, 0.0],
                [f64::MIN_POSITIVE, 1.0],
                -f64::from_bits(1),
            ),
            ([-0.0, 0.0], [1.0, -1.0], -0.0),
            ([f64::MAX / 2.0, -f64::MAX / 2.0], [1.0, 1.0], 1.0),
            ([f64::MAX, 0.0], [2.0, 1.0], 1.0),
        ];
        for (row, values, bound) in cases {
            let matrix = [row.to_vec()];
            let bounds = [bound];
            let problem = ProjectedComponent {
                matrix: &matrix,
                bounds: &bounds,
                quadratic: vec![1.0; 2],
                linear: values.map(|value| -value).to_vec(),
            };
            let filtered = problem.certify(&values, &[Wide::from(0.0)]);
            let wide = problem.certify(&values, &[Wide::from(f64::from_bits(1))]);
            assert_eq!(
                filtered[0], wide[0],
                "row={row:?}, values={values:?}, bound={bound}"
            );
        }
    }

    #[test]
    fn binary64_polygon_requires_the_extended_certificate() {
        let matrix = vec![
            vec![1.0, 0.0],
            vec![-1.0, 0.0],
            vec![0.0, 1.0],
            vec![0.0, -1.0],
        ];
        let (values, certificate) =
            ProjectedComponent::solve(&matrix, &[1.0; 4], &[2.0, 3.0], &[2.0, 3.0]).unwrap();
        assert_eq!(values, [1.0, 1.0]);
        assert_eq!(certificate.status, "certified binary64 movement polygon");
        assert!(certificate.primal_residual <= 2e-9);
        assert!(certificate.stationarity_residual <= 2e-8);
        assert!(certificate.complementarity_residual <= 2e-8);
    }

    #[test]
    fn binary64_polygon_preserves_candidate_bits_at_small_scales_and_ties() {
        // Frozen before specializing the two-coordinate sums. The repeated
        // boundary also checks that equally good candidates keep their order.
        let cases = [
            (
                vec![
                    vec![1.0, 1e-150],
                    vec![-1.0, 0.0],
                    vec![1e-150, 1.0],
                    vec![0.0, -1.0],
                ],
                vec![1e-150; 4],
                vec![1.0, 2.0],
                vec![-2e-150, 3e-150],
                [0x20da2fe76a3f9475, 0xa0ca2fe76a3f9475],
            ),
            (
                vec![
                    vec![1.0, 1.0],
                    vec![1.0, 1.0],
                    vec![1.0, -1.0],
                    vec![-1.0, 1.0],
                    vec![-1.0, -1.0],
                ],
                vec![0.5; 5],
                vec![0.1, 0.3],
                vec![-1.0, -1.0],
                [0x3fd8000000000000, 0x3fc0000000000000],
            ),
        ];
        for (matrix, bounds, quadratic, linear, expected) in cases {
            let problem = ProjectedComponent {
                matrix: &matrix,
                bounds: &bounds,
                quadratic,
                linear,
            };
            let (candidate, _) = problem.polygon_candidate().unwrap();
            assert_eq!(
                candidate.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                expected
            );
        }
    }

    #[test]
    fn nearly_parallel_polygon_escalates_precision_for_amplified_duals() {
        let matrix = vec![vec![1.0, 1.0], vec![-1.0, -1.0 + 1e-10]];
        let problem = ProjectedComponent {
            matrix: &matrix,
            bounds: &[0.0; 2],
            quadratic: vec![1.0; 2],
            linear: vec![-0.123, -1.0],
        };
        let (candidate, dual) = problem.polygon_candidate().unwrap();
        assert_eq!(candidate[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(candidate[1].to_bits(), 0.0_f64.to_bits());
        assert!(!ProjectedComponent::certified(
            problem.certify(&candidate, &dual)
        ));
        let (values, certificate) =
            ProjectedComponent::solve(&matrix, &[0.0; 2], &[1.0; 2], &[0.123, 1.0]).unwrap();
        assert_eq!(values, [0.0; 2]);
        assert_eq!(certificate.status, "enumerated movement polygon");
        assert!(certificate.primal_residual <= 2e-9);
        assert!(certificate.stationarity_residual <= 2e-8);
        assert!(certificate.complementarity_residual <= 2e-8);
    }

    #[test]
    fn frozen_thin_regions_and_amplified_multipliers_keep_the_original_certificates() {
        let fixtures: Vec<serde_json::Value> =
            serde_json::from_str(include_str!("fixtures/movement-regressions.json")).unwrap();
        for fixture in fixtures {
            let name = fixture["name"].as_str().unwrap();
            if fixture.get("matrix").is_some() {
                let matrix: Vec<Vec<f64>> =
                    serde_json::from_value(fixture["matrix"].clone()).unwrap();
                let bounds: Vec<f64> = serde_json::from_value(fixture["bounds"].clone()).unwrap();
                let quadratic: Vec<f64> =
                    serde_json::from_value(fixture["quadratic"].clone()).unwrap();
                let desired: Vec<f64> = serde_json::from_value(fixture["desired"].clone()).unwrap();
                let (values, certificate) =
                    ProjectedComponent::solve(&matrix, &bounds, &quadratic, &desired)
                        .unwrap_or_else(|error| panic!("{name}: {error}"));
                assert!(values.iter().all(|v| v.is_finite()), "{name}");
                assert!(
                    certificate.status.contains("polished"),
                    "{name}: {}",
                    certificate.status
                );
                assert!(certificate.primal_residual <= 2e-9, "{name}");
                assert!(certificate.stationarity_residual <= 2e-8, "{name}");
                assert!(certificate.complementarity_residual <= 2e-8, "{name}");
            } else {
                let proposal: Vec<[f64; 2]> =
                    serde_json::from_value(fixture["proposal"].clone()).unwrap();
                let indices: Vec<[i64; 2]> =
                    serde_json::from_value(fixture["indices"].clone()).unwrap();
                let weights: Vec<f64> = serde_json::from_value(fixture["weights"].clone()).unwrap();
                let points: Vec<usize> = serde_json::from_value(fixture["points"].clone()).unwrap();
                let normals: Vec<[f64; 2]> =
                    serde_json::from_value(fixture["normals"].clone()).unwrap();
                let bounds: Vec<f64> = serde_json::from_value(fixture["bounds"].clone()).unwrap();
                let scale = fixture["scale"].as_f64().unwrap();
                let mut dofs =
                    MovementDofs::new(&records(&vec![[0.0; 2]; proposal.len()]), false).unwrap();
                dofs.index = indices
                    .iter()
                    .map(|axes| axes.map(|x| usize::try_from(x).ok()))
                    .collect();
                dofs.weights = weights;
                let rows: Vec<_> = points
                    .into_iter()
                    .zip(normals)
                    .zip(bounds)
                    .map(|((point, normal), bound)| Halfplane {
                        point,
                        normal,
                        bound,
                    })
                    .collect();
                let (step, certificate) = dofs
                    .project(&proposal, &rows, scale)
                    .unwrap_or_else(|error| panic!("{name}: {error}"));
                assert!(step.iter().flatten().all(|v| v.is_finite()), "{name}");
                assert!(certificate.primal_residual <= 2e-9 * scale, "{name}");
                assert!(certificate.stationarity_residual <= 2e-8, "{name}");
                assert!(certificate.complementarity_residual <= 2e-8, "{name}");
                assert!(certificate.components > 1, "{name}");
            }
        }
    }
}
