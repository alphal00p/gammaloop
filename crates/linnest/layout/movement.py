"""Axis-wise native movement coordinates and certified sparse convex projection.

Clarabel and SciPy are pinned in requirements-movement.txt under scratch deps. The
objective gives each drawing coordinate unit weight: a grouped coordinate has
weight equal to its member count and receives their mean raw force. Route knot
charges affect force construction, not this Euclidean movement metric.
"""

from __future__ import annotations

import warnings

import clarabel
import numpy as np
from scipy import sparse
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve, qr
from scipy.sparse.csgraph import connected_components


def _pair(value):
    if value is None:
        return np.zeros(2)
    if isinstance(value, dict):
        if set(value) != {"x", "y"}:
            raise ValueError("Native shifts must contain x and y")
        value = [value["x"], value["y"]]
    result = np.asarray(value, dtype=float)
    if result.shape != (2,) or not np.isfinite(result).all():
        raise ValueError("Native coordinates must be finite pairs")
    return result


class MovementDOFs:
    """Eliminate affine native constraints without moving their reference pins."""

    def __init__(self, keys, records=None, *, nodes_fixed=False):
        self.keys = tuple(keys)
        count = len(keys)
        ids = {key: i for i, key in enumerate(keys)}
        parent = list(range(2 * count))

        def root(index):
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        def join(first, second):
            parent[root(first)] = root(second)

        fixed = np.zeros((count, 2), dtype=bool)
        self.reference = np.zeros((count, 2))
        self.relations = []
        self.signs = []
        if records is not None:
            records = list(records)
            if len(records) != count or {row["key"] for row in records} != set(keys):
                raise ValueError(
                    "Native constraint records must own each position exactly once"
                )
            by_id = {}
            for row in records:
                if row["kind"] in ("node", "edge"):
                    identity = (row["kind"], row["id"])
                    if identity in by_id:
                        raise ValueError(f"Duplicate native identity: {identity}")
                    by_id[identity] = row
            for row in records:
                index = ids[row["key"]]
                reference = np.asarray(row["initial_position"], dtype=float)
                shift = _pair(row.get("shift"))
                if (
                    reference.shape != (2,)
                    or shift.shape != (2,)
                    or not np.isfinite([reference, shift]).all()
                ):
                    raise ValueError(
                        "Native reference positions and shifts must be finite pairs"
                    )
                self.reference[index] = reference
                fixed[index] = row["fixed_axes"]
                for axis, name in enumerate(("x", "y")):
                    constraint = row["constraints"][name]
                    if constraint == "Free" or constraint == "Fixed":
                        if constraint == "Fixed":
                            fixed[index, axis] = True
                        continue
                    if not isinstance(constraint, dict) or set(constraint) != {
                        "Grouped"
                    }:
                        raise ValueError(
                            f"Unknown native axis constraint: {constraint}"
                        )
                    target, direction = constraint["Grouped"]
                    if len(target) != 1 or direction not in (
                        "Any",
                        "PositiveOnly",
                        "NegativeOnly",
                    ):
                        raise ValueError(f"Malformed native group: {constraint}")
                    kind, number = next(iter(target.items()))
                    other = by_id.get((kind.lower(), number))
                    if other is None:
                        raise ValueError(f"Unknown native group reference: {target}")
                    # Native layout-nodes=fixed ignores node group declarations.
                    if nodes_fixed and row["kind"] == "node":
                        continue
                    other_index = ids[other["key"]]
                    other_shift = _pair(other.get("shift"))
                    join(2 * index + axis, 2 * other_index + axis)
                    self.relations.append(
                        (index, other_index, axis, shift[axis] - other_shift[axis])
                    )
                    # Native fixed-node references also suppress the edge sign reflection.
                    if direction != "Any" and not (
                        nodes_fixed and other["kind"] == "node"
                    ):
                        sign = 1.0 if direction == "PositiveOnly" else -1.0
                        initial = sign * (reference[axis] - shift[axis])
                        if initial <= 0:
                            raise ValueError(
                                f"Native directed reference has the wrong sign: {row['key']}.{name}"
                            )
                        floor = min(
                            initial * 0.5,
                            1e-9 * max(1.0, abs(reference[axis]), abs(shift[axis])),
                        )
                        self.signs.append((index, axis, sign, shift[axis], floor))
        frozen = {root(index) for index in np.flatnonzero(fixed.ravel())}
        components = {}
        self.index = np.full((count, 2), -1, dtype=int)
        for coordinate in range(2 * count):
            representative = root(coordinate)
            if representative not in frozen:
                dof = components.setdefault(representative, len(components))
                self.index.ravel()[coordinate] = dof
        self.count = len(components)
        self.fixed = fixed
        self.weights = np.bincount(
            self.index[self.index >= 0], minlength=self.count
        ).astype(float)
        self.has_records = records is not None
        # Cycles are valid affine equivalences, provided their shifted references
        # agree; contradictory cycles and directed signs fail here, before a run.
        if self.has_records:
            self.validate(self.reference)

    def validate(self, p, tolerance=1e-9):
        if p.shape != self.index.shape or not np.isfinite(p).all():
            raise ValueError(
                "Movement coordinates must be finite and match native keys"
            )
        if self.has_records:
            error = np.max(
                np.abs(p[self.fixed] - self.reference[self.fixed]), initial=0.0
            )
            if error > tolerance:
                raise ArithmeticError(f"Native fixed reference changed by {error}")
            for first, second, axis, offset in self.relations:
                if abs(p[first, axis] - p[second, axis] - offset) > tolerance:
                    raise ArithmeticError(
                        "Native grouped coordinates disagree with their shifts"
                    )
            for index, axis, sign, shift, _floor in self.signs:
                if sign * (p[index, axis] - shift) <= 0:
                    raise ArithmeticError(
                        "Native directed coordinate crossed its sign boundary"
                    )

    def average(self, proposal):
        values = np.bincount(
            self.index[self.index >= 0],
            weights=proposal[self.index >= 0],
            minlength=self.count,
        )
        return values / self.weights

    def expand(self, values):
        result = np.zeros(self.index.shape)
        free = self.index >= 0
        result[free] = values[self.index[free]]
        return result

    def sign_rows(self, p):
        return [
            (
                index,
                [(-sign if axis == 0 else 0.0), (-sign if axis == 1 else 0.0)],
                sign * (p[index, axis] - shift) - floor,
            )
            for index, axis, sign, shift, floor in self.signs
        ]


def project_halfplanes(proposal, dofs, point_ids, normals, bounds, *, scale):
    """Nearest feasible displacement, with explicit primal and KKT certificates.

    The input rows mean normal dot movement[point_id] <= bound. Zero must be
    feasible. Equality elimination is exact; Clarabel handles only inequalities.
    A failed solve or residual check raises, and never substitutes zero motion.
    """
    if not np.isfinite(proposal).all() or not np.isfinite(bounds).all():
        raise ValueError("Nonfinite force or movement bounds")
    if np.min(bounds, initial=0.0) < -1e-12 * scale:
        raise ArithmeticError("Old drawing is outside a movement halfplane")
    bounds = np.maximum(bounds, 0.0) / scale
    columns = dofs.index[point_ids]
    rows = np.repeat(np.arange(len(bounds)), 2)
    values = normals.ravel()
    free = columns.ravel() >= 0
    matrix = sparse.coo_matrix(
        (values[free], (rows[free], columns.ravel()[free])),
        shape=(len(bounds), dofs.count),
    ).tocsc()
    matrix.eliminate_zeros()
    active = np.asarray(matrix.getnnz(axis=1) > 0).ravel()
    matrix, bounds = matrix[active].tocsc(), bounds[active]
    if not dofs.count:
        return np.zeros_like(proposal), {
            "status": "all coordinates fixed",
            "primal_residual": 0.0,
            "stationarity_residual": 0.0,
            "complementarity_residual": 0.0,
            "iterations": 0,
            "rows": 0,
        }
    desired = dofs.average(proposal) / scale
    # A row joins the X/Y DOFs of one point. Solve the connected components of
    # this incidence graph independently: unrelated forces and narrow corridors
    # must not determine each other's objective scaling or KKT conditioning.
    pattern = matrix.copy()
    pattern.data[:] = 1.0
    component_count, labels = connected_components(pattern.T @ pattern, directed=False)
    values = np.zeros(dofs.count)
    evidence = {
        "status": "independent convex movement components",
        "primal_residual": 0.0,
        "stationarity_residual": 0.0,
        "complementarity_residual": 0.0,
        "iterations": 0,
        "rows": len(bounds),
        "components": component_count,
    }
    for component in range(component_count):
        columns = np.flatnonzero(labels == component)
        block = matrix[:, columns]
        used = np.asarray(block.getnnz(axis=1) > 0).ravel()
        block = block[used].tocsc()
        candidate, certificate = _project_component(
            block, bounds[used], dofs.weights[columns], desired[columns]
        )
        values[columns] = candidate
        for name in (
            "primal_residual",
            "stationarity_residual",
            "complementarity_residual",
        ):
            evidence[name] = max(evidence[name], certificate[name])
        evidence["iterations"] += certificate["iterations"]
    evidence["primal_residual"] *= scale
    return dofs.expand(values * scale), evidence


def _project_small_polygon(matrix, bounds, quadratic, linear):
    """1D/2D convex projection by enumerating boundary feet and vertices.

    The closest point of a convex polygon is the desired point, a perpendicular
    boundary foot, or a polygon vertex. Enumerating these candidates is bounded
    by the local separator count and avoids ill-conditioned KKT solves in thin
    wedges. Extended precision retains opposing nearly parallel constraints.
    The original-coordinate certificate below still controls acceptance.
    """
    matrix = np.asarray(matrix.toarray(), dtype=np.longdouble)
    bounds = np.asarray(bounds, dtype=np.longdouble)
    quadratic = np.asarray(quadratic, dtype=np.longdouble)
    linear = np.asarray(linear, dtype=np.longdouble)
    desired = -linear / quadratic
    dual = np.zeros(len(bounds), dtype=np.longdouble)
    if len(desired) == 1:
        coefficient = matrix[:, 0]
        positive, negative = (
            np.flatnonzero(coefficient > 0),
            np.flatnonzero(coefficient < 0),
        )
        upper_row = positive[np.argmin(bounds[positive] / coefficient[positive])]
        lower_row = negative[np.argmax(bounds[negative] / coefficient[negative])]
        lower, upper = (
            bounds[lower_row] / coefficient[lower_row],
            bounds[upper_row] / coefficient[upper_row],
        )
        value = np.clip(desired[0], lower, upper)
        active = lower_row if desired[0] < lower else upper_row
        if value != desired[0]:
            dual[active] = -(quadratic[0] * value + linear[0]) / coefficient[active]
        return np.array([value], dtype=float), dual
    root_weights = np.sqrt(quadratic)
    transformed = matrix / root_weights
    lengths = np.sqrt(np.sum(transformed * transformed, axis=1))
    normals = transformed / lengths[:, None]
    offsets = bounds / lengths
    target = desired * root_weights
    distance = normals @ target - offsets
    target_roundoff = (
        128 * np.finfo(np.longdouble).eps * max(1.0, np.max(np.abs(target)))
    )
    if np.all(distance <= target_roundoff):
        return np.asarray(desired, dtype=float), dual
    feet = target - distance[:, None] * normals
    first, second = np.triu_indices(len(bounds), 1)
    det = (
        normals[first, 0] * normals[second, 1] - normals[first, 1] * normals[second, 0]
    )
    nonparallel = det != 0
    first, second, det = first[nonparallel], second[nonparallel], det[nonparallel]
    vertices = np.column_stack(
        (
            (offsets[first] * normals[second, 1] - normals[first, 1] * offsets[second])
            / det,
            (normals[first, 0] * offsets[second] - offsets[first] * normals[second, 0])
            / det,
        )
    )
    candidates = np.concatenate((feet, vertices))
    violation = normals @ candidates.T - offsets[:, None]
    roundoff = (
        128
        * np.finfo(np.longdouble).eps
        * np.maximum(1.0, np.max(np.abs(candidates), axis=1))
    )
    feasible = np.all(violation <= roundoff[None, :], axis=0)
    if not np.any(feasible):
        raise ArithmeticError(
            "No certified boundary candidate in a feasible movement polygon"
        )
    # Filter by the normal cone as well as primal feasibility. Near-coincident
    # vertices can have indistinguishable objective values, but a negative
    # multiplier identifies the wrong corner regardless of that roundoff tie.
    delta = target - vertices
    first_dual = (
        (delta[:, 0] * normals[second, 1] - normals[second, 0] * delta[:, 1])
        / det
        / lengths[first]
    )
    second_dual = (
        (normals[first, 0] * delta[:, 1] - delta[:, 0] * normals[first, 1])
        / det
        / lengths[second]
    )
    nonnegative_dual = np.concatenate(
        (distance / lengths >= -2e-9, (first_dual >= -2e-9) & (second_dual >= -2e-9))
    )
    feasible &= nonnegative_dual
    if not np.any(feasible):
        raise ArithmeticError(
            "No primal/dual certified candidate in a movement polygon"
        )
    # Drop the constant target-squared term to retain comparisons under strong
    # raw forces; the movement cap bounds candidates, not the desired point.
    cost = 0.5 * np.sum(candidates * candidates, axis=1) - candidates @ target
    cost[~feasible] = np.inf
    selected = int(np.argmin(cost))
    value = candidates[selected]
    if selected < len(bounds):
        dual[selected] = distance[selected] / lengths[selected]
    else:
        vertex = selected - len(bounds)
        dual[first[vertex]] = first_dual[vertex]
        dual[second[vertex]] = second_dual[vertex]
    return np.asarray(value / root_weights, dtype=float), dual


def _polish_active_face(matrix, bounds, quadratic, linear, values, dual):
    """Refine a conic candidate's active face without accepting a failed solve.

    Near-parallel active constraints can amplify tiny conic slack errors. Factor
    their small KKT system once, then refine its residual in extended precision.
    A negative multiplier excludes that face. Attempts only remove predicted
    active constraints and are bounded by their original count; the caller must
    still certify the result against every original inequality.
    """
    if (
        values.shape != quadratic.shape
        or dual.shape != bounds.shape
        or not np.isfinite(values).all()
        or not np.isfinite(dual).all()
    ):
        return None
    dense = matrix.toarray()
    active = np.flatnonzero((dual > 1e-7) & (np.abs(dense @ values - bounds) < 1e-7))
    for _attempt in range(len(active) + 1):
        if len(active):
            _, triangular, pivot = qr(dense[active].T, mode="economic", pivoting=True)
            rank = np.count_nonzero(
                np.abs(np.diag(triangular))
                > 1e-12 * max(1.0, np.max(np.abs(triangular)))
            )
            basis = active[pivot[:rank]]
        else:
            rank, basis = 0, active
        constraints = dense[basis]
        system = np.block(
            [[np.diag(quadratic), constraints.T], [constraints, np.zeros((rank, rank))]]
        )
        right = np.concatenate((-linear, bounds[basis]))
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", LinAlgWarning)
                factor = lu_factor(system)
                solution = lu_solve(factor, right).astype(np.longdouble)
                extended_system = system.astype(np.longdouble)
                extended_right = right.astype(np.longdouble)
                for _ in range(6):
                    residual = extended_right - extended_system @ solution
                    solution += lu_solve(
                        factor, np.asarray(residual, dtype=float)
                    ).astype(np.longdouble)
        except (ValueError, np.linalg.LinAlgError, LinAlgWarning):
            return None
        if not np.isfinite(solution).all():
            return None
        candidate = np.asarray(solution[: len(quadratic)], dtype=float)
        multipliers = np.zeros(len(bounds), dtype=np.longdouble)
        multipliers[basis] = solution[len(quadratic) :]
        if np.min(multipliers, initial=0.0) < -2e-9:
            active = active[active != np.argmin(multipliers)]
            continue
        return candidate, multipliers
    return None


def _project_component(matrix, bounds, quadratic, desired):
    """Solve one genuinely coupled block, preserving the unit-coordinate metric."""
    linear = -quadratic * desired
    objective_scale = max(1.0, float(np.max(np.abs(linear), initial=0.0)))
    quadratic, linear = quadratic / objective_scale, linear / objective_scale
    # Certify the unconstrained optimum directly when it already satisfies all
    # bounds. This is an exact fast path, not a numerical failure fallback.
    if np.all(matrix @ desired <= bounds):
        return desired.copy(), {
            "status": "unconstrained optimum",
            "primal_residual": 0.0,
            "stationarity_residual": 0.0,
            "complementarity_residual": 0.0,
            "iterations": 0,
            "rows": len(bounds),
        }

    def certify(values, dual):
        if values is None or dual is None:
            return (np.inf,) * 4
        values, dual = np.asarray(values), np.asarray(dual)
        if (
            values.shape != desired.shape
            or dual.shape != (len(bounds),)
            or not np.isfinite(values).all()
            or not np.isfinite(dual).all()
        ):
            return (np.inf,) * 4
        extended = matrix.astype(np.longdouble)
        values, dual = values.astype(np.longdouble), dual.astype(np.longdouble)
        excess = extended @ values - bounds.astype(np.longdouble)
        return (
            float(np.max(excess, initial=0.0)),
            float(
                np.max(
                    np.abs(quadratic * values + linear + extended.T @ dual), initial=0.0
                )
            ),
            float(np.max(np.abs(dual * excess), initial=0.0)),
            float(max(0.0, -np.min(dual, initial=0.0))),
        )

    def certified(evidence):
        return all(
            value <= limit for value, limit in zip(evidence, (2e-9, 2e-8, 2e-8, 2e-9))
        )

    if len(desired) <= 2:
        values, dual = _project_small_polygon(matrix, bounds, quadratic, linear)
        status, iterations = "enumerated movement polygon", 0
        evidence = certify(values, dual)
    else:
        settings = clarabel.DefaultSettings()
        settings.verbose = False
        settings.max_threads = 1
        settings.max_iter = 200
        settings.direct_solve_method = "qdldl"
        settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-12
        # The frozen K5 thin-region regression requires large opposing multipliers.
        # Default 1e-8 KKT regularization stalls above its feasible corridor. Use
        # dynamic pivot regularization and iterative refinement for the original QP;
        # neither its geometry bounds nor its acceptance tolerances are changed.
        settings.static_regularization_constant = 0.0
        settings.dynamic_regularization_eps = 1e-15
        settings.dynamic_regularization_delta = 1e-10
        settings.iterative_refinement_abstol = 1e-15
        settings.iterative_refinement_reltol = 1e-15
        settings.iterative_refinement_max_iter = 30
        solution = clarabel.DefaultSolver(
            sparse.diags(quadratic, format="csc"),
            linear,
            matrix,
            bounds,
            [clarabel.NonnegativeConeT(len(bounds))],
            settings,
        ).solve()
        values, dual = np.asarray(solution.x), np.asarray(solution.z)
        status = f"Clarabel {solution.status}"
        iterations = int(solution.iterations)
        evidence = certify(values, dual)
    if len(desired) > 2 and not certified(evidence):
        polished = _polish_active_face(matrix, bounds, quadratic, linear, values, dual)
        if polished is not None:
            candidate, multipliers = polished
            candidate_evidence = certify(candidate, multipliers)
            if certified(candidate_evidence):
                values, dual, evidence = candidate, multipliers, candidate_evidence
                status += " (active-face polished)"
    primal, stationarity, complementarity, dual_violation = evidence
    if not certified(evidence):
        raise ArithmeticError(
            f"Uncertified joint movement ({status}): primal={primal}, stationarity={stationarity}, complementarity={complementarity}, dual_sign={dual_violation}"
        )
    return values, {
        "status": status,
        "primal_residual": primal,
        "stationarity_residual": stationarity,
        "complementarity_residual": complementarity,
        "iterations": iterations,
        "rows": len(bounds),
    }
