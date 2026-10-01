"""Experimental ImPrEd-family planar relaxation; no production API/default changes.

Uses paper force-exponent schedules and midpoint separating halfplanes (the
less restrictive convex variant), extended to finite segments. Each subdivided
segment applies ImPrEd attraction. Route samples share fixed charge budgets,
and flexible external routes do not repel themselves. These repulsion
extensions, balanced external flow, and the conservative movement schedule
differ from the paper's original force model.
All moved points share old-state separators; route changes are independently checked.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import time
from collections import Counter
from concurrent.futures import CancelledError
from itertools import combinations, pairwise
from pathlib import Path

import numpy as np
from geometry import validate_geometry
from movement import MovementDOFs, project_halfplanes
from seed import initialize


def model_arrays(state, *, flexible_edges=None, parallel_attraction_balance=0.0):
    """Build fixed physical-edge charge budgets and flexible repulsion masks.

    External tips and knots share one unit charge; internal knots share 0.35.
    Flexible dangling routes exclude all pairs within their own connection so
    inserting a knot cannot reintroduce owner/tip repulsion. Internal routes
    retain ImPrEd's consecutive-only exemption, including self-loops.
    Flexibility is an explicit edge label, independent of subdivision count.
    """
    flexible_edges = set(
        state.get("flexible_edges", ()) if flexible_edges is None else flexible_edges
    )
    keys = list(state["positions"])
    ids = {k: i for i, k in enumerate(keys)}
    p = np.array([state["positions"][k] for k in keys], dtype=float)
    routes = [[ids[k] for k in route] for route in state["routes"].values()]
    segments = np.array(
        [(a, b) for route in routes for a, b in pairwise(route)], dtype=int
    )
    charges = np.array([0.35 if k.startswith("b:") else 1.0 for k in keys])
    for edge, route in state["routes"].items():
        endpoint = state.get("external_ids", {}).get(edge)
        samples = [key for key in route if key.startswith("b:") or key == endpoint]
        if samples:
            total = 1.0 if endpoint is not None else 0.35
            charges[[ids[key] for key in samples]] = total / len(samples)
    weights = charges[:, None] * charges[None, :]
    np.fill_diagonal(weights, 0)
    for edge, route in zip(state["routes"], routes):
        if edge in flexible_edges:
            pairs = (
                combinations(route, 2)
                if edge in state.get("external_ids", {})
                else pairwise(route)
            )
            for first, second in pairs:
                weights[first, second] = weights[second, first] = 0
    # Physical multiplicity ignores bends, direction and external connections.
    endpoint_pairs = {
        edge: tuple(sorted((route[0], route[-1])))
        for edge, route in state["routes"].items()
        if edge not in state.get("external_ids", {}) and route[0] != route[-1]
    }
    multiplicities = Counter(endpoint_pairs.values())
    attraction_weights = np.array(
        [
            multiplicities[endpoint_pairs[edge]] ** (-parallel_attraction_balance)
            if edge in endpoint_pairs
            else 1.0
            for edge, route in state["routes"].items()
            for _ in pairwise(route)
        ]
    )
    return keys, ids, p, routes, segments, weights, attraction_weights


def point_segments(p, segments):
    a, b = p[segments[:, 0]], p[segments[:, 1]]
    ab = b - a
    t = np.einsum("nsi,si->ns", p[:, None, :] - a, ab) / np.maximum(
        np.sum(ab * ab, axis=1), 1e-20
    )
    t = np.clip(t, 0, 1)
    diff = p[:, None, :] - (a + t[:, :, None] * ab)
    d = np.linalg.norm(diff, axis=2)
    normal = diff / np.maximum(d[:, :, None], 1e-15)
    mask = (np.arange(len(p))[:, None] != segments[None, :, 0]) & (
        np.arange(len(p))[:, None] != segments[None, :, 1]
    )
    return t, d, normal, mask


def node_edge_forces(segments, data, *, target, clearance, strength, exponent):
    """Repel from open segment interiors, retaining balanced endpoint reactions.

    ImPrEd's force applies only to perpendicular projections inside the edge.
    The same finite-segment data keeps endpoint distances for topology guards.
    """
    t, dist, normal, valid = data
    gamma = target * clearance
    magnitude = np.where(
        valid & (t > 0) & (t < 1) & (dist < gamma),
        strength * target * (1 - dist / gamma).clip(0) ** exponent,
        0,
    )
    pairs = normal * magnitude[:, :, None]
    force = np.sum(pairs, axis=1)
    np.add.at(force, segments[:, 0], -np.sum(pairs * (1 - t)[:, :, None], axis=0))
    np.add.at(force, segments[:, 1], -np.sum(pairs * t[:, :, None], axis=0))
    return force


def segment_attraction_forces(p, segments, *, target, strength, exponent):
    """Apply ImPrEd attraction independently to each subdivided segment.

    Tension A l (l/δ)^a gives longer segments more pull. Splitting adds serial
    compliance; unequal adjacent tensions move a knot toward equal spacing.
    Strength is scalar or one coefficient per segment, without changing this law.
    """
    force = np.zeros_like(p)
    vectors = p[segments[:, 1]] - p[segments[:, 0]]
    lengths = np.linalg.norm(vectors, axis=1)
    values = strength * (lengths / target) ** exponent
    attraction = values[:, None] * vectors
    np.add.at(force, segments[:, 0], attraction)
    np.add.at(force, segments[:, 1], -attraction)
    return force


def separating_halfplanes(p, segments, pair_data=None):
    """Return old-state point/segment and point/point movement inequalities.

    Each row is normal dot displacement[point] <= available room. The stable
    perpendicular and clearance are shared by both sides of each separator.
    These constraints apply to the entire simultaneous straight-line step.
    """
    point_rows, normal_rows, bound_rows = [], [], []
    if len(segments):
        t, d, normal, valid = (
            pair_data if pair_data is not None else point_segments(p, segments)
        )
        a, b = p[segments[:, 0]], p[segments[:, 1]]
        vectors = b - a
        lengths = np.linalg.norm(vectors, axis=1)
        perpendicular = np.column_stack((-vectors[:, 1], vectors[:, 0]))
        perpendicular /= np.maximum(lengths[:, None], 1e-30)
        signed = np.einsum("nsi,si->ns", p[:, None, :] - a, perpendicular)
        interior = (t > 0) & (t < 1)
        stable = (
            perpendicular[None, :, :] * np.where(signed >= 0, 1.0, -1.0)[:, :, None]
        )
        normal = np.where(interior[:, :, None], stable, normal)
        d = np.where(interior, np.abs(signed), d)
        clearance = 1e-6 * max(float(np.median(lengths)), 1.0)
        room = 0.45 * np.maximum(d - clearance, 0)
        point_ids, segment_ids = np.nonzero(valid)
        point_rows.append(point_ids)
        normal_rows.append(-normal[valid])
        bound_rows.append(room[valid])
        for endpoint in (0, 1):
            vertex_ids = segments[:, endpoint]
            offset = np.einsum(
                "nsi,nsi->ns", normal, p[:, None, :] - p[vertex_ids][None, :, :]
            )
            # Interior offsets equal d exactly; avoid subtractive cancellation.
            available = np.where(interior, room, np.maximum(offset - d, 0) + room)
            point_rows.append(vertex_ids[segment_ids])
            normal_rows.append(normal[valid])
            bound_rows.append(available[valid])
    else:
        clearance = 1e-6
    differences = p[:, None, :] - p[None, :, :]
    distances = np.linalg.norm(differences, axis=2)
    normals = differences / np.maximum(distances[:, :, None], 1e-30)
    valid = ~np.eye(len(p), dtype=bool)
    point_rows.append(np.nonzero(valid)[0])
    normal_rows.append(-normals[valid])
    bound_rows.append((0.45 * np.maximum(distances - clearance, 0))[valid])
    return (
        np.concatenate(point_rows),
        np.concatenate(normal_rows),
        np.concatenate(bound_rows),
    )


def protect(
    p, segments, proposal, pair_data=None, *, dofs=None, cap=None, certificate=None
):
    """Project raw movement jointly onto native constraints and safe halfplanes.

    Grouped axes share a displacement DOF; the other axis remains independent.
    The objective aggregates raw forces before the per-point trust region is
    applied. A regular 16-gon inscribed in the Euclidean movement cap makes the
    solve a convex QP (maximum directional conservatism 1-cos(pi/16), 1.92%).
    The returned fractions are diagnostics, not per-point scalar multipliers.
    """
    if not len(p):
        return proposal.copy(), np.ones(0)
    dofs = MovementDOFs(range(len(p))) if dofs is None else dofs
    dofs.validate(p)
    if cap is None:
        cap = max(float(np.linalg.norm(proposal, axis=1).max()), 1e-12)
    if not math.isfinite(cap) or cap <= 0:
        raise ValueError("Movement cap must be finite and positive")
    point_ids, normals, bounds = separating_halfplanes(p, segments, pair_data)
    # Unit-normal rows outside the trust disk cannot bind. Keep the full arrays
    # for the independent certificate after solving the reduced sparse problem.
    near = bounds < cap * (1 + 1e-12)
    solve_ids, solve_normals, solve_bounds = (
        point_ids[near],
        normals[near],
        bounds[near],
    )
    angles = (np.arange(16) + 0.5) * (2 * np.pi / 16)
    directions = np.column_stack((np.cos(angles), np.sin(angles)))
    solve_ids = np.concatenate((solve_ids, np.repeat(np.arange(len(p)), 16)))
    solve_normals = np.concatenate((solve_normals, np.tile(directions, (len(p), 1))))
    solve_bounds = np.concatenate(
        (solve_bounds, np.full(len(p) * 16, cap * math.cos(math.pi / 16)))
    )
    signs = dofs.sign_rows(p)
    if signs:
        solve_ids = np.concatenate((solve_ids, [row[0] for row in signs]))
        solve_normals = np.concatenate((solve_normals, [row[1] for row in signs]))
        solve_bounds = np.concatenate((solve_bounds, [row[2] for row in signs]))
    step, evidence = project_halfplanes(
        proposal, dofs, solve_ids, solve_normals, solve_bounds, scale=cap
    )
    separator_error = float(
        np.max(np.einsum("ij,ij->i", normals, step[point_ids]) - bounds, initial=0.0)
    )
    cap_error = max(0.0, float(np.linalg.norm(step, axis=1).max()) - cap)
    if separator_error > 2e-9 * cap or cap_error > 2e-9 * cap:
        raise ArithmeticError(
            f"Joint movement geometry certificate failed: separators={separator_error}, cap={cap_error}"
        )
    dofs.validate(p + step)
    if certificate is not None:
        certificate.update(
            evidence,
            separator_residual=separator_error,
            cap_residual=cap_error,
            dofs=dofs.count,
        )
    norms = np.linalg.norm(proposal, axis=1)
    fractions = np.minimum(1.0, np.linalg.norm(step, axis=1) / np.maximum(norms, 1e-30))
    fractions[norms == 0] = 1.0
    return step, fractions


def adapt(
    state,
    target,
    serial,
    *,
    edges=None,
    max_points=None,
    remove=True,
    protected_points=(),
    split_length_ratio=3.0,
    contract_chord_ratio=2.0,
):
    """ImPrEd §5.3 refinement, with default split 3δ and contraction chord 2δ.

    ``max_points`` is an optional user cap, not a length or angle criterion.
    Physical anchors, constrained knots, and the existing crossing witness stay
    unchanged even when the paper's local geometric conditions permit removal.
    """
    if not (
        math.isfinite(split_length_ratio)
        and math.isfinite(contract_chord_ratio)
        and 0 < contract_chord_ratio < split_length_ratio
    ):
        raise ValueError("Refinement ratios must be finite with 0 < contract < split")
    positions = state["positions"]
    protected = set(protected_points) | set(state.get("anchor_ids", {}).values())
    initial_crossings = validate_geometry(positions, state["routes"])["crossings"]
    added = removed = 0
    for edge, route in list(state["routes"].items()):
        if edges is not None and edge not in edges:
            continue
        expanded = [route[0]]
        room = math.inf if max_points is None else max_points - (len(route) - 2)
        for a, b in pairwise(route):
            if (
                math.dist(positions[a], positions[b]) > target * split_length_ratio
                and room > 0
            ):
                key = f"b:{edge}:auto{serial}"
                serial += 1
                while key in positions:
                    key = f"b:{edge}:auto{serial}"
                    serial += 1
                positions[key] = [
                    (positions[a][k] + positions[b][k]) / 2 for k in (0, 1)
                ]
                expanded.append(key)
                added += 1
                room -= 1
            expanded.append(b)
        state["routes"][edge] = expanded
    if not remove:
        return serial, added, removed
    for edge in state["routes"]:
        if edges is not None and edge not in edges:
            continue
        route = state["routes"][edge]
        j = 1
        while j < len(route) - 1:
            corner_ids = route[j - 1 : j + 2]
            a, b, c = [np.array(positions[k]) for k in corner_ids]
            chord = float(np.linalg.norm(c - a))
            # Separate insertion/removal thresholds prevent immediate churn.
            removable = (
                route[j] not in protected
                and chord > 1e-7
                and chord < contract_chord_ratio * target
            )
            if removable:
                # Endpoint crossing checks alone miss jumping an isolated vertex
                # or whole component across the swept corner triangle. Require
                # its interior and boundary to contain no other drawing point.
                low = np.minimum(np.minimum(a, b), c)
                high = np.maximum(np.maximum(a, b), c)
                tolerance = 1e-10 * max(target, 1.0)
                for point_id, coordinates in positions.items():
                    if point_id in corner_ids:
                        continue
                    point = np.array(coordinates)
                    if np.any(point < low - tolerance) or np.any(
                        point > high + tolerance
                    ):
                        continue
                    orientations = []
                    for first, second in ((a, b), (b, c), (c, a)):
                        edge_vector, delta = second - first, point - first
                        orientations.append(
                            edge_vector[0] * delta[1] - edge_vector[1] * delta[0]
                        )
                    if (
                        min(orientations) >= -tolerance
                        or max(orientations) <= tolerance
                    ):
                        removable = False
                        break
            if removable:
                key = route.pop(j)
                # The deleted knot must not remain as an isolated obstacle in
                # the replacement geometry's contact/crossing audit.
                candidate_positions = positions.copy()
                candidate_positions.pop(key)
                geometry = validate_geometry(candidate_positions, state["routes"])
                if (
                    not any(geometry[k] for k in ("contacts", "overlaps", "degenerate"))
                    and geometry["crossings"] == initial_crossings
                ):
                    positions.pop(key)
                    removed += 1
                    continue
                route.insert(j, key)
            j += 1
    return serial, added, removed


def solve(
    initial,
    steps=1000,
    target=2.4,
    flexible=True,
    pull=0.45,
    edge_clearance=0.4,
    node_edge_strength=4.0,
    project=None,
    repulsion=2.5,
    attraction=2.5,
    parallel_attraction_balance=0.0,
    pull_scales=None,
    pull_balance=1.0,
    external_max_points=0,
    split_length_ratio=3.0,
    contract_chord_ratio=2.0,
    cancelled=None,
    progress=None,
):
    if not all(
        math.isfinite(value) and value >= 0
        for value in (
            repulsion,
            attraction,
            pull,
            pull_balance,
            parallel_attraction_balance,
        )
    ):
        raise ValueError("Force multipliers must be finite and nonnegative")
    if not (
        math.isfinite(split_length_ratio)
        and math.isfinite(contract_chord_ratio)
        and 0 < contract_chord_ratio < split_length_ratio
    ):
        raise ValueError("Refinement ratios must be finite with 0 < contract < split")
    if external_max_points not in (0, 1, 2, 3):
        raise ValueError("External routes support at most three intermediate points")
    if external_max_points and (project is None or flexible):
        raise ValueError(
            "External refinement requires the native projector and fixed internal routes"
        )
    if cancelled is not None and cancelled():
        raise CancelledError()
    if not initial["positions"]:
        return copy.deepcopy(initial)
    if project is not None and flexible:
        raise ValueError(
            "Native constraint sessions currently require fixed route identities"
        )
    started = time.monotonic()
    state = copy.deepcopy(initial)
    initial_geometry = validate_geometry(state["positions"], state["routes"])
    if any(initial_geometry[k] for k in ("contacts", "overlaps", "degenerate")):
        raise ValueError(f"Invalid initial carrier geometry: {initial_geometry}")
    initial_crossings = initial_geometry["crossings"]
    crossed = {part[0] for crossing in initial_crossings for part in crossing}
    external_edges = set(state.get("external_ids", {}))
    # Internal native point identities may be fixed while bends remain movable.
    # Choose labels once; contracting a flexible carrier never makes it rigid.
    seed_flexible = {
        edge
        for edge, route in state["routes"].items()
        if edge not in external_edges and len(route) > 2
    }
    flexible_edges = set(state.get("flexible_edges", seed_flexible))
    if flexible:
        flexible_edges.update(state["routes"])
    if external_max_points:
        flexible_edges.update(external_edges)
    flexible_edges -= crossed
    if not flexible_edges <= set(state["routes"]):
        raise ValueError("Flexible edge identities must name existing routes")
    state["flexible_edges"] = sorted(flexible_edges)
    added = removed = serial = 0
    histories = []
    checked = 0
    min_scale = 1.0
    max_constraint_residual = 0.0
    keys, ids, p, _routes, segs, weights, attraction_weights = model_arrays(
        state,
        flexible_edges=flexible_edges,
        parallel_attraction_balance=parallel_attraction_balance,
    )
    dofs = project.movement_constraints(keys) if project is not None else None
    movement_certificates = {
        "method": "joint convex QP in native axis DOFs",
        "metric": "unit weight per drawing coordinate; grouped raw forces averaged before trust bounds",
        "trust_region": "16-gon inscribed in the scheduled Euclidean movement cap",
        "acceptance": "finite primal, dual-sign, stationarity and complementarity certificates; no zero-step failure fallback",
        "normalized_primal_tolerance": 2e-9,
        "normalized_stationarity_tolerance": 2e-8,
        "normalized_complementarity_tolerance": 2e-8,
        "coordinate_tolerance": "normalized primal tolerance times scheduled movement cap",
        "solves": 0,
        "status_counts": {},
        "max_iterations": 0,
        "max_primal_residual": 0.0,
        "max_stationarity_residual": 0.0,
        "max_complementarity_residual": 0.0,
        "max_separator_residual": 0.0,
        "max_cap_residual": 0.0,
    }
    if project is not None:
        projected = project(state["positions"])
        residual = max(math.dist(projected[k], state["positions"][k]) for k in keys)
        if residual > 1e-9:
            raise ValueError(
                "Apply native constraints to the seed before topology protection"
            )
    incoming = set(state["incoming_ids"])
    outgoing = set(state["outgoing_ids"])
    if pull_scales is not None and (
        set(pull_scales) != incoming | outgoing
        or any(not math.isfinite(value) or value <= 0 for value in pull_scales.values())
    ):
        raise ValueError(
            "Native pull scales must cover all external endpoints with positive finite values"
        )
    # Native topology and coordinate groups stay fixed as routes refine. Evaluate
    # once, including the zero-pull case before exponentiation: group clearance
    # can produce scales above one, so a finite exponent can overflow.
    try:
        endpoint_pull = {
            key: 0.0
            if pull == 0
            else pull
            * target
            * (pull_scales[key] ** pull_balance if pull_scales is not None else 1.0)
            for key in incoming | outgoing
        }
    except OverflowError as error:
        raise ValueError(
            "Topology balancing produces a non-finite external force"
        ) from error
    if any(not math.isfinite(value) for value in endpoint_pull.values()):
        raise ValueError("Topology balancing produces a non-finite external force")
    physical_ids = [ids[key] for key in state["node_ids"].values()]
    refinement_events = []
    if not len(segs):
        return state
    for epoch in range(steps):
        if epoch % 25 == 0 and cancelled is not None and cancelled():
            raise CancelledError()
        u = epoch / max(steps - 1, 1)
        rep_exp = 2 + 2 * u
        att_exp = 1 - 0.6 * u
        d = p[:, None, :] - p[None, :, :]
        r = np.linalg.norm(d, axis=2)
        coefficients = weights * (target / np.maximum(r, target * 0.02)) ** rep_exp
        force = repulsion * np.sum(d * coefficients[:, :, None], axis=1)
        force += segment_attraction_forces(
            p,
            segs,
            target=target,
            strength=attraction * attraction_weights,
            exponent=att_exp,
        )
        data = point_segments(p, segs)
        force += node_edge_forces(
            segs,
            data,
            target=target,
            clearance=edge_clearance,
            strength=node_edge_strength,
            exponent=rep_exp,
        )
        reaction = np.zeros(2)
        if incoming and outgoing:
            for collection, sign in ((incoming, -1), (outgoing, 1)):
                for key in collection:
                    outward = sign * endpoint_pull[key]
                    force[ids[key], 0] += outward
                    reaction[0] += outward
        else:
            center = np.mean([p[ids[k]] for k in state["node_ids"].values()], axis=0)
            for key in incoming | outgoing:
                vector = p[ids[key]] - center
                length = np.linalg.norm(vector)
                if length > 1e-9:
                    outward = endpoint_pull[key] * vector / length
                    force[ids[key]] += outward
                    reaction += outward
        if pull_scales is not None and physical_ids:
            force[physical_ids] -= reaction / len(physical_ids)
        proposal = force * 0.045
        cap = target * (0.10 * (1 - u) + 0.002)
        evidence = {}
        step, scales = protect(
            p, segs, proposal, data, dofs=dofs, cap=cap, certificate=evidence
        )
        movement_certificates["solves"] += 1
        statuses = movement_certificates["status_counts"]
        statuses[evidence["status"]] = statuses.get(evidence["status"], 0) + 1
        movement_certificates["max_iterations"] = max(
            movement_certificates["max_iterations"], evidence["iterations"]
        )
        for field in (
            "primal_residual",
            "stationarity_residual",
            "complementarity_residual",
            "separator_residual",
            "cap_residual",
        ):
            name = f"max_{field}"
            movement_certificates[name] = max(
                movement_certificates[name], evidence[field]
            )
        min_scale = min(min_scale, float(scales.min()))
        p += step
        if project is None:
            p -= np.mean(p, axis=0)
        if epoch % 25 == 0 or epoch == steps - 1:
            state["positions"] = {k: list(map(float, p[i])) for i, k in enumerate(keys)}
            geometry = validate_geometry(state["positions"], state["routes"])
            checked += 1
            if any(geometry[k] for k in ("contacts", "overlaps", "degenerate")) or (
                geometry["crossings"] != initial_crossings
            ):
                raise ArithmeticError(
                    f"Topology changed at iteration {epoch}: {geometry}"
                )
            if epoch % 100 == 0 or epoch == steps - 1:
                if project is not None:
                    projected = project(state["positions"])
                    residual = max(
                        math.dist(projected[k], state["positions"][k]) for k in keys
                    )
                    max_constraint_residual = max(max_constraint_residual, residual)
                    if residual > 1e-9:
                        raise ArithmeticError(
                            f"Constraint changed at iteration {epoch}: {residual}"
                        )
                histories.append(
                    {
                        "iteration": epoch,
                        "positions": copy.deepcopy(state["positions"]),
                        "routes": copy.deepcopy(state["routes"]),
                        "max_move": float(np.linalg.norm(step, axis=1).max()),
                    }
                )
                if progress is not None:
                    progress(
                        {
                            **copy.deepcopy(histories[-1]),
                            "iteration": epoch + 1,
                            "total": steps,
                        }
                    )
            if flexible and epoch > 0:
                serial, a, r = adapt(
                    state,
                    target,
                    serial,
                    edges=flexible_edges,
                    split_length_ratio=split_length_ratio,
                    contract_chord_ratio=contract_chord_ratio,
                )
                added += a
                removed += r
                keys, ids, p, _routes, segs, weights, attraction_weights = model_arrays(
                    state,
                    flexible_edges=flexible_edges,
                    parallel_attraction_balance=parallel_attraction_balance,
                )
                physical_ids = [ids[key] for key in state["node_ids"].values()]
                if histories and histories[-1]["iteration"] == epoch:
                    histories[-1]["positions"] = copy.deepcopy(state["positions"])
                    histories[-1]["routes"] = copy.deepcopy(state["routes"])
            if (
                external_max_points
                and epoch > 0
                and (epoch % 25 == 0 or epoch == steps - 1)
            ):
                # Contract first, then register it before any insertion. A new
                # midpoint belongs to the resulting native route, not an old
                # segment whose unnecessary corner has just been removed.
                for operation in ("contract", "insert"):
                    if operation == "contract":
                        eligible = external_edges - crossed
                        protected = {
                            row["key"]
                            for row in project.ready["constraints"]
                            if any(row["fixed_axes"])
                            or any(
                                value != "Free" for value in row["constraints"].values()
                            )
                        }
                        before_routes = copy.deepcopy(state["routes"])
                        serial, inserted, contracted = adapt(
                            state,
                            target,
                            serial,
                            edges=eligible,
                            split_length_ratio=split_length_ratio,
                            contract_chord_ratio=contract_chord_ratio,
                            max_points=0,
                            protected_points=protected,
                        )
                    else:
                        eligible = external_edges - crossed
                        before_routes = copy.deepcopy(state["routes"])
                        serial, inserted, contracted = adapt(
                            state,
                            target,
                            serial,
                            edges=eligible,
                            split_length_ratio=split_length_ratio,
                            contract_chord_ratio=contract_chord_ratio,
                            max_points=external_max_points,
                            remove=False,
                        )
                    if inserted or contracted:
                        geometry = validate_geometry(
                            state["positions"], state["routes"]
                        )
                        if (
                            any(
                                geometry[key]
                                for key in ("contacts", "overlaps", "degenerate")
                            )
                            or geometry["crossings"] != initial_crossings
                        ):
                            raise ArithmeticError(
                                "External route update changed carrier topology"
                            )
                        project.update_external_routes(state)
                        added += inserted
                        removed += contracted
                        keys, ids, p, _routes, segs, weights, attraction_weights = (
                            model_arrays(
                                state,
                                flexible_edges=flexible_edges,
                                parallel_attraction_balance=parallel_attraction_balance,
                            )
                        )
                        physical_ids = [ids[key] for key in state["node_ids"].values()]
                        # Rebuild membership after deletion; original pins and
                        # surviving constrained coordinates stay native-owned.
                        dofs = project.movement_constraints(keys)
                        if histories and histories[-1]["iteration"] == epoch:
                            histories[-1]["positions"] = copy.deepcopy(
                                state["positions"]
                            )
                            histories[-1]["routes"] = copy.deepcopy(state["routes"])
                        refinement_events.append(
                            {
                                "iteration": epoch + 1,
                                "edges": sorted(
                                    e
                                    for e in state["routes"]
                                    if state["routes"][e] != before_routes[e]
                                ),
                                "inserted_points": inserted,
                                "removed_points": contracted,
                            }
                        )
    state["history"] = histories
    state["report"] = {
        **initial["report"],
        "method": "experimental ImPrEd-family: variable powers + convex midpoint halfplanes + open-segment forces",
        "steps": steps,
        "target_segment_length": target,
        "flexible_edges": sorted(flexible_edges),
        "split_length_ratio": split_length_ratio,
        "contract_chord_ratio": contract_chord_ratio,
        "flexible": flexible,
        "edge_clearance_ratio": edge_clearance,
        "pull": pull,
        "node_edge_strength": node_edge_strength,
        "repulsion": repulsion,
        "attraction": attraction,
        "parallel_attraction_balance": parallel_attraction_balance,
        "external_pull_scales": pull_scales,
        "external_pull_balance": pull_balance if pull_scales is not None else None,
        "external_max_points": external_max_points,
        "external_refinement_events": refinement_events,
        "bend_insertions": added,
        "bend_removals": removed,
        "geometry_checks": checked,
        "constraints_during_relaxation": project is not None,
        "max_constraint_residual": max_constraint_residual,
        "initial_crossings": len(initial_crossings),
        "minimum_movement_fraction": min_scale,
        "movement_certificate": movement_certificates,
        "seconds": time.monotonic() - started,
        "geometry": validate_geometry(state["positions"], state["routes"]),
    }
    return state


def cubic_routes(state, rounded=True):
    """Round corners only within their empty corner triangles.

    Every curve's convex hull stays inside the checked corner triangle; other
    routes are outside. Native renderer uses these cubics verbatim. A sampled
    independent audit additionally checks the resulting complete carriers.
    """
    positions = state["positions"]
    routes = state["routes"]
    all_segments = [(a, b) for route in routes.values() for a, b in pairwise(route)]
    curves = {}
    corners = 0
    for edge, route in routes.items():
        points = [np.array(positions[k]) for k in route]
        pieces = []
        current = points[0]
        for i in range(1, len(points) - 1):
            prev, center, nxt = points[i - 1 : i + 2]
            left = prev - center
            right = nxt - center
            radius = (
                min(np.linalg.norm(left), np.linalg.norm(right)) * 0.35
                if rounded
                else 0
            )
            # Distance to every nonincident segment bounds a disk around corner.
            for a, b in all_segments:
                if route[i] in (a, b):
                    continue
                pa, pb = np.array(positions[a]), np.array(positions[b])
                v = pb - pa
                t = np.clip(np.dot(center - pa, v) / max(np.dot(v, v), 1e-20), 0, 1)
                radius = min(radius, 0.35 * np.linalg.norm(center - (pa + t * v)))
            entry = center + radius * left / max(np.linalg.norm(left), 1e-15)
            leave = center + radius * right / max(np.linalg.norm(right), 1e-15)
            if np.linalg.norm(current - entry) > 1e-10:
                pieces.append(
                    [
                        current,
                        current + (entry - current) / 3,
                        current + 2 * (entry - current) / 3,
                        entry,
                    ]
                )
            if radius > 1e-9:
                pieces.append(
                    [
                        entry,
                        entry + (center - entry) * 2 / 3,
                        leave + (center - leave) * 2 / 3,
                        leave,
                    ]
                )
                corners += 1
            current = leave
        end = points[-1]
        if np.linalg.norm(end - current) > 1e-10:
            pieces.append(
                [
                    current,
                    current + (end - current) / 3,
                    current + 2 * (end - current) / 3,
                    end,
                ]
            )
        curves[edge] = [[point.tolist() for point in cubic] for cubic in pieces]
    return curves


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("diagram", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--fixed", action="store_true")
    args = parser.parse_args()
    initial = initialize(json.loads(args.diagram.read_text()), scale=2.4)
    result = solve(initial, steps=args.steps, flexible=not args.fixed)
    if result["report"]["planar"]:
        result["route_cubics"] = cubic_routes(result)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["report"], indent=2))
