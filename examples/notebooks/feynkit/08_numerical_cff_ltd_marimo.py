# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo==0.24.0", "docstring-to-markdown>=0.17,<1", "ty>=0.0.75", "numpy", "matplotlib"]
# ///

# ruff: noqa: B018, PLR1711 -- marimo uses final expressions and explicit cell returns.

import marimo

__generated_with = "0.24.0"
app = marimo.App(
    width="medium",
    app_title="From a loop diagram to a numerical integral",
)


@app.cell(hide_code=True)
def _():
    from math import prod
    from time import perf_counter

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from symbolica import (
        E,
        Expression,
        Float,
        FunctionDefinition,
        NumericalIntegrator,
        Replacement,
        S,
        Symbol,
    )
    from symbolica.community import hepkit as hep
    from symbolica.community.hepkit import oneloop

    return (
        E,
        Expression,
        Float,
        FunctionDefinition,
        NumericalIntegrator,
        Replacement,
        S,
        Symbol,
        hep,
        mo,
        np,
        oneloop,
        perf_counter,
        plt,
        prod,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # From a loop diagram to a numerical integral

    One scalar triangle, two energy representations, four sampling strategies.
    Hepkit generates **CFF and LTD**; Symbolica evaluates and integrates them;
    OneLOop provides an independent reference.

    Inspect the orientations, trees and surfaces below, then change the
    kinematics to explore cancellation and convergence. The two- and three-loop
    gallery also lets you explore larger families and their overlapping regions.
    """)
    return


@app.cell
def _(hep):
    model = hep.Model.phi3()
    diagram = (
        model.process(["phi"], ["phi", "phi"])
        .generate_diagrams(loops=1, max_vertices=3, maximum_bridges=0, progress=None)
        .diagrams[0]
    )
    cff = diagram.integrate_energy(method="cff")
    ltd = diagram.integrate_energy(method="ltd")
    diagram.render(momenta=True)
    return cff, diagram, ltd


@app.cell
def _(cff):
    cff
    return


@app.cell
def _(ltd):
    ltd
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Polynomial energy numerators

    CFF also integrates an explicit polynomial numerator. Bounds are inferred in
    the **physical edge energies**, before routing. A family now includes its
    evaluated numerator and coefficient; hover a contribution to inspect the
    energy substitutions. Powers denote multiplicity of the same surface.
    """)
    return


@app.cell
def _(S, diagram):
    _Q, _cind = S("gammalooprs::Q", "spenso::cind")
    numerator_edge = diagram.internal_edges[0].id
    quadratic_cff = diagram.integrate_energy(
        method="cff", numerator=_Q(numerator_edge, _cind(0)) ** 2
    )
    quadratic_cff
    return (quadratic_cff,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Two- and three-loop surfaces

    Choose a topology, then open its CFF and LTD explorers. **Shift-click CFF
    surface factors** to overlay several regions; switch families to compare
    nested and disjoint circlings. In LTD, click a cut to navigate trees, click
    a tree edge to switch its surface pair, and hover a cut to trace its pole sign.

    The boxes have four external legs; the tetrahedron is a vacuum graph.
    These examples use symbolic energies. The numerical section below evaluates
    the one-loop triangle.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    multiloop_topology = mo.ui.dropdown(
        {
            "Double box · 2 loops": """digraph double_box {
                edge [particle="phi"];
                ext [style=invis];
                ext -> a; ext -> d; c -> ext; f -> ext;
                a -> b; b -> c; d -> e; e -> f;
                a -> d; b -> e [lmb_id=0]; c -> f [lmb_id=1];
            }""",
            "Crossed double box · 2 loops": """digraph crossed_double_box {
                edge [particle="phi"];
                ext [style=invis];
                ext -> x; ext -> y; u -> ext; v -> ext;
                a -> x; x -> b [lmb_id=0];
                a -> y; y -> b [lmb_id=1];
                a -> u; u -> v; v -> b;
            }""",
            "Triple box · 3 loops": """digraph triple_box {
                edge [particle="phi"];
                ext [style=invis];
                ext -> a; ext -> e; d -> ext; h -> ext;
                a -> b; b -> c; c -> d;
                e -> f; f -> g; g -> h;
                a -> e; b -> f [lmb_id=0];
                c -> g [lmb_id=1]; d -> h [lmb_id=2];
            }""",
            "Tetrahedron · 3 loops": """digraph tetrahedron {
                edge [particle="phi"];
                a -> b; a -> c; a -> d;
                b -> c [lmb_id=0]; b -> d [lmb_id=1]; c -> d [lmb_id=2];
            }""",
        },
        value="Double box · 2 loops",
        label="Diagram",
        full_width=True,
    )
    multiloop_topology
    return (multiloop_topology,)


@app.cell
def _(hep, multiloop_topology):
    multiloop_diagram = hep.FeynmanDiagram.from_dot(
        hep.Model.phi3(), multiloop_topology.value
    )
    multiloop_cff = multiloop_diagram.integrate_energy(method="cff")
    multiloop_ltd = multiloop_diagram.integrate_energy(method="ltd")
    return multiloop_cff, multiloop_ltd


@app.cell
def _(multiloop_cff):
    multiloop_cff.orientations[5]
    return


@app.cell
def _(multiloop_ltd):
    multiloop_ltd
    return


@app.cell(hide_code=True)
def _(
    E,
    Expression,
    FunctionDefinition,
    Replacement,
    S,
    Symbol,
    cff,
    diagram,
    ltd,
    mo,
):
    loop_variables = S("triangle::kx", "triangle::ky", "triangle::kz")
    _m = S("triangle::m")
    _p = S("triangle::p0", "triangle::px", "triangle::py", "triangle::pz")
    _q = S("triangle::q0", "triangle::qx", "triangle::qy", "triangle::qz")
    _momentum, _cind = S("gammalooprs::Q", "spenso::cind")
    _external = [[a + b for a, b in zip(_p, _q)], _p, _q]
    _on_shell = cff.on_shell_energies
    _energies = [energy.symbol for energy in _on_shell.values()]
    _energy_replacements = [
        Replacement(_momentum(edge.id, _cind(0)), _external[edge.external_index][0])
        for edge in diagram.external_edges
    ]
    _measure = (2 * Symbol.PI) ** -3
    # Both energy representations use dq0/(2*pi*i). The showcase uses
    # i*d4k/(2*pi)^4, hence the common minus sign and spatial measure.
    cff_kernel = (-cff.to_expression() * _measure).replace_multiple(
        _energy_replacements
    )
    ltd_kernels = [
        (-cut * _measure).replace_multiple(_energy_replacements)
        for cut in [residue.to_expression() for residue in ltd.residues]
    ]
    assert (cff_kernel - sum(ltd_kernels)).together().is_zero()

    _k, _p_symbol = S("gammalooprs::K", "gammalooprs::P")
    _basis = diagram.loop_momentum_basis
    _coordinates = [
        Replacement(_k(0, _cind(axis)), value)
        for axis, value in enumerate(loop_variables, start=1)
    ] + [
        Replacement(
            _p_symbol(_basis.external_edges.index(edge.id), _cind(axis)),
            _external[edge.external_index][axis],
        )
        for edge in diagram.external_edges
        for axis in range(1, 4)
    ]
    shift_coefficients, _energy_functions = [], []
    for _edge in diagram.internal_edges:
        _signature = _basis.edge_signatures[_edge.id]
        assert len(_signature.loops) == 1 and abs(_signature.loops[0]) == 1
        # The same shifts locate numerical integration channels.
        shift_coefficients.append(
            [c * _signature.loops[0] for c in _signature.external]
        )
        _energy = _on_shell[_edge.id]
        _definition = _energy.to_expression().replace_multiple(_coordinates)
        _definition = _definition.replace(_edge.particle.mass, _m)
        _energy_functions.append(FunctionDefinition(_energy.symbol, [], _definition))
    evaluator_parameters = [*loop_variables, _m, *_p, *_q]
    _options = {
        "params": evaluator_parameters,
        "functions": _energy_functions,
        "n_cores": 1,
    }
    cff_evaluator = cff_kernel.evaluator(**_options)
    comparison_evaluator = Expression.evaluator_multiple(
        [cff_kernel, *ltd_kernels], **_options
    )
    energy_evaluator = Expression.evaluator_multiple(_energies, **_options)

    # Inspect native surfaces; S and -S have the same zero contour.
    _surfaces = {}
    for _surface in ltd.surfaces:
        _expression = (
            _surface.to_expression().replace_multiple(_energy_replacements).expand()
        )
        _first = next(int(c) for c in _surface.energy_coefficients.values() if c != 0)
        _canonical = (_expression if _first > 0 else -_expression).expand()
        _surfaces[_canonical] = _surface.kind
    surface_expressions = sorted(_surfaces, key=str)
    surface_kinds = [_surfaces[surface] for surface in surface_expressions]
    surface_evaluator = Expression.evaluator_multiple(surface_expressions, **_options)
    mo.md(r"""
    **Symbolica verifies LTD − CFF = 0 exactly.** The numerical kernels include
    the spatial measure and on-shell energy factors, with numerator one.
    We use $I=\int i\,d^4k/[(2\pi)^4D_0D_1D_2]$, so the OneLOop reference is
    $-C_0/(16\pi^2)$.
    """)
    return (
        cff_evaluator,
        comparison_evaluator,
        energy_evaluator,
        shift_coefficients,
        surface_evaluator,
        surface_kinds,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Kinematics and numerical stability
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    regime = mo.ui.dropdown(
        ["On-shell external legs", "Spacelike external legs"],
        value="Spacelike external legs",
        label="Kinematics",
    )
    mass = mo.ui.slider(0, 2, step=0.05, value=1, label="Internal mass m")
    mo.hstack([regime, mass])
    return mass, regime


@app.cell(hide_code=True)
def _(hep, mass, mo, np, regime, shift_coefficients):
    spacelike = regime.value == "Spacelike external legs"
    p = hep.FourMomentum(*([0, 0, 0, 5] if spacelike else [0.5, 0, 0, 0.5]))
    q = hep.FourMomentum(*([1, 4, 3, 2] if spacelike else [0.5, 0, 0, -0.5]))
    integrable = spacelike or mass.value > 0.5
    channel_centers = -(
        np.asarray(shift_coefficients)
        @ np.array([(p + q).components(), p.components(), q.components()])
    )[:, 1:]
    _constants = np.array([mass.value, *p.components(), *q.components()])

    def kinematic_inputs(points):
        points = np.atleast_2d(points)
        return np.column_stack(
            [points, np.broadcast_to(_constants, (len(points), len(_constants)))]
        )

    mo.ui.table(
        [
            {"invariant": label, "value": value}
            for label, value in [
                ("p²", p.mass_squared),
                ("q²", q.mass_squared),
                ("(p+q)²", (p + q).mass_squared),
                ("m²", mass.value**2),
            ]
        ],
        pagination=False,
        selection=None,
    )
    return channel_centers, integrable, kinematic_inputs, p, q, spacelike


@app.cell(hide_code=True)
def _(
    Float,
    cff_evaluator,
    comparison_evaluator,
    kinematic_inputs,
    mo,
    np,
    plt,
):
    _radii = np.geomspace(0.03, 1e7, 90)
    _direction = np.array([0.31, 0.47, 0.83])
    _inputs = kinematic_inputs(
        _radii[:, None] * (_direction / np.linalg.norm(_direction))
    )
    # Increase precision before evaluating the on-shell roots, not afterwards.
    _reference = np.array(
        [
            float(
                cff_evaluator.evaluate_with_prec(
                    [Float(value, decimal_digits=50) for value in point], 50
                )[0]
            )
            for point in _inputs
        ]
    )
    _values = comparison_evaluator.evaluate(_inputs)
    _figure, _axes = plt.subplots(1, 2, figsize=(10, 3.4), layout="constrained")
    for _i, _cut in enumerate(_values[:, 1:].T):
        _axes[0].loglog(_radii, np.abs(_cut), alpha=0.65, label=f"LTD cut {_i}")
    _axes[0].loglog(
        _radii, np.abs(_reference), color="black", label="50-digit CFF reference"
    )
    for _result, _label in [
        (_values[:, 1:].sum(axis=1), "LTD sum"),
        (_values[:, 0], "CFF"),
    ]:
        _axes[1].loglog(
            _radii,
            np.maximum(np.abs((_result - _reference) / _reference), 1e-17),
            label=_label,
        )
    for _axis, _ylabel in zip(
        _axes, ["Absolute integrand", "Relative error · float64"]
    ):
        _axis.set(xlabel="|k| along a fixed ray", ylabel=_ylabel)
        _axis.legend(fontsize=8)
        _axis.grid(alpha=0.2)
    plt.close(_figure)
    mo.vstack(
        [
            _figure,
            mo.md(
                "The same Symbolica evaluator computes the reference at 50 digits, including "
                "the square roots. Only the final comparison is rounded to float64. "
                "The plot floor is $10^{-17}$."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Where the denominators vanish
    """)
    return


@app.cell(hide_code=True)
def _(
    kinematic_inputs,
    mass,
    mo,
    np,
    plt,
    spacelike,
    surface_evaluator,
    surface_kinds,
):
    _extent = 7 if spacelike else 2
    _grid = np.linspace(-_extent, _extent, 180)
    _xx, _zz = np.meshgrid(_grid, _grid)
    _points = np.column_stack([_xx.ravel(), np.zeros(_xx.size), _zz.ravel()])
    _surfaces = surface_evaluator.evaluate(kinematic_inputs(_points))
    _figure, _axes = plt.subplots(1, 2, figsize=(10, 3.4), layout="constrained")
    for _axis, _kind, _title in zip(
        _axes, ["H", "E"], ["H-surfaces · cancelling poles", "E-surfaces · thresholds"]
    ):
        _visible = 0
        for _surface, _surface_kind in zip(_surfaces.T, surface_kinds):
            if _surface_kind == _kind and _surface.min() < 0 < _surface.max():
                _axis.contour(
                    _xx,
                    _zz,
                    _surface.reshape(_xx.shape),
                    levels=[0],
                    colors=["#b7771b" if _kind == "H" else "#16827d"],
                    linewidths=1.8,
                )
                _visible += 1
        _axis.set(
            title=_title,
            xlabel="kx",
            ylabel="kz",
            aspect="equal",
            xlim=(-_extent, _extent),
            ylim=(-_extent, _extent),
        )
        _axis.grid(alpha=0.15)
        if not _visible:
            if spacelike and _kind == "E":
                _reason = "Spacelike external invariants:\nno real E-surface."
            elif not spacelike and _kind == "H" and mass.value > 0:
                _reason = "Lightlike external legs with m > 0:\nno real H-surface."
            elif not spacelike and _kind == "E" and mass.value > 0.5:
                _reason = (
                    "Below threshold: 2m > √s = 1.\nUse m = 0.25 to see the E-surface."
                )
            elif not spacelike and _kind == "E" and mass.value == 0.5:
                _reason = "At threshold: the E-surface shrinks to a point."
            else:
                _reason = "No sign-changing zero contour in this slice."
            _axis.text(
                0.5,
                0.5,
                _reason,
                transform=_axis.transAxes,
                ha="center",
                va="center",
                fontsize=9,
            )
    plt.close(_figure)
    mo.vstack(
        [
            _figure,
            mo.md(
                "Native LTD denominator factors on the $k_y=0$ slice. "
                "**H-contours:** spacelike legs, $m=1$ (the default). "
                "**E-contour:** on-shell legs, $m=0.25$. "
                "A degenerate zero exactly at threshold need not form a contour."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Numerical integration

    Compare uniform Cartesian and spherical maps with Symbolica's adaptive
    spherical grid and discrete channels. Each method gets the same sample
    budget; Symbolica accumulates the estimates and standard errors.

    The Cartesian map is $k_a=s\log[x_a/(1-x_a)]$. The spherical map uses
    $r=sx_0/(1-x_0)$, $\cos\theta=2x_1-1$ and $\phi=2\pi x_2$.
    Channels center that map at each propagator and partition unity with
    $w_i=E_i^{-\alpha}/\sum_j E_j^{-\alpha}$.
    """)
    return


@app.cell(hide_code=True)
def _(Expression, S, Symbol, prod):
    _x = S("map::x0", "map::x1", "map::x2")
    _scale = S("map::scale")
    _cartesian = [_scale * (x / (1 - x)).log() for x in _x]
    _r = _scale * _x[0] / (1 - _x[0])
    _cos, _phi = 2 * _x[1] - 1, 2 * Symbol.PI * _x[2]
    _sin = (1 - _cos**2).sqrt()
    _spherical = [_r * _sin * _phi.cos(), _r * _sin * _phi.sin(), _r * _cos]
    # Differentiate the radial maps in Symbolica; the spherical solid angle is 4π.
    coordinate_maps = {
        "Cartesian": Expression.evaluator_multiple(
            [*_cartesian, prod(k.derivative(x) for k, x in zip(_cartesian, _x))],
            [*_x, _scale],
            n_cores=1,
        ),
        "Spherical": Expression.evaluator_multiple(
            [*_spherical, 4 * Symbol.PI * _r**2 * _r.derivative(_x[0])],
            [*_x, _scale],
            n_cores=1,
        ),
    }
    return (coordinate_maps,)


@app.cell(hide_code=True)
def _(mo):
    integration_settings = mo.ui.dictionary(
        {
            "batch_size": mo.ui.dropdown(
                [1000, 3000, 10000], value=3000, label="Points per batch"
            ),
            "batches": mo.ui.slider(4, 20, value=10, label="Batches"),
            "scale": mo.ui.slider(0.1, 5, step=0.1, value=1, label="Map scale"),
            "alpha": mo.ui.slider(0, 4, step=0.5, value=2, label="Channel exponent α"),
            "seed": mo.ui.number(0, 1_000_000, value=2026, label="Random seed"),
        }
    ).form(submit_button_label="Run integration comparison")
    integration_settings
    return (integration_settings,)


@app.cell(hide_code=True)
def _(
    NumericalIntegrator,
    cff_evaluator,
    channel_centers,
    coordinate_maps,
    energy_evaluator,
    integrable,
    integration_settings,
    kinematic_inputs,
    mass,
    mo,
    np,
    oneloop,
    p,
    perf_counter,
    q,
):
    mo.stop(
        not integrable,
        mo.md(
            "Choose **m > 1/2** or **spacelike legs** to integrate. Physical thresholds "
            "need contour deformation; the real-axis sampler here does not provide it."
        ),
    )
    mo.stop(
        integration_settings.value is None,
        mo.md("Select **Run integration comparison** to sample."),
    )
    _settings = integration_settings.value
    _finite, _pole, _double_pole = oneloop.c0(
        p.mass_squared,
        q.mass_squared,
        (p + q).mass_squared,
        *([mass.value**2] * 3),
        backend="expression",
    )
    assert abs(_pole) < 1e-12 and abs(_double_pole) < 1e-12
    assert abs(_finite.imag) < 1e-10 * max(1, abs(_finite.real))
    reference = -_finite.real / (16 * np.pi**2)
    methods = [
        "Uniform Cartesian",
        "Uniform spherical",
        "Adaptive spherical",
        "Adaptive channels",
    ]
    histories, summaries, distributions = [], [], {}
    with mo.status.spinner(title="Sampling with Symbolica"):
        for _method_index, _method in enumerate(methods):
            _adaptive, _channels = _method_index >= 2, _method_index == 3
            _grids = [
                NumericalIntegrator.continuous(
                    3, n_bins=32 if _adaptive else 1, min_probability_density=0.05
                )
                for _ in range(len(channel_centers) if _channels else 1)
            ]
            _grid = NumericalIntegrator.discrete(_grids) if _channels else _grids[0]
            _rng = NumericalIntegrator.rng(_settings["seed"], _method_index)
            _map = coordinate_maps["Cartesian" if _method_index == 0 else "Spherical"]
            _weights = []
            _start = perf_counter()
            for _iteration in range(_settings["batches"]):
                _samples = _grid.sample(_settings["batch_size"], _rng)
                _coordinates = np.array([sample.c for sample in _samples])
                _mapped = _map.evaluate(
                    np.column_stack(
                        [_coordinates, np.full(len(_samples), _settings["scale"])]
                    )
                )
                _points, _jacobian = _mapped[:, :3].copy(), _mapped[:, 3].copy()
                if _channels:
                    _selected = np.array([sample.d[0] for sample in _samples])
                    _points += channel_centers[_selected]
                    _partition = energy_evaluator.evaluate(
                        kinematic_inputs(_points)
                    ) ** (-_settings["alpha"])
                    _jacobian *= _partition[
                        np.arange(len(_samples)), _selected
                    ] / _partition.sum(axis=1)
                _values = (
                    cff_evaluator.evaluate(kinematic_inputs(_points)).ravel()
                    * _jacobian
                )
                if not np.all(np.isfinite(_values)):
                    raise FloatingPointError(f"Non-finite integrand in {_method}")
                _grid.add_training_samples(_samples, _values.tolist())
                _rate = 1.5 if _adaptive else 0.0
                _mean, _error, _chi2 = _grid.update(_rate, _rate)
                _, _, _, _negative, _positive, _count = _grid.get_live_estimate()
                # Symbolica applies the full inverse density, including the
                # discrete layer. Keep these weights only for the histogram.
                _weights.append(
                    _values * np.array([sample.weights[0] for sample in _samples])
                )
                histories.append(
                    {
                        "method": _method,
                        "evaluations": _count,
                        "mean": _mean,
                        "error": _error,
                    }
                )
            distributions[_method] = np.concatenate(_weights)
            summaries.append(
                {
                    "method": _method,
                    "estimate": _mean,
                    "standard error": _error,
                    "evaluations": _count,
                    "deviation / error": (_mean - reference) / _error,
                    "max |weight|": max(abs(_negative), _positive),
                    "seconds": perf_counter() - _start,
                }
            )
    mo.vstack(
        [
            mo.md(f"**OneLOop reference:** {reference:.12g}"),
            mo.ui.table(summaries, pagination=False, selection=None),
        ]
    )
    return distributions, histories, methods, reference


@app.cell(hide_code=True)
def _(distributions, histories, methods, mo, np, plt, reference):
    _figure, _axes = plt.subplots(1, 3, figsize=(12, 3.4), layout="constrained")
    for _method in methods:
        _history = [row for row in histories if row["method"] == _method]
        _n = [row["evaluations"] for row in _history]
        _axes[0].errorbar(
            _n,
            [row["mean"] / reference - 1 for row in _history],
            yerr=[row["error"] / abs(reference) for row in _history],
            marker=".",
            capsize=2,
            label=_method,
        )
        _axes[1].loglog(
            _n, [row["error"] / abs(reference) for row in _history], marker="."
        )
        _weights = np.abs(distributions[_method]) / abs(reference)
        _axes[2].hist(
            np.log10(_weights[_weights > 0]),
            bins=45,
            density=True,
            histtype="step",
            label=_method,
        )
    _axes[0].axhline(0, color="black", linewidth=0.8)
    _axes[0].set(xlabel="Evaluations", ylabel="Estimate / reference − 1")
    _axes[0].legend(fontsize=7)
    _axes[1].set(xlabel="Evaluations", ylabel="Relative standard error")
    _axes[2].set(xlabel="log₁₀(|weight| / |reference|)", ylabel="Density")
    plt.close(_figure)
    mo.vstack(
        [
            _figure,
            mo.md(
                "All batches contribute; adaptive grids update between batches. Error bars "
                "are statistical estimates. Repeat with another seed or map scale, then "
                "try **massless spacelike** kinematics."
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
