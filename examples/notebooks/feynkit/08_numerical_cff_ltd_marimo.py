# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo==0.24.0", "numpy", "matplotlib"]
# ///

# ruff: noqa: B018, PLR1711 -- marimo uses final expressions and explicit cell returns.

import marimo

__generated_with = "0.24.0"
app = marimo.App(
    width="medium", app_title="From a loop diagram to a numerical integral"
)


@app.cell
def _():
    from decimal import Decimal, localcontext
    from io import BytesIO
    from time import perf_counter

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from symbolica import E, NumericalIntegrator, S, Symbol
    from symbolica.community import hepkit as hep
    from symbolica.community.hepkit import oneloop

    return (
        BytesIO,
        Decimal,
        E,
        NumericalIntegrator,
        S,
        Symbol,
        hep,
        localcontext,
        mo,
        np,
        oneloop,
        perf_counter,
        plt,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # From a loop diagram to a numerical integral

    **The same integral, different representations, different numerical behavior.**

    A scalar triangle is enough to explore loop-tree duality (LTD), cross-free
    families (CFF), physical thresholds, and adaptive Monte Carlo integration.
    Hepkit generates the graph and CFF; Symbolica evaluates and samples it;
    OneLOop supplies an independent reference.

    We study the scalar integral with numerator one, stripping the model's
    couplings and graph weights. The scalar model supplies the topology;
    external momenta below are independent off-shell inputs.

    This notebook needs a **full native Symbolica community host** containing
    `hepkit` and `hepkit.oneloop`, plus Marimo, NumPy, and Matplotlib. Run it in
    that environment; an ordinary Symbolica wheel is not sufficient.
    """)
    return


@app.cell
def _(hep):
    model = hep.Model.phi3()
    generated = model.process(["phi"], ["phi", "phi"]).generate_diagrams(
        loops=1, max_vertices=3, maximum_bridges=0, progress=None
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    cff = diagram.build_cff()
    diagram.render(momenta=True)
    return cff, diagram


@app.cell
def _(cff, diagram, mo):
    mo.vstack(
        [
            mo.md("## Read the graph's routing and generated surfaces"),
            mo.ui.table(
                [
                    {
                        "edge": edge.id,
                        "momentum": diagram.loop_momentum_basis.edge_signatures[
                            edge.id
                        ].format_momentum(),
                    }
                    for edge in diagram.internal_edges
                ],
                pagination=False,
                selection=None,
            ),
            mo.md(
                f"**{len(cff.orientations)} acyclic orientations** contribute to this triangle."
            ),
            cff.to_expression(expand_surfaces=True).formatted(),
        ]
    )
    return


@app.cell
def _(E, S, Symbol, cff, diagram):
    energies = S("triangle_demo::E0", "triangle_demo::E1", "triangle_demo::E2")
    p_energy, q_energy = S("triangle_demo::p0", "triangle_demo::q0")
    ose, momentum, cind = S("gammalooprs::OSE", "gammalooprs::Q", "spenso::cind")
    external_energies = [p_energy + q_energy, p_energy, q_energy]
    cff_denominators = cff.to_expression(expand_surfaces=True)
    shift_coefficients = []
    for _edge, _energy in zip(diagram.internal_edges, energies):
        _signature = diagram.loop_momentum_basis.edge_signatures[_edge.id]
        assert len(_signature.loops) == 1 and abs(_signature.loops[0]) == 1
        # D(-k-r) = D(k+r): normalize each scalar propagator to +k.
        shift_coefficients.append(
            [_coefficient * _signature.loops[0] for _coefficient in _signature.external]
        )
        cff_denominators = cff_denominators.replace(ose(_edge.id), _energy)
    for _edge in diagram.external_edges:
        cff_denominators = cff_denominators.replace(
            momentum(_edge.id, cind(0)), external_energies[_edge.external_index]
        )
    energy_shifts = [
        sum((_c * _e for _c, _e in zip(_row, external_energies)), E("0"))
        for _row in shift_coefficients
    ]

    # Define a positive below-threshold scalar integral explicitly. The generated
    # CFF here is bare: supply 1/prod(2E) and the spatial measure exactly once.
    cff_integrand = cff_denominators / (
        8 * energies[0] * energies[1] * energies[2] * (2 * Symbol.PI) ** 3
    )
    ltd_cuts = []
    for _i in range(3):
        _cut = 1 / (2 * energies[_i] * (2 * Symbol.PI) ** 3)
        for _j in range(3):
            if _i != _j:
                _cut /= (
                    energies[_i] + energy_shifts[_j] - energy_shifts[_i]
                ) ** 2 - energies[_j] ** 2
        ltd_cuts.append(_cut)
    # Exact algebraic certificate, before floating-point kinematics or sampling.
    assert (cff_integrand - sum(ltd_cuts)).together() == E("0")
    evaluator_parameters = [*energies, p_energy, q_energy]
    cff_evaluator = cff_integrand.evaluator(evaluator_parameters, n_cores=1)
    ltd_evaluators = [
        cut.evaluator(evaluator_parameters, n_cores=1) for cut in ltd_cuts
    ]
    return (
        cff_evaluator,
        cff_integrand,
        energies,
        evaluator_parameters,
        ltd_cuts,
        ltd_evaluators,
        shift_coefficients,
    )


@app.cell(hide_code=True)
def _(ltd_cuts, mo):
    mo.vstack(
        [
            mo.md(r"""
        ## Three residues, six causal terms

        Write each propagator as $D_i=(k+r_i)^2-m^2+i0$ and
        $E_i=\sqrt{|\mathbf k+\mathbf r_i|^2+m^2}$.
        Closing the energy contour below the real axis selects
        $k^0=E_i-r_i^0$. Each row below is one LTD residue.

        The cell above verifies **LTD − CFF = 0 exactly** with Symbolica.
        Numerical stability is a separate question: the individual residues
        can be much larger than their sum.

        Our normalization is
        $I=\int d^4k\,i/[(2\pi)^4D_0D_1D_2]$,
        with spatial measure included in the displayed kernels.
        The OneLOop comparison therefore uses $I=-C_0/(16\pi^2)$.
        """),
            mo.accordion(
                {f"LTD cut {i}": cut.formatted() for i, cut in enumerate(ltd_cuts)}
            ),
        ]
    )
    return


@app.cell
def _(mo):
    regime = mo.ui.dropdown(
        ["On-shell external legs", "Spacelike external legs"],
        value="On-shell external legs",
        label="Kinematics",
    )
    mass = mo.ui.slider(0, 2, step=0.05, value=1, label="Internal mass m")
    map_scale = mo.ui.slider(0.1, 5, step=0.1, value=1, label="Coordinate-map scale")
    channel_power = mo.ui.slider(0, 4, step=0.5, value=2, label="Channel exponent α")
    mo.vstack(
        [
            mo.md("## Choose the numerical experiment"),
            mo.hstack([regime, mass]),
            mo.hstack([map_scale, channel_power]),
            mo.md(
                "For the massless Euclidean example, select **Spacelike external legs** and set **m = 0**."
            ),
        ]
    )
    return channel_power, map_scale, mass, regime


@app.cell
def _(
    Decimal,
    cff_evaluator,
    cff_integrand,
    evaluator_parameters,
    localcontext,
    ltd_evaluators,
    np,
    shift_coefficients,
):
    class Triangle:
        """Numerical kinematics and cube maps for this generated scalar graph."""

        def __init__(self, spacelike, mass):
            self.mass = mass
            self.p = np.array([0.0, 0.0, 0.0, 5.0] if spacelike else [0.5, 0, 0, 0.5])
            self.q = np.array([1.0, 4.0, 3.0, 2.0] if spacelike else [0.5, 0, 0, -0.5])
            self.shifts = np.asarray(shift_coefficients) @ np.array(
                [self.p + self.q, self.p, self.q]
            )
            self.integrable = spacelike or mass > 0.5

        @staticmethod
        def square(vector):
            return float(vector[0] ** 2 - vector[1:] @ vector[1:])

        def inputs(self, points):
            points = np.atleast_2d(points)
            shifted = points[:, None, :] + self.shifts[None, :, 1:]
            on_shell = np.sqrt(np.sum(shifted**2, axis=2) + self.mass**2)
            return np.column_stack(
                [
                    on_shell,
                    np.full(len(points), self.p[0]),
                    np.full(len(points), self.q[0]),
                ]
            )

        def evaluate(self, points):
            return np.asarray(cff_evaluator.evaluate(self.inputs(points))).reshape(-1)

        def cuts(self, points):
            values = self.inputs(points)
            return np.column_stack(
                [
                    np.asarray(evaluator.evaluate(values)).reshape(-1)
                    for evaluator in ltd_evaluators
                ]
            )

        def precise(self, point):
            # Recompute the square roots at high precision too; promoting rounded
            # binary64 energies would not provide an independent stability check.
            with localcontext() as context:
                context.prec = 60
                on_shell = [
                    (
                        sum(
                            (Decimal(str(k)) + Decimal(str(r))) ** 2
                            for k, r in zip(point, shift[1:])
                        )
                        + Decimal(str(self.mass)) ** 2
                    ).sqrt()
                    for shift in self.shifts
                ]
                values = [*on_shell, Decimal(str(self.p[0])), Decimal(str(self.q[0]))]
                result = cff_integrand.evaluate(
                    dict(zip(evaluator_parameters, values)), decimal_digit_precision=50
                )
                return float(result.real)

        def map(self, coordinates, scale, cartesian=False, channels=None, alpha=2):
            x = np.asarray(coordinates)
            if np.any((x <= 0) | (x >= 1)):
                raise ValueError(
                    "The coordinate maps require points strictly inside (0, 1)^3"
                )
            if cartesian:
                points = scale * np.log(x / (1 - x))
                jacobian = np.prod(scale / (x * (1 - x)), axis=1)
            else:
                radius = scale * x[:, 0] / (1 - x[:, 0])
                cosine, azimuth = 2 * x[:, 1] - 1, 2 * np.pi * x[:, 2]
                sine = np.sqrt(1 - cosine**2)
                points = radius[:, None] * np.column_stack(
                    [sine * np.cos(azimuth), sine * np.sin(azimuth), cosine]
                )
                jacobian = 4 * np.pi * radius**2 * scale / (1 - x[:, 0]) ** 2
            if channels is not None:
                points -= self.shifts[channels, 1:]
                on_shell = self.inputs(points)[:, :3]
                factors = on_shell ** (-alpha)
                jacobian *= factors[np.arange(len(points)), channels] / np.sum(
                    factors, axis=1
                )
            return self.evaluate(points) * jacobian

    return (Triangle,)


@app.cell
def _(Triangle, mass, mo, regime):
    triangle = Triangle(regime.value == "Spacelike external legs", mass.value)
    mo.ui.table(
        [
            {"invariant": "p²", "value": triangle.square(triangle.p)},
            {"invariant": "q²", "value": triangle.square(triangle.q)},
            {"invariant": "(p+q)²", "value": triangle.square(triangle.p + triangle.q)},
            {"invariant": "m²", "value": triangle.mass**2},
        ],
        pagination=False,
        selection=None,
    )
    return (triangle,)


@app.cell
def _(BytesIO, mo, np, plt, triangle):
    _radii = np.geomspace(0.03, 1e7, 90)
    _direction = np.array([0.31, 0.47, 0.83])
    _points = _radii[:, None] * (_direction / np.linalg.norm(_direction))
    _reference = np.array([triangle.precise(point) for point in _points])
    _cuts = triangle.cuts(_points)
    _cff = triangle.evaluate(_points)
    _ltd = _cuts.sum(axis=1)
    _figure, _axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for _i in range(3):
        _axes[0].loglog(_radii, np.abs(_cuts[:, _i]), alpha=0.65, label=f"LTD cut {_i}")
    _axes[0].loglog(
        _radii, np.abs(_reference), color="black", label="Sum · 50-digit reference"
    )
    _axes[0].set(xlabel="|k| along a fixed ray", ylabel="Absolute integrand")
    for _values, _label in [(_ltd, "LTD sum · float64"), (_cff, "CFF · float64")]:
        _relative = np.abs((_values - _reference) / _reference)
        _axes[1].loglog(_radii, np.maximum(_relative, 1e-17), label=_label)
    _axes[1].set(xlabel="|k| along the same ray", ylabel="Relative error")
    for _axis in _axes:
        _axis.legend(fontsize=8)
        _axis.grid(alpha=0.2)
    _image = BytesIO()
    _figure.savefig(_image, format="png", dpi=140)
    _plot = mo.image(_image.getvalue())
    plt.close(_figure)
    mo.vstack(
        [
            mo.md("## Equivalent algebra, different floating-point behavior"),
            _plot,
            mo.md(
                "The reference recomputes energies at 60 digits and evaluates CFF at 50 digits. "
                "Errors below $10^{-17}$ are placed on the plot floor; the reference is rounded "
                "to float64 only for this comparison. Inspect the growing cancellation between cuts."
            ),
        ]
    )
    return


@app.cell
def _(BytesIO, mo, np, plt, triangle):
    _extent = 7 if triangle.p[0] == 0 else 2
    _grid = np.linspace(-_extent, _extent, 180)
    _xx, _zz = np.meshgrid(_grid, _grid)
    _points = np.column_stack([_xx.ravel(), np.zeros(_xx.size), _zz.ravel()])
    _energies = triangle.inputs(_points)[:, :3]
    _figure, _axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    _counts = [0, 0]
    for _i in range(3):
        for _j in range(3):
            if _i == _j:
                continue
            _shift = triangle.shifts[_j, 0] - triangle.shifts[_i, 0]
            for _kind, _surface in enumerate(
                [
                    _energies[:, _i] - _energies[:, _j] + _shift,
                    _energies[:, _i] + _energies[:, _j] + _shift,
                ]
            ):
                if _surface.min() < 0 < _surface.max():
                    _axes[_kind].contour(
                        _xx, _zz, _surface.reshape(_xx.shape), levels=[0]
                    )
                    _counts[_kind] += 1
    for _axis, _title, _count in zip(
        _axes,
        ["H-surfaces · cancel in the sum", "E-surfaces · physical thresholds"],
        _counts,
    ):
        _axis.set(title=_title, xlabel="kx", ylabel="kz", aspect="equal")
        if not _count:
            _axis.text(
                0.5,
                0.5,
                "No zero contour in this slice",
                transform=_axis.transAxes,
                ha="center",
            )
    _image = BytesIO()
    _figure.savefig(_image, format="png", dpi=140)
    _plot = mo.image(_image.getvalue())
    plt.close(_figure)
    mo.vstack(
        [
            mo.md("## Where do the denominators vanish?"),
            _plot,
            mo.md(
                "Slice: $k_y=0$. Spacelike kinematics expose cancelling H-surfaces. "
                "With on-shell external legs, lower the mass through $m=1/2$ to see "
                "the physical threshold open. Degenerate zeros exactly at threshold "
                "need not appear as contour lines."
            ),
        ]
    )
    return


@app.cell
def _(NumericalIntegrator, np, perf_counter):
    class IntegrationExperiment:
        """Compare independent production samples after optional grid training."""

        methods = (
            "Uniform Cartesian",
            "Uniform spherical",
            "Adaptive spherical",
            "Adaptive channels",
        )

        @classmethod
        def run(cls, triangle, scale, alpha, batch_size=3000, batches=10, seed=2026):
            if not triangle.integrable:
                raise ValueError("This experiment requires threshold-free kinematics")
            histories, summaries, distributions = [], [], {}
            for method_index, method in enumerate(cls.methods):
                adaptive, channels = method_index >= 2, method_index == 3
                grids = [
                    NumericalIntegrator.continuous(
                        3, n_bins=32, min_probability_density=0.05
                    )
                    for _ in range(3 if channels else 1)
                ]
                grid = NumericalIntegrator.discrete(grids) if channels else grids[0]
                rng = NumericalIntegrator.rng(seed, method_index)
                warmup = 2 if adaptive else 0
                count, mean, m2, maximum = 0, 0.0, 0.0, 0.0
                retained = []
                start = perf_counter()
                for iteration in range(batches):
                    samples = grid.sample(batch_size, rng)
                    coordinates = np.array([sample.c for sample in samples])
                    selected = (
                        np.array([sample.d[0] for sample in samples])
                        if channels
                        else None
                    )
                    values = triangle.map(
                        coordinates,
                        scale,
                        cartesian=method_index == 0,
                        channels=selected,
                        alpha=alpha,
                    )
                    if not np.all(np.isfinite(values)):
                        raise FloatingPointError(f"Non-finite integrand in {method}")
                    if iteration < warmup:
                        grid.add_training_samples(samples, values.tolist())
                        grid.update(1.5, 1.5)
                        continue
                    # Weights are cumulative from each layer downwards:
                    # weights[0] includes the full discrete/continuous density.
                    # The physical channel partition was applied in triangle.map.
                    weights = values * np.array(
                        [sample.weights[0] for sample in samples]
                    )
                    batch_mean = float(weights.mean())
                    batch_m2 = float(np.sum((weights - batch_mean) ** 2))
                    delta = batch_mean - mean
                    total = count + batch_size
                    m2 += batch_m2 + delta**2 * count * batch_size / total
                    mean += delta * batch_size / total
                    count = total
                    error = np.sqrt(m2 / (count * (count - 1)))
                    maximum = max(maximum, float(np.max(np.abs(weights))))
                    retained.append(weights)
                    histories.append(
                        {
                            "method": method,
                            "evaluations": (iteration + 1) * batch_size,
                            "production": count,
                            "mean": mean,
                            "error": float(error),
                        }
                    )
                distributions[method] = np.concatenate(retained)
                summaries.append(
                    {
                        "method": method,
                        "estimate": mean,
                        "standard error": float(error),
                        "evaluations": batches * batch_size,
                        "pilot": warmup * batch_size,
                        "production": count,
                        "max |weight|": maximum,
                        "relative spread": float(np.sqrt(count) * error / abs(mean)),
                        "seconds": perf_counter() - start,
                    }
                )
            return histories, summaries, distributions

    return (IntegrationExperiment,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Spend the same sample budget four ways

    The Cartesian map uses $k_a=s\log[x_a/(1-x_a)]$. The spherical map uses
    $r=sx_0/(1-x_0)$, $\cos\theta=2x_1-1$, and $\phi=2\pi x_2$.
    Both include their Jacobians.

    For multiple channels, center one spherical map at each propagator's
    spatial shift and split the integrand with
    $w_i=E_i^{-\alpha}/\sum_j E_j^{-\alpha}$, so $\sum_i w_i=1$.
    Symbolica learns both the continuous grids and the discrete channel
    probabilities. Coincident centers are retained here to keep the direct
    connection to the three propagators visible.

    **Each method gets the same total number of evaluations.** Adaptive methods
    spend two batches training, then freeze their grids. Reported estimates and
    standard errors use only independent production samples; training cost is
    included on the horizontal axis. One-standard-error bars are statistical
    estimates, not rigorous error bounds.
    """)
    return


@app.cell
def _(mo):
    integration_settings = mo.ui.dictionary(
        {
            "batch_size": mo.ui.dropdown(
                [1000, 3000, 10000], value=3000, label="Points per batch"
            ),
            "batches": mo.ui.slider(4, 20, value=10, label="Batches"),
            "seed": mo.ui.number(0, 1_000_000, value=2026, label="Random seed"),
        }
    ).form(submit_button_label="Run integration comparison")
    integration_settings
    return (integration_settings,)


@app.cell
def _(
    IntegrationExperiment,
    channel_power,
    integration_settings,
    map_scale,
    mo,
    np,
    oneloop,
    triangle,
):
    mo.stop(
        not triangle.integrable,
        mo.md(
            "**Choose m > 1/2 or spacelike external legs to integrate.** The real-axis "
            "sampler here does not regularize physical thresholds; the geometry and "
            "stability plots above remain available."
        ),
    )
    mo.stop(
        integration_settings.value is None,
        mo.md("Select **Run integration comparison** to sample."),
    )
    _settings = integration_settings.value
    _finite, _pole, _double_pole = oneloop.c0(
        triangle.square(triangle.p),
        triangle.square(triangle.q),
        triangle.square(triangle.p + triangle.q),
        triangle.mass**2,
        triangle.mass**2,
        triangle.mass**2,
        backend="expression",
    )
    assert abs(_pole) < 1e-12 and abs(_double_pole) < 1e-12
    assert abs(_finite.imag) < 1e-10 * max(1, abs(_finite.real))
    reference = -_finite.real / (16 * np.pi**2)
    with mo.status.spinner(
        title="Sampling uniform, adaptive, and multichannel integrands"
    ):
        histories, summaries, distributions = IntegrationExperiment.run(
            triangle, map_scale.value, channel_power.value, **_settings
        )
    for _row in summaries:
        _row["deviation / error"] = (_row["estimate"] - reference) / _row[
            "standard error"
        ]
    mo.vstack(
        [
            mo.md(f"**Independent OneLOop reference:** {reference:.12g}"),
            mo.ui.table(summaries, pagination=False, selection=None),
        ]
    )
    return distributions, histories, reference, summaries


@app.cell
def _(BytesIO, IntegrationExperiment, distributions, histories, mo, np, plt, reference):
    _figure, _axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    for _method in IntegrationExperiment.methods:
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
        _weights = _weights[_weights > 0]
        _axes[2].hist(
            np.log10(_weights), bins=45, density=True, histtype="step", label=_method
        )
    _axes[0].axhline(0, color="black", linewidth=0.8)
    _axes[0].set(
        xlabel="Evaluations, including training", ylabel="Estimate / reference − 1"
    )
    _axes[0].legend(fontsize=7)
    _axes[1].set(
        xlabel="Evaluations, including training", ylabel="Relative standard error"
    )
    _axes[2].set(xlabel="log₁₀(|sample weight| / |reference|)", ylabel="Density")
    _image = BytesIO()
    _figure.savefig(_image, format="png", dpi=140)
    _plot = mo.image(_image.getvalue())
    plt.close(_figure)
    mo.vstack(
        [
            mo.md("## Convergence and weight distributions"),
            _plot,
            mo.md(
                "Repeat with another seed, change the map scale, then try massless spacelike "
                "kinematics. Which improvement survives all three changes? A faster kernel "
                "and a lower-variance sampling distribution improve different parts of the calculation."
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
