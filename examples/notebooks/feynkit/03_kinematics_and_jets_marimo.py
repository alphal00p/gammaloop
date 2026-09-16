# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo==0.24.0",
#     "symbolica==3.0.0",
#     "typst==0.15.0",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="FeynKit: Kinematics and jets")


@app.cell
def _():
    from functools import partial

    import marimo as mo

    table = partial(
        mo.ui.table,
        pagination=False,
        selection=None,
        show_download=False,
    )
    return mo, table


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Relativistic kinematics and jet clustering

    FeynKit's Python kinematics use double precision, caller-defined energy units, component order $(E,p_x,p_y,p_z)$, and the mostly-minus metric $(+,-,-,-)$.
    """)
    return


@app.cell
def _():
    from math import cos, cosh, sin, sinh

    import symbolica.community.feynkit as fk

    return cos, cosh, fk, sin, sinh


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Three- and four-momenta

    Construct an on-shell four-vector from a spatial momentum and a mass. The invariant mass check below makes the metric convention explicit.
    """)
    return


@app.cell
def _(fk):
    spatial = fk.ThreeMomentum(3.0, 4.0, 0.0)
    p = spatial.on_shell(12.0)
    {
        "components": p.components(),
        "pT": p.pt,
        "mass_squared": p.mass_squared,
        "mass": p.mass,
        "rapidity": p.rapidity,
    }
    return (p,)


@app.cell
def _(fk, p):
    q = fk.FourMomentum(5.0, 1.0, 2.0, 3.0)
    expected = p.energy * q.energy - p.px * q.px - p.py * q.py - p.pz * q.pz
    {"dot_product": p.dot(q), "component_formula": expected}
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Lorentz transformations

    Boosts require a dimensionless three-velocity with norm below one. Rotations and boosts provide explicit inverse operations, which makes invariant checks straightforward.
    """)
    return


@app.cell
def _(fk):
    rest = fk.FourMomentum(5.0, 0.0, 0.0, 0.0)
    boost = fk.Boost(fk.ThreeMomentum(0.6, 0.0, 0.0))
    boosted = boost.apply(rest)
    recovered = boost.apply_inverse(boosted)
    {
        "rest": rest.components(),
        "boosted": boosted.components(),
        "recovered": recovered.components(),
        "invariant_mass_squared": boosted.mass_squared,
    }
    return


@app.cell
def _(fk):
    unit_x = fk.ThreeMomentum(1.0, 0.0, 0.0)
    unit_y = fk.Rotation.quarter_turn(fk.Axis.Z).apply_three(unit_x)
    unit_y
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Angular distances

    Azimuth is wrapped consistently, and `delta_r` uses rapidity for four-momenta (pseudorapidity for three-momenta).
    """)
    return


@app.cell
def _(fk):
    a = fk.FourMomentum(20.0, 10.0, 0.0, 10.0)
    b = fk.FourMomentum(20.0, 0.0, 10.0, -10.0)
    delta_phi = a.delta_phi(b)
    delta_r = a.delta_r(b)
    {
        "delta_phi": delta_phi,
        "delta_R": delta_r,
    }
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Generalized-$k_T$ jets

    The sequential-recombination algorithms use E-scheme four-vector addition. Results are sorted by decreasing transverse momentum, and constituent indices refer to positions in the input list.
    """)
    return


@app.cell
def _(cos, cosh, fk, sin, sinh):
    particles = [
        fk.FourMomentum(
            _pt * cosh(_rapidity),
            _pt * cos(_phi),
            _pt * sin(_phi),
            _pt * sinh(_rapidity),
        )
        for _pt, _rapidity, _phi in [
            (80.0, 0.20, 0.10),
            (25.0, 0.25, 0.18),
            (35.0, -1.10, 2.40),
        ]
    ]
    return (particles,)


@app.cell
def _(fk, particles):
    jet_definition = fk.JetDefinition.anti_kt(radius=0.6, minimum_pt=5.0)
    clustered = jet_definition.cluster(particles)
    [
        {
            "constituents": jet.constituent_indices,
            "pT": jet.pt,
            "rapidity": jet.rapidity,
            "phi": jet.phi,
        }
        for jet in clustered.jets
    ]
    return (clustered,)


@app.cell
def _(fk, particles):
    definitions = {
        "kt": fk.JetDefinition.kt(0.6, minimum_pt=5.0),
        "cambridge_aachen": fk.JetDefinition.cambridge_aachen(0.6, minimum_pt=5.0),
        "anti_kt": fk.JetDefinition.anti_kt(0.6, minimum_pt=5.0),
    }
    algorithm_constituents = {
        algorithm: [
            jet.constituent_indices for jet in definition.cluster(particles).jets
        ]
        for algorithm, definition in definitions.items()
    }
    algorithm_constituents
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Practical checklist

    - Keep one unit system throughout an event.
    - Check `mass_squared` after constructing or transforming momenta.
    - Choose the jet radius and minimum $p_T$ in the context of the analysis.
    - Preserve the input ordering when interpreting `constituent_indices`.
    - Catch `KinematicsError` around user-provided configurations or non-finite event data.
    """)
    return


@app.cell(hide_code=True)
def _(clustered, mo, particles):
    mo.Html(
        f'<p data-notebook-ready="03_kinematics_and_jets_marimo">Clustered {len(particles)} input momenta into {len(clustered.jets)} anti-kt jets.</p>'
    )
    return


if __name__ == "__main__":
    app.run()
