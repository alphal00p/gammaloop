# Generic sampling, channel selection, and cut/left/right composition

Research addendum, 12 September 2026. This updates the API and coverage proposal
in [REPORT.md](REPORT.md); it does not change GammaLoop code. In particular, the
earlier `channel_mode`/`channel_partition` names and assumption that ordinary LMB
channels must remain selected are superseded below.

## Bounded GL638 weights: the precise claim

At the measured regular H/Z corner, use local normal coordinates u=H, v=alpha Z,
with a fixed positive metric alpha, and R=sqrt(u^2+v^2). The measured leading
behavior is F~C(theta,y)/R. On a compact patch where the chart Jacobian and C are
bounded and the tangential density is bounded below, either q_uv~1/R or
q_uv~1/sqrt(|uv|) gives bounded leading event weights. For the second choice,
sqrt(|uv|)/R <= 1/sqrt(2). These densities are locally normalizable.

A mixture of separate half-power H and half-power Z channels still scales as
R^(-1/2) along ordinary rays and leaves unbounded R^(-1/2) weights. The two
normals must be targeted jointly. The same argument applies to the exact native
A-star pullback described in the main report.

This establishes a remedy for the identified regular singularity. It does not
prove a global bound for GL638, even after excluding tangencies: other star
images, hierarchical soft limits, shrinking fibers, center transitions and UV
tails require coverage and appropriate densities. Finite numerical samples do
not establish uniform boundedness of the angular coefficient. Channel and
adaptive-grid probabilities must also avoid removing the required density.

## Generic geometry and amplitudes

The framework should be generic from its first implementation. GL638's analytic
ellipsoid and affine-star formulas are optional accelerations, not required
process-specific branches.

* Standard positive-energy surfaces use the active routing matrix, with its
  kernel moved into spectator coordinates, an interior center and numerical
  radial roots on nonempty regular convex fibers.
* Independent blocks use product maps; acyclic dependencies use conditional
  maps with triangular Jacobians.
* Several normals in one block, or nonlinear projected targets, require generic
  implicit charts. For residual vector S, choose normal-coordinate pivots and
  solve S(x_normal,x_tangent)=s. The local Cartesian determinant is
  1/|det(partial S/partial x_normal)|. Multiple patches and explicit branches
  handle pivot changes and global topology.

The last layer is substantial work: root finding alone does not certify support,
invertibility, branch probabilities or projection derivatives. Rank-deficient
and pinched strata need appropriate additional charts. A generic implementation
cannot assume that every nonlinear pullback has a global spherical chart.

For m regular independent normals with an integrable isotropic leading power
R^(-p), p<m, a local density q~R^(-p) can be sampled by
R=Rmax U^(1/(m-p)). Product normal powers with 0<beta_j<1 and sum beta_j>=p
also bound this radial leading power when angular coefficients are bounded.
Inferring the correct power for an arbitrary full integrand is a separate
problem from constructing this generic family of densities.

Amplitudes and cross sections already share the runtime sampling schema and
graph-channel interface. An amplitude supplies its fixed external momenta;
a cross-section target may supply prepared cut kinematics; a star target also
supplies its exact projection/center context. These are contexts for the same
geometry machinery. Automatic and manual modes must both support amplitudes,
including intersecting thresholds in a single amplitude.

## A simpler single channel-selection interface

The following syntax is proposed, not currently implemented:

```toml
[sampling]
sampling_multichanneling = true
sampling_channels = "monte_carlo"
sampling_channel_weight = "map_density"
channel_selection = { GL638 = ["auto:surfaces"] }
```

The requested renames are `lmb_multichanneling` to `sampling_multichanneling`,
`lmb_channels` to `sampling_channels`, and `lmb_channel_weight` to
`sampling_channel_weight`. The second remains the summed/Monte-Carlo execution
mode. Channel selection is a separate list over one named catalogue:

| Selection for a graph | Meaning |
|---|---|
| `["lmb:optimized"]` | Ordinary channels only |
| `["auto:surfaces"]` | Automatically constructed surface channels only |
| `["lmb:optimized", "auto:surfaces"]` | Both families |
| `["H", "HZ"]` | Exactly these user-defined channels |

Defining a channel does not activate it. Explicit ordinary-basis selectors and
named definitions belong to the same catalogue; the old basis-ID filter should
not become a competing selection mechanism. TOML, Python, CLI labels and saved
settings should use the same names and resolved identities.

The new map-density estimator is a semantic change, not an alias for the current
LU inverse-Jacobian partition. It evaluates all channel densities in common raw
coordinates, with branch/support accounting, and multiplies the full physical
integrand. It must retain the existing cut sum and apply adaptive/discrete
sampling weights exactly once. Subtraction groups and their centers do not
change when a sampling channel is selected.

### Surface-only coverage

Surface-only selection must be supported explicitly; the implementation must
not silently append an excluded ordinary LMB channel. Single-surface radial
maps can have normalized support over all r>=0, including ordinary tails.

For example, compactify r with x=r/(r+s), s>0, and let x_star be the image of
r_star. A normalized full-interval profile is

q_x(x) = (1-beta) |x-x_star|^(-beta) /
         [x_star^(1-beta)+(1-x_star)^(1-beta)],  0<x<1, 0<beta<1.

The radial PDF is q_r=q_x s/(r+s)^2. Together with normalized angular and
complement densities, this defines a full-support surface channel. This example
has radial tail r^(-2); tail exponents may need adjustment for a particular UV
behavior. Empty-fiber branches require an explicitly normalized ordinary radial
profile, since there is no shell to focus on. Resolved inspection must show this
behavior. A selected set consisting only of local intersection patches must
provide a declared covering construction or fail coverage validation.

These tail/empty-fiber components are part of the declared channel density,
not additional selected LMB channels. Complete momentum routing is still needed
to represent any of the channels.

### Concise custom definitions

```toml
[sampling.channel_definitions.GL638.HZ]
around = "intersect(cut(2,4,12), cut(3,10,13))"
on_cut = [2,6,10]
```

This names the direct H/Z geometry; exact CT-star images need their own
unambiguous catalogue selectors. `on_cut` supplies a frame, not a restriction
of the physical cut sum. The compiler chooses and exports the complete sampling
parent, signed active cycles, complements, profile and normalization. It must
diagnose ambiguous energy shifts, orientations, host frames or star instances.
Choosing a sampling parent does not change the mandatory supplied parent of any
threshold subtraction metadata.

The shorthand `cut(3,4,7) x cut(2,11) x complement(4,7)` is valid when its factors
are independent. `x` must not silently become an intersection or conditional
map. Experts retain structured controls for exact routing, dependencies, target
frames, energy versus radial widths, powers and branch behavior. Short and long
forms compile into one canonical representation.

## One channel for a physical cut and both sides

Using schematic catalogue labels C, T_L and T_R, propose:

```toml
[sampling.channel_definitions.G.cut_and_thresholds]
around = "phase_space(C) -> (left(T_L) x right(T_R))"
```

The host cut is constructed first. At fixed physical cut data, choose internal
left and right cycles preserving all cut-edge momenta. Their amplitude threshold
maps are then conditionally independent. Their dimensions must be present: a
tree side cannot supply a nontrivial internal loop-threshold block.

In the standard LU rescaling representation, let Q(y,x_L,x_R) denote a point on
the host cut hypersurface C(Q)=0 and restore the raw scale by K=lambda Q. On a
regular radial patch with a unique cut intersection, the cone Jacobian is

d^D K = lambda^(D-1) |Q dot n_C| d lambda d Sigma_C.

Here d Sigma_C is the Euclidean hypersurface measure, itself parameterized by
the cut phase-space and left/right coordinates. Its measure must not be confused
with dQ delta(C(Q)), which additionally contains 1/|grad C|. Existing physical
LU residue factors remain in the integrand; the sampling-map Jacobian is counted
once. For an ordinary cut with normalized positive LU insertion h(t), the
remaining auxiliary radial factor is h(1/lambda)/lambda^2: sample t from h and
set lambda=1/t to match it exactly. The source supplies t_star^(3L) h(t_star)
and the inverse radial cut derivative, which yield this factor after the cone
change. Raised cuts require a positive envelope for derivative-dependent h
terms instead of this simple exact-matching prescription.

The host cut energy residual after its own LU preparation is identically zero.
It is therefore not a third independent physical threshold normal. The physical
cut chart organizes phase space and the raw auxiliary direction; left/right
maps focus genuine threshold distances within that phase space. This is why
`phase_space(C)` is distinct from a generic `surface(C)` shell target.

Other prepared-cut or star targets may couple these blocks and require a
conditional or joint chart. The compiler checks this from the actual routing
and projection derivatives. Selecting a host-cut channel still evaluates the
complete requested cut sum at the resulting raw point.

## Automatic construction and inspection

Automatic mode should build graph- and kinematics-dependent channels for both
amplitudes and cross sections: ordinary shell targets, feasible intersections,
existing exact-star targets, and cut/left/right compositions. It should examine
all requested orientations, deduplicate equivalent geometry, eliminate impossible
intersections using mass/energy bounds and rank checks, and impose a channel
budget. Automatic sampling does not construct or modify the IR subtraction
metadata deferred by the user.

An optional finite pilot can tune widths and probabilities and prune candidates
using variance, extreme weights and evaluation cost while preserving support.
Neither geometric enumeration nor finite pilots prove globally optimal sampling
or globally bounded weights. The selected geometry/profile catalogue should be
frozen before production and recorded for exact replay; ordinary unbiased grid
and channel-probability adaptation may continue.

Proposed `display sampling --candidates`, `display sampling --resolved` and
`inspect sampling --channel NAME` expose generated names, support, dependencies,
inverse densities and full routing. Users may select a generated name, replace
the selection, or export and edit the resolved definitions. The Python API
should accept the same settings structure rather than a second sampling DSL.

Validation must include amplitudes and cut/left/right cross-section channels:
analytic volumes and nonconstant moments; finite-difference Jacobians; forward/
inverse and branch tests; raw-density partitions; summed versus sampled channel
equivalence; support and empty-fiber tests; and unchanged physical cut sums.
GL638 bounded-weight diagnostics additionally need joint H/Z and exact-star
rays, angular/hierarchical scans, all orientations and independent integration
pilots. Unit volume alone cannot verify an arbitrary map.

## Evidence and review

The genericity/bounded-weight argument and API were independently reviewed by
two research agents, including follow-up audits for amplitudes and cut/left/right
composition. Source checks used runtime.rs:1095 and the shared GraphTerm channel
methods in amplitude/mod.rs:1122 and cross_section/mod.rs:1590. The main report
contains the existing numerical experiments and primary-source references,
including Soper's [2001 integration-point paper](https://arxiv.org/pdf/hep-ph/0103262)
and Ohl's [multichannel adaptation paper](https://arxiv.org/pdf/hep-ph/9806432).
Implicit chart atlases and the proposed API are our design proposals, not claims
that those papers provide an implementation of this complete framework.
