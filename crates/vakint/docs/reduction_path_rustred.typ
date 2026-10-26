#set document(title: "The RustRed reduction path", author: "ValentinHirschi")
#set page(paper: "a4", margin: (x: 17mm, top: 16mm, bottom: 16mm),
  footer: context align(right, text(size: 8pt, fill: rgb("687584"))[RustRed reduction path · #counter(page).display("1")]))
#set text(font: "Libertinus Serif", size: 10.4pt, fill: rgb("243343"))
#set par(justify: false, leading: 0.62em, spacing: 0.75em)
#set heading(numbering: none)
#show heading.where(level: 1): set text(size: 24pt, weight: "bold", fill: rgb("173f65"))
#show heading.where(level: 2): set text(size: 17pt, weight: "bold", fill: rgb("173f65"))
#show heading.where(level: 3): set text(size: 11.6pt, weight: "bold", fill: rgb("173f65"))
#let blue = rgb("eaf2fa")
#let green = rgb("eaf5ef")
#let amber = rgb("fff3df")
#let gray = rgb("f1f4f7")
#let ink = rgb("173f65")
#let panel(body, fill: blue) = block(width: 100%, inset: 9pt, radius: 4pt, fill: fill, body)
#let small(body) = text(size: 8.1pt, fill: rgb("5b6876"), body)
#let stage(n, title, body) = grid(columns: (9mm, 1fr), gutter: 3mm,
  align(center, text(size: 18pt, weight: "bold", fill: ink, str(n))),
  [*#title* #linebreak() #body])

= The RustRed reduction path

#text(size: 11.8pt)[Where equal integrals meet, what “fully reduced” means, and what a cache remembers.]

#small[Source snapshot: Native `4190f2c6` and Vakint `d29eca063` · 26 October 2026.
This describes the existing default path. The examples are schematic exact identities, not measured reductions.]

== There are intermediate merges

RustRed does not postpone every addition until the end. It merges equal keys at
several boundaries, but each boundary sees a different object. The distinction is
between *equal immediate children* and *equal descendants hidden inside separate
branches or roots*.

#panel[
  #stage(1, [Admit and collect the input roots], [
    Validate every incoming key and coefficient, retain its conditions, then add
    coefficients of identical roots. Exact zero entries disappear.
    Thus $a I + b I arrow.r (a+b) I$ before recursion, while admission and denominator
    checks still apply even if $a+b=0$.
  ])
  #align(center, text(fill: ink)[↓])
  #stage(2, [Evaluate one rule; collect its immediate children], [
    Specialize the selected formula, check its original denominators and descent,
    then combine equal RHS child keys. A local RHS $U-U+V$ becomes $V$;
    that occurrence of $U$ needs no recursive expansion.
  ])
  #align(center, text(fill: ink)[↓])
  #stage(3, [Finish children by depth-first traversal; combine their results], [
    A child is solved or supplied by a completed cache entry. Once all surviving
    children are ready, back-substitute their terminal maps, multiply by the edge
    coefficients, and add equal terminal keys progressively. Exact zeros disappear
    at this level too. Store the completed parent result.
  ])
  #align(center, text(fill: ink)[↓])
  #stage(4, [Collect one completed root into the request output], [
    Multiply that root's terminal map by its external input coefficient, then add
    it to the running output. Continue with the next root. Cross-root cancellation
    happens here, after each contributing root has been terminal-expanded.
  ])
]

=== What “fully reduces a root” means

For a successful strict reduction, all required descendants have been replaced by
declared Native stopping terminals or sealed zeros, with exact parameter-dependent
coefficients and the required nonzero conditions. A completed cached child can
replace its entire recursive computation. “Fully” describes the *result*, not an
instruction to recompute every descendant, and not evaluation of a Laurent series
or a numerical master integral.
Coefficients may retain a factorized internal representation: terminal completion
does not require expanding every rational-polynomial product.

The parent combination waits for its children to finish. It does not continuously
pool the not-yet-reduced descendants of every branch. Optional preserve mode can
leave undeclared integrals explicit; a ready partial result is not strict terminal
completion and is not eligible as a distributed completed identity.

#small[*Source anchors.* Native `core/solver/candidate_reduction/evaluator.rs:250–302`;
`core/solver/candidate_reduction/routed/application/engine.rs:877–1000,1070–1124,1387–1417`.
Here and below `core/` means `crates/rustred-core/src/`; `app/` means `crates/rustred-app/src/`.]

#pagebreak()
== A cancellation can be real and still arrive late

Let $M_1, dots, M_4$ denote four declared stopping objects in this example. Assume
the selected exact rules are

#panel[
  #grid(columns: (1fr, 1fr), column-gutter: 10mm, row-gutter: 5pt,
    [$X arrow.r 2 U + V$], [$Y arrow.r -2 U + W$],
    [$U arrow.r M_1 + M_2$], [$V arrow.r M_1 - M_3$],
    [$W arrow.r -M_1 + M_4$], [input: $X+Y$])
]

=== Current traversal: completed child results first

#table(columns: (21mm, 1fr, 1fr), inset: 7pt, stroke: 0.4pt + rgb("d5dee7"),
  table.header([*Step*], [*Work / available value*], [*Collected result*]),
  fill: (x, y) => if y == 0 { blue } else if calc.odd(y) { gray } else { white },
  [Solve $X$], [Finish $U$ and $V$; combine $2 U+V$.],
    [$3 M_1 + 2 M_2 - M_3$],
  [Solve $Y$], [Reuse completed $U$ if still resident; finish $W$; combine $-2 U+W$.],
    [$-3 M_1 - 2 M_2 + M_4$],
  [Collect], [Add the two fully expanded root contributions.],
    [$-M_3+M_4$],
)

The coefficients of $M_1$ and $M_2$ cancel exactly. Nevertheless $U$ was needed
while computing $X$. Under the stated residency assumption, it is computed *once*
and later reused: cancellation of its final contribution does not mean its
expansion was avoided. With finite eviction, a later miss may recompute it;
an admitted persistent hit can also supply the completed value.

#grid(columns: (1fr, 1fr), gutter: 4mm,
  panel([
    *Local cancellation — implemented* #linebreak()
    $R arrow.r U-U+V arrow.r V$ #linebreak()
    The identical child keys meet in one rule's RHS collection. No recursive
    work for these cancelled occurrences of $U$ is scheduled. Checks attached to
    the selected formula are still retained.
  ], fill: green),
  panel([
    *A global weighted frontier — hypothetical* #linebreak()
    $X+Y arrow.r 2 U+V-2 U+W$ #linebreak()
    Pool these descendants first: $2U-2U=0$, then solve only $V+W$.
    The current default traversal does not perform this pooling.
  ], fill: amber),
)

The hypothetical path would avoid $U$'s work, but would then have no completed
unit mapping for $U$, $X$ or $Y$ to cache without additional work. Actually
completed $V$ and $W$ mappings could be cached. The current path instead produces
reusable mappings for $U$, $X$ and $Y$; a weighted partial/global sum must never
be stored as one of those individual identities.

=== A diamond within one root has the same distinction

If $R arrow.r A+B$, $A arrow.r U+V$ and $B arrow.r -U+W$, $A$ and $B$ are
different immediate children. RustRed finishes their results before combining
them at $R$. Shared $U$ can be reused from cache, but the $+U$ and $-U$ hidden
inside different branches are not first merged into a common unreduced frontier.

=== The notation hides the operational keys

Our $I$, $U$ and $M_i$ suppress family and routing detail. Native actually memoizes
*Route(key)* and *Apply(owner, key)* separately. A same-owner child remains Apply;
a proper pinched support can return to Route for transport to its owner. An equal
vector of integral powers alone does not make those two states interchangeable.
Every stored expansion is for the *unit node*: the external coefficient $a$ of
$a I$ belongs to the calling edge/root, not to the cached identity for $I$.

#small[*Source anchors.* Native `core/solver/candidate_reduction/cache.rs:149–179`;
`core/solver/candidate_reduction/routed/trace.rs:12–16,89–108`;
`core/solver/candidate_reduction/routed/application/engine.rs:983–1000,1085–1124`.
The example assumes no eviction of the completed shared node between its uses.]

#pagebreak()
== What is remembered, and when it becomes a file

#table(columns: (30mm, 1fr, 1fr), inset: 7pt, stroke: 0.4pt + rgb("d5dee7"),
  table.header([*Layer*], [*How it is generated*], [*What it means*]),
  fill: (x, y) => if y == 0 { blue } else if calc.odd(y) { gray } else { white },
  [*In-memory* #linebreak() native memoization],
    [DFS creates completed terminal maps for leaves and combined parents;
     later roots can reuse completed descendants.],
    [Exact unit-node expansion and its guards. Bounded policies may evict
     completed entries; live parent dependencies are pinned.],
  [*Persistent* #linebreak() backing / capture],
    [Preparation/export or an attached capture sink validates eligible completed
     records. A backing supplies previously generated records on demand.],
    [Exact unit-node coefficients and guards, bound to operational keys and the
     authenticated package. No pending stack or weighted partial request sum.],
  [*Running* #linebreak() physical-example cache],
    [When attached, normal reductions forward completed descendants to a bounded
     writer while the larger request runs. Backing reads and captures are separate.],
    [A durable prefix can survive an unfinished root/workload. Flush/snapshot/close
     establish I/O barriers; a closed writer does not certify workload completion.],
  [*Publication* #linebreak() compressed cache],
    [A separate export/bundle/archive step packages selected or cumulative
     records for transport.],
    [Compression and selection create no new identities. Hosting the artifact
     is separate from producing it and validating the workload.],
)

An entry means $I = sum_i c_i M_i$ under retained conditions $g_j != 0$.
The persistent binding covers the source package, routes, entry scope, zero
authority and terminal substitutions, plus family, arity and solver layer.
Loads check the binding and current policy; coefficients and conditions are kept
together. A coefficient cancelling to zero does not erase the original guards
of an evaluated formula or imported subtree.

*Unresolved leaves are excluded from distributed completed records.* Preserve
mode can memoize a partial result locally. Capture/export checks its remaining
keys and skips it if any is undeclared to the bound Native application. An exactly cancelled parent can itself be
complete and exportable even when a partial child's identity is excluded.

=== From Vakint's numerator to Native terminals and the master catalogue

#panel[
  #align(center)[
    physical numerator → scalar lowering → weighted scalar keys #linebreak()
    ↓ #linebreak()
    Native Route/Apply → terminal map → saved substitutions #linebreak()
    ↓ #linebreak()
    optional terminal-product tensor projection / lower-factor reduction #linebreak()
    ↓ #linebreak()
    exact surviving symbols → external catalogue → retain / substitute
  ]
]

Scalar lowering creates weighted scalar keys before Native reduction. The
optional `factor_projection` callback is later: it projects tensor numerators
in terminal product charts, retaining projector denominator conditions *before*
materializing and simplifying their sum. Native stopping terminals and saved
substitutions are internal conventions; their final family-qualified output map
is not the external catalogue itself.

After exact combination of the resulting master symbols, Vakint distinguishes
known coefficients, declared-but-uncomputed masters, bare Native announcements,
and undeclared residuals. The 640-entry five-loop campaign declares external
masters with uncomputed coefficients; it is not a claim that Native has 640
leaves or that those masters have been evaluated. This classification precedes
optional master-coefficient substitution. Vakint's topology matching at ingress
is a different use of a catalogue.

#small[*Source anchors.* Native `core/solver/candidate_reduction/routed/application/persistent.rs:71–101,291–344,401–460,504–546`;
`app/application/saved_application/mod.rs:433–476,532–607`; `app/application/saved_application/running_cache.rs:435–456,689–718`.
Vakint `src/rustred_evaluation/factor_projection.rs:90–113`;
`cache_five_loop.rs:299–419`; `cache_running.rs:20–33,234–256`;
`five_loop.rs:488–527,889–967,1006–1021`; `masters.rs:513–563,644–655`;
`data.rs:1149–1178`.
Revisions: Native `4190f2c6f3158e3fc3cbcb1937e12e5971da9175`;
Vakint `d29eca063a33690c8395a8e0c34dca869847f987`.]

#pagebreak()
== Integral descent is stronger than template ordering

RustRed checks a hierarchy of *instantiated integral keys*. Every surviving child
of an applied formula must be strictly smaller than its target under the saved
owner's mathematical ordering. Between owners, leaving the current support
requires a proper pinch; verified transports preserve the number of active
denominators. The detailed local order can include physical weighted degrees,
numerator and dot excess, and coordinate tie breakers. *Rank alone is not a
general progress measure.*

These checks give a conditional termination argument for finite concrete inputs
and finite branching, assuming each local algebra/transport operation finishes.
They do not themselves establish that every nonmaster has an applicable rule,
or that a single sweep through rule templates suffices.

=== Why a fixed exhaustive template sweep can fail

Consider the following schematic rewrite system, unrelated to the actual IBPs:

#panel[
  #grid(columns: (1fr, 1fr), column-gutter: 8mm,
    [$A_n arrow.r B_(n-1)$ for $n >= 1$],
    [$B_n arrow.r A_n$ for $n >= 0$])
  #align(center)[$A_0$ is terminal; $mu(A_n)=2n$, $mu(B_n)=2n+1$.]
]

Every step strictly decreases $mu$, yet full reduction requires alternating
the two templates:

$ A_2 arrow.r B_1 arrow.r A_1 arrow.r B_0 arrow.r A_0. $

Exhausting the first template and then the second leaves $A_1$ reducible.
Reversing the template order leaves $B_1$ reducible. The graph of concrete
instances is acyclic; the template graph has a cycle.

One exhaustive sweep needs an additional triangularity property: *later templates
must never create matches for earlier templates*. A dependency analysis could
identify ordered phases and mutually dependent groups requiring local iteration.
No such one-sweep template ordering has been established here for the campaign.

=== The useful order is over the live weighted integrals

Suppose routing and canonicalization have been incorporated into a common
well-founded order, and each permitted rewrite satisfies

$ K arrow.r sum_j c_j K_j, quad K_j < K. $

A global frontier can then select its *greatest live key*, apply one valid rule,
and immediately merge the children with all other live contributions. A key is
not expanded while a larger parent capable of contributing to it remains.
Consequently, each distinct key needs expansion at most once per request.

#panel(fill: green)[
  *Cancellation before work.* In the example on page 2, both $X$ and $Y$ are
  greater than $U$. Process them before $U$: its coefficient becomes
  $2-2=0$, so its complete descendant calculation is avoided. Exhausting one
  template recursively could reach $U$ before the other parent contributes.
]

Exact combination requires canonical identities and compatible validity
conditions. Branching can increase the number of live terms even while every
child is smaller; the active collection's multiset order supplies the descent
argument. A common order across owners is an implementation obligation, not
something provided by sorting only on rank.

Completed unit memoization already avoids much repeated work. A frontier adds
the opportunity to skip globally cancelled branches and their large terminal
maps, at the possible cost of a wider live expression. Contextual cancellation
still supplies no reusable unit identity for the cancelled integral. Neither
descent nor the frontier alone proves completeness of the master catalogue.

#small[*Hierarchy audit.* Native `core/solver/candidate_reduction/evaluator.rs:161–189,266–301`;
`core/solver/candidate_reduction/routed/trace.rs:89–109`;
`core/sector/ordering/comparison.rs:44–111`;
`docs/research/exact_finite_traversal_termination_2026-10-15.md:7–76`.
This page states the conditional frontier argument; it does not describe an
implemented global frontier or certify the campaign's template dependency graph.]
