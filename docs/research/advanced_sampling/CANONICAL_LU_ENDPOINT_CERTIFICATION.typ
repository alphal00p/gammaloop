= Canonical LU root certification at an exhausted bracket
<canonical-lu-root-certification-at-an-exhausted-bracket>
The exact GL638 baseline source that interrupted confirmation now completes ordinary Double and forced-Arb physical evaluation. Both return six valid cuts with the original canonical point, Jacobian, partition and complete Sample unchanged. The existing native accuracy audit passes all 25 checks. The #link("GL638_CONFIRMATION_ROOT.typ")[original interruption] remains preserved and unpooled; a fresh complete confirmation is still required.

The captured source reaches canonical physical overlap preparation for cut group 3, physical cut 0, before body evaluation. Its Arb1000 root bracket has no representable interior midpoint. The last/lower endpoint fails the existing bracket-resolution allowance, while the other already represented endpoint satisfies it:

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Quantity], [Last/lower endpoint], [Alternate/upper endpoint],),
    table.hline(),
    [Native residual], [−2.866985836e−298], [+1.911323891e−298],
    [Original equation residual, directed enclosure], [≈−9.206668058e−299], [≈+9.549216133e−301],
    [Native strict budget], [1.866527237e−298], [1.866527237e−298],
    [Existing bracket-resolution budget], [2.796743259e−298], [2.796743259e−298],
    [Native resolution acceptance], [Fails], [Passes],
    [Original equation within strict budget], [Passes], [Passes],
  )]
  , kind: table
  )

The directed interval widths are approximately 1.15e−613 and 1.01e−613; both exclude zero. They establish finite residual accuracy, not an exact root. The LU settings stay at tolerance 1, Ecm 1000 GeV, 2000 iterations and 64 bracket expansions. This owner already used the routed ray and diagnostic wrapper.

`RadialRootDiagnostics::solve` now checks the alternate endpoint after both existing accepted return paths. It uses the same residual allowance `epsilon*tolerance*Ecm + abs(derivative)*bracket_width`, with finite values, positive radius/derivative, a valid exhausted bracket and fresh four-probe consistency checks for that candidate. Acceptance records its own observation; rejection retains the original error and precision history. Sampling laws and physical expressions remain unchanged. Each consistency check uses four probes; checking the alternate can add another set.

The regression pins the complete old strict error and original-equation bounds. An adversarial callback changes only the alternate derivative by a factor of eight: its residual fits the allowance, but its own slope checks reject it. Thirty-five broader tests pass in 470.160 s. The final endpoint fixture passes in 1.244 s after two needless-borrow cleanups and a stronger iteration assertion; it overlaps one of those 35 tests. Production solver source is identical across both runs. Final format, check and clippy pass, with 52 pre-existing diagnostics and none on changed lines.

Optimized build7 passes in 443.487 s with unchanged source. Build6 was deliberately interrupted after 38.638 s for the test-only cleanup; it supplies no successful build or physical evidence. The retained source identity is base `3e68f16588b8272c0f88513c0c43305b5d67a059` plus patch `3f5669ae1de62c09b26d4aa3888cecd47058a786f00e0893bac7bfd191910588`.

The rebuilt one-worker replay takes 74.056 s including loading. Both calls are `Stable(2)` and leave all 35 state hashes unchanged. Their complete-sample total has relative complex-norm difference 1.36411e−13, below the unchanged combined Double/Arb budget 1.000001e−6; all six cut comparisons pass. The one own-component relative miss is cut-2 Re: Double zero versus Arb 3.87479e−565. Its negligible absolute error passes the existing total-scale cut criterion. It is retained as an underflow diagnostic and does not make the ordinary lane unstable.

The #link("gl638_hosted_joint_gate/mc_confirmation_root/endpoint_repair/artifact_hashes.json")[50-record archive] contains the original endpoint log/source, audited Decimal700 interpretation, source patch and both fixture snapshots, exact gate commands, interrupted and completed build provenance, complete replay records and unchanged metric owner. The endpoint parsing snippet was not separately retained; its raw log and fixture preserve the reproduction inputs. See the #link("gl638_hosted_joint_gate/mc_confirmation_root/endpoint_repair/source_audit.json")[source audit];, #link("gl638_hosted_joint_gate/mc_confirmation_root/endpoint_repair/gates/summary.json")[gate records];, #link("gl638_hosted_joint_gate/mc_confirmation_root/endpoint_repair/replay/accuracy_audit.json.gz")[physical accuracy audit] and #link("gl638_hosted_joint_gate/mc_confirmation_root/endpoint_repair/provenance.json")[current provenance];. Earlier source-only audit snapshots still say build/physics pending; the later records and current provenance establish completion. No generated-process state payload or executable is copied, and no confirmation statistics are inferred from this single-source repair gate.
