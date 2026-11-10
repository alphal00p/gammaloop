= Immutable fixture manifest

Archive: `fixtures.tar` (USTAR, sorted groups: cases, references, licenses;
regular files only, mode 0644, uid/gid/mtime zero, empty owner names).
SHA256: `8f5cf9c5ba70d466f69b800882552aa0b41171efd96cd8fe73b05fb585a99e1f`.
The archive contains only 24 inputs, two supplied SVG hash lists and four
license files, not frozen sources or packages. `run.py --dry-run` performs
deterministic extraction and writes readable per-file hashes to metadata.

Reconciled baseline archive: `references.tar` (USTAR, regular files only,
mode 0644, uid/gid/mtime zero). SHA256:
`633b5b736b31d9f1378cf9b4108c3785a036e0d874558b90b469297b3be3142c`.
It contains 24 PNGs at 150 ppi and four generated reports: source/environment
metadata, raw timings, medians, and paired snapshot PNG hashes.

The following diagram groups each have `-momenta.typ` and `-particles.typ`,
except gluonic/fermionic whose second suffix is `-notebook.typ`:

- gluonic: four-loop numerator notebook FK0, 4 loops, 8 vertices/13 edges,
  2 external legs, 13 gluons.
- fermionic: same notebook, 4 loops, 8/13, 2 legs, 8 b-quarks/5 gluons.
- paper-b1: double triangle soft IR, 3 loops, 6/10, 2 legs,
  6 d-quarks/2 gluons/2 photons.
- higgs-sunset: 2 loops, 2/5, 2 legs, Higgs.
- crossed-higgs-double-box: 2 loops, 6/11, 4 legs, Higgs.
- nonplanar-higgs-k33: non-planar K3,3, 4 loops, 6/13, 4 legs, Higgs.
- multi-leg-six-gluon-incoming: 1 loop, 6/12, 6 legs, gluons.
- multi-leg-six-higgs-theta: 2 loops, 8/15, 6 legs, Higgs.
- multi-leg-eight-gluon-ladder: 3 loops, 8/18, 8 legs, gluons.
- xs-extra-higgs-scattering: cross section, 2 loops, 4/7, 2 initial states.
- xs-extra-higgs-radiative: cross section, 3 loops, 4/7, 1 initial state.
- xs-extra-higgs-two-loop-virtual: cross section, 4 loops, 6/10, 1 initial state.

All except the notebook pair originate in Linnest native layout
`crates/linnest/layout/native/fixtures/seeds.json`. Notebook source:
symbolica-community `examples/hep/four_loop_numerator.py`.
Cross sections split initial-state particles into incoming/outgoing legs.

== Exact input SHA256
```text
42e82c8ccaa2618e767293aa3e9876b5e47ca8e13cebe1adf633bbb00843b86c crossed-higgs-double-box-momenta.typ
bc58b3f5925b7977b78380fed60e61f338ca779445dc0d45e9866bf3cb44ee9d crossed-higgs-double-box-particles.typ
8c31d5a1db3596fc30b12af9c87e18a8277ff4c9fd7b409e3e813a1117fa1a65 fermionic-momenta.typ
081e4b62b35b77b79a3f27fa8df1d511e2fdea27cf2e5cc023e6e25998e355e2 fermionic-notebook.typ
eddbb12c9162c1c06257d83ec3d4e7b3c9b62e4efa3f3568232dede55c600d8c gluonic-momenta.typ
a4b40a372213de06cc2f734ea11f3b1b730db0ecf0775134e1746717d4c54814 gluonic-notebook.typ
e9563ee88992b29798b6d28303c305745d551b226e51e99ddb668b92bedbbcf8 higgs-sunset-momenta.typ
50168926bb66b2ebca1a0ba2b82736ca3a1edc4957be793a61eaeeb696b215dc higgs-sunset-particles.typ
2b45222ae49777ba139c1f3316cf537b5d3ed55ba77318f1d2d198acf04e17e3 multi-leg-eight-gluon-ladder-momenta.typ
fe7e306e88a9226bd0f0dab18aead9c1e4652510c0f1d4d24c8b09b5cd71961d multi-leg-eight-gluon-ladder-particles.typ
1fb3755c1034663e0e2540403a5a44934aca5b1d4c56a645a8a11d1749189444 multi-leg-six-gluon-incoming-momenta.typ
8bdc0c4fe7cdbce5e201d56a6506c19d1edcb4f417f23db0f048608a3dcea3cf multi-leg-six-gluon-incoming-particles.typ
9579dbe89ae0446bddc323a924793e35948509c2043f5a67442dce82d9f52379 multi-leg-six-higgs-theta-momenta.typ
75944645db39598a993a145a976ca4ae38e3d3c305c6f7b1189a65a5d543f9bc multi-leg-six-higgs-theta-particles.typ
a31b12dd45816d4ecbce4b2904deaeec4a9cf2f0ea87d789330f3fdb72bf8d33 nonplanar-higgs-k33-momenta.typ
3f9d00631e5a73ed7bf55b5fa7e66817ad3f8b24b37142c4a9934f6b816cba6a nonplanar-higgs-k33-particles.typ
e0d3a09c1562a29ff051d9ed124290e48142dbfdb641f019d7c5be0aa7004ac4 paper-b1-momenta.typ
c03f515e76b8e01afdfc764c4c26d5e2dcdf5353119725b30e12bd20c2013f8d paper-b1-particles.typ
3bc40cf01cec6f0beda139fd66c3e4caa39aa5626ab72ee3455c7abdd73e7e03 xs-extra-higgs-radiative-momenta.typ
b97e4010af4711af60bccdd7ccf3e9992790093e3e30082385bada89e8479752 xs-extra-higgs-radiative-particles.typ
8d59ee3a60ce73e3b90c717b8e9fe4021d1af42e2209ecbca0872ce4090ac57b xs-extra-higgs-scattering-momenta.typ
f944df1dbbacd7993e13f1757547d893c42eb6156aa3ec83ada8215389abc664 xs-extra-higgs-scattering-particles.typ
ea7a0b6040846593eedcde49f99c624ae116232082d77f1cf8a58c01e65c6082 xs-extra-higgs-two-loop-virtual-momenta.typ
9e3912b74da08127614e9408d7bcf1495590f433ec189ec370076eb86dc6425d xs-extra-higgs-two-loop-virtual-particles.typ
```
