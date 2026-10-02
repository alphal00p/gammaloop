# LUDY calculations with FeynKit's symbolic stack

An executable **partial port through scalar integrands and integral-index
reduction**, with the full-integration gaps recorded in
[investigation.typ](investigation.typ).

```sh
python examples/ludy/compute.py --check --output /tmp/ludy.json
```

Use a current Symbolica host with FeynKit, Spenso and Idenso. Hosts that package
FeynKit as `symbolica.community.hep` need `--module hep`.
The private reference checkout is not required to run the port.
