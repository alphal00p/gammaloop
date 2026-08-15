# Citing GammaLoop workspace software

Zenodo archives this repository as six independent software families. Each
family has one stable concept DOI for the project as a whole and a distinct
version DOI for every archived release. Cite the version DOI when reporting a
calculation so that the cited source and files are immutable. Use the concept
DOI only when referring to a project without selecting a version.

## Citable families

| Family | Current automated archive scope | Version and tag |
| --- | --- | --- |
| `gammaloop monorepo` | Complete tracked repository | GammaLoop release version, `vX.Y.Z` |
| `gammaloop` | `gammalooprs`, `gammaloop-api`, and the Python source distribution | `gammaloop-api` version, `vX.Y.Z` |
| `spenso` | `spenso`, `spenso-macros`, `spenso-hep-lib`, and `spynso3` | `spenso` version, `spenso-vX.Y.Z` |
| `idenso` | `idenso` | `idenso` version, `idenso-vX.Y.Z` |
| `vakint` | `vakint` | `vakint` version, `vakint-vX.Y.Z` |
| `linnet` | `linnet`, `clinnet`, `linnet-py`, and `linnest` | `linnet` version, `linnet-vX.Y.Z` |

The Linnet family continues the existing
[concept DOI `10.5281/zenodo.15494393`](https://doi.org/10.5281/zenodo.15494393).
Its adopted `0.17.0` record keeps the creators and CC-BY-4.0 metadata with which
it was originally published; the migration does not rewrite that immutable
history. It also predates this repository's shared grant and community policy,
so that baseline record does not acknowledge the grant or belong to the
`alphaloop` community. New monorepo-managed Linnet versions use Lucien Huber as
creator, the repository's MIT license, the shared grant, and the community.
The other concept DOIs are recorded in
[`records.json`](../.zenodo/records.json) after their first production deposits.
Until then, their citation files intentionally omit a DOI rather than
publishing a placeholder identifier.

Spenso and Idenso start independent release series at versions `0.6.0` and
`0.3.0`. Both are derived from the frozen
[combined legacy concept DOI `10.5281/zenodo.15913113`](https://doi.org/10.5281/zenodo.15913113).

## Historical archives

The initial archive migration backfills:

- full-monorepo versions `0.2.0`, `0.2.1`, `0.3.0`, `0.3.3`, and `0.3.4` from
  their Git tags;
- GammaLoop versions `0.2.0`, `0.2.1`, `0.3.0`, `0.3.3`, and `0.3.4` from Git
  tags, plus PyPI-only versions `0.0.1`, `0.1.0`, `0.1.1`, `0.3.1`, and `0.3.2`,
  which have no matching repository tag;
- Vakint versions `0.1.0` and `0.1.2`.

Earlier Spenso and Idenso releases remain represented by their combined legacy
record. Linnet continues its existing Zenodo version history instead of
creating a replacement concept. Historical GammaLoop and monorepo records
preserve the layouts and source artifacts that existed at their respective
versions; they do not claim the current workspace split or artifact set.

## Archived files and relations

Each synchronized automated family version contains a deterministic source
archive, its applicable canonical `.crate` archives, a release-specific
`CITATION.cff`, the MIT license, `PROVENANCE.json`, and `SHA256SUMS`. The
GammaLoop family also contains the Python source distribution. Historical
backfill records instead contain the source artifacts enumerated for that
release. Python wheels remain available from PyPI and the GitHub release but
are deliberately excluded from Zenodo, keeping the citable archive
platform-independent and source-based.

The monorepo version links to the exact component version DOIs with `hasPart`.
Each component version links back to the monorepo concept DOI with `isPartOf`.
Automation reuses an existing version only when its immutable files and stable
metadata match, so retrying a release does not create a duplicate citation.

The one-time historical migration is deliberately one-directional: backfilled
component records link to the monorepo concept, but already-published historical
monorepo records are not rewritten later to add component DOIs. Exact `hasPart`
links begin with the synchronized automated releases.

## Citation files

GitHub displays the root [`CITATION.cff`](../CITATION.cff) for the GammaLoop
product. The monorepo and component citation files are:

- [gammaloop monorepo](../.zenodo/citation/gammaloop-monorepo.cff)
- [spenso](../crates/spenso/CITATION.cff)
- [idenso](../crates/idenso/CITATION.cff)
- [vakint](../crates/vakint/CITATION.cff)
- [linnet](../crates/linnet/CITATION.cff)

The checked-in files describe the current software families without embedding
a mutable release version. During a release, automation generates an archived
copy containing the exact version DOI, concept DOI, version, and release date.
Zenodo also provides BibTeX and other formatted citations on every record page.

## License and Symbolica

Code authored in this repository is available under the MIT License. Modern
official GammaLoop Python distributions, normal Cargo builds from the repository
root, and the default Nix package activate Symbolica under GammaLoop's OEM
license. Consumers compiling modern `gammalooprs` or `gammaloop-api` releases
from crates.io must set the public compile-time selector
`SYMBOLICA_OEM_LICENSE=SYMBOLICA_OEM_GAMMALOOP`; it is not a personal license or
secret. The OEM grant does not transfer to standalone libraries: users of
Spenso with Symbolica integration, Idenso, or Vakint must configure their own
Symbolica license. Linnet's default functionality does not require Symbolica,
but its optional `symbolica` feature requires a user-provided license.
Historical versions retain their original Symbolica integration and activation
behavior; consult the archived source for the cited version. Symbolica and
other third-party software remain governed by their own terms.

New records for all six families acknowledge Swiss National Science Foundation
grant `PCEFP2_203335` and are submitted to the `alphaloop` Zenodo community. The
adopted Linnet `0.17.0` baseline is the documented historical exception.
