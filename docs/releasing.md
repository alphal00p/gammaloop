# Releasing GammaLoop

Releases use a reviewed version/changelog pull request followed by one protected
workflow. The workflow publishes the public Rust crates in dependency order,
publishes the `gammaloop` Python distribution, creates annotated component tags
protected by a repository ruleset, creates the immutable GitHub release, and
archives the six citable software families on Zenodo.

## One-time repository setup

1. Allow GitHub Actions to create pull requests. If pull requests created with
   `GITHUB_TOKEN` do not start the required pull-request checks, give the
   release preparation workflow a narrowly scoped GitHub App token instead.
2. In the repository's **Settings** page, enable **Release immutability** in
   the **Releases** section before publishing the first automated release. It
   protects the canonical assets and associated tag consumed by PyPI and
   Zenodo, and applies only to releases created after it is enabled. See
   [GitHub's immutable release settings](https://docs.github.com/en/code-security/how-tos/secure-your-supply-chain/establish-provenance-and-integrity/prevent-release-changes).
3. Under **Settings > Rules > Rulesets**, create an active tag ruleset for
   `*-v*` with **Restrict updates** and **Restrict deletions** enabled. Leave
   **Restrict creations** disabled so the release workflow can create each new
   annotated component tag. GitHub release immutability protects only the
   canonical `vX.Y.Z` tag associated with the release; this
   [tag ruleset](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/creating-rulesets-for-a-repository)
   protects the component tags.
4. Create the `release-tag`, `crates-io`, `pypi`, `github-release`, `zenodo`, and
   `zenodo-sandbox` environments. Restrict them to the protected release path.
   Required reviewers add a deliberate manual gate; omit them from publishing
   environments when releases should run unattended after the tag is created.
5. Configure crates.io Trusted Publishing for every existing public crate,
   using `.github/workflows/release.yml` and the `crates-io` environment.
6. Configure PyPI Trusted Publishing for the `gammaloop` project, using
   `.github/workflows/release.yml` and the `pypi` environment.
7. Configure the repository's `SYMBOLICA_LICENSE` secret with a maintainer-owned
   key for standalone `spenso`, `idenso`, and `vakint` CI. The secret is not
   passed to official GammaLoop distribution builds.
8. Confirm that the GammaLoop OEM agreement covers every published platform and
   artifact type. Official Python wheels and source distributions and the
   default Nix package activate the GammaLoop OEM path without an end-user
   `SYMBOLICA_LICENSE`. The credential-free release gate also builds an external
   `gammalooprs` and `gammaloop-api` consumer with the public compile-time selector
   `SYMBOLICA_OEM_LICENSE=SYMBOLICA_OEM_GAMMALOOP`. Standalone
   Symbolica-dependent libraries do not inherit this entitlement and require
   their users to configure their own license.
9. Create separate Zenodo personal access tokens on
   [production](https://zenodo.org/) and
   [sandbox](https://sandbox.zenodo.org/). Give each token only the
   `deposit:write` and `deposit:actions` scopes documented by the
   [Zenodo deposition API](https://developers.zenodo.org/). Store each as an
   environment secret named `ZENODO_TOKEN`: the production token in `zenodo`
   and the sandbox token in `zenodo-sandbox`. Do not reuse either token on the
   other service or expose it as a repository-wide secret.
10. In the `alphaloop` community settings on both services, open **Submission
   policy** and select **Allow curators, managers and owners to publish without
   review**. The token owner must be an owner of both communities. This
   [Zenodo review policy](https://help.zenodo.org/docs/communities/manage-community-settings/review-policy/)
   lets the workflow publish into the community without a separate UI approval.

## Zenodo sandbox and historical migration

Run the `Zenodo sandbox` workflow from `main` before the first production
deposit and after changing Zenodo metadata or automation. It packages the
current sources, publishes run-suffixed test versions on Zenodo Sandbox, and
audits the resulting files, metadata, DOI relations, and community membership.
Sandbox and production are separate services, so a successful sandbox run does
not create or change a production record.

After a successful sandbox run, run the production `Zenodo backfill` workflow
once from a reviewed `main` commit. The workflow follows the immutable sources
and checksums in [`.zenodo/backfill.json`](../.zenodo/backfill.json), creates the
historical versions in chronological order, initializes independent Spenso and
Idenso lineages, and adopts the existing Linnet lineage. Check the five new
family pages in the production `alphaloop` community and inspect Linnet's
existing concept page separately after it finishes. The adopted Linnet `0.17.0`
baseline predates the shared grant/community policy and remains outside that
community; future Linnet versions use both. Normal releases refuse to publish
until every production baseline in
[`.zenodo/records.json`](../.zenodo/records.json) exists, so the backfill must
finish before the first automated release.

Download the backfill result artifact and copy the five newly assigned concept
record IDs into `.zenodo/records.json` in a follow-up pull request. Add the same
concept DOIs to the corresponding checked-in citation files. Linnet is already
pinned to concept record `15494393`. The stable family markers let retries find
the new lineages before that follow-up is merged, while recording the IDs makes
future lineage checks and GitHub's citation UI more explicit.

## Normal release

1. Review the draft release pull request opened by `Prepare release`. Confirm
   its changelogs and versions, and make it ready for review. The
   `gammaloop-api` and `gammalooprs` versions must be identical. Because the
   Python version is derived from `gammaloop-api`, there is no separate Python
   version to edit.
2. Merge the release pull request after CI passes.
3. Run the `Release` workflow on `main` with the new version, without the `v`
   prefix. After the registry-credential-free Rust and Python gates pass, the
   workflow creates the annotated `vX.Y.Z` tag, publishes the crates.io
   packages, preserves the gated artifacts in the package tags and GitHub
   release, publishes those canonical Python files to PyPI, and publishes or
   reuses the matching Zenodo family versions.
4. Approve protected environments if required by repository policy.

An annotated `vX.Y.Z` tag pushed by a maintainer also starts the same workflow.
The manual workflow is preferred because it creates the tag only after the
package gates pass. Once the tag exists, all registry, GitHub, and Zenodo
publication is automatic; there is no per-release Zenodo form to complete.
Never move or reuse a release tag or a published version.

The publish stages are retry-safe: the GitHub release preserves the gated
Python distributions before PyPI receives them, so a partial four-file PyPI
upload can resume with the original bytes. An existing registry file must
either match the canonical release artifact or have matching crate contents
after Cargo's generated VCS and lock metadata is removed. A version collision
with different contents stops the release and requires a version bump.

## First release migration

The existing releases predate this process and must not be moved. The first
automated release must advance every family beyond its production baseline:

| Family version | Must be greater than | Smallest patch release |
| --- | --- | --- |
| GammaLoop and monorepo | `0.3.4` | `0.3.5` |
| Spenso | `0.6.0` | `0.6.1` |
| Idenso | `0.3.0` | `0.3.1` |
| Vakint | `0.1.2` | `0.1.3` |
| Linnet | `0.17.0` | `0.17.1` |

Every changed crate whose current version already exists on crates.io must also
receive a new version in the release pull request. Release validation requires
canonical MIT text in `LICENSE.md` and requires every workspace package to
declare or inherit MIT.

Five package names are not registered on crates.io yet:

- `gammaloop-tracing-filter-macros`
- `symbolica-utils`
- `gammaloop-tracing-filter`
- `gammalooprs`
- `gammaloop-api`

After merging the first release pull request, use a maintainer token from that
exact clean `main` commit to publish every crate version in the following order
that is not already on crates.io:

1. `clinnet`
2. `gammaloop-tracing-filter-macros`
3. `spenso-macros`
4. `symbolica-utils`
5. `vakint`
6. `linnet`
7. `spenso`
8. `idenso`
9. `gammaloop-tracing-filter`
10. `spenso-hep-lib`
11. `spynso3`
12. `gammalooprs`
13. `gammaloop-api`

Run `cargo hakari publish -p <name>` and wait for each new version to appear in
the index before publishing its dependants. Publishing the whole changed chain
manually avoids asking crates.io Trusted Publishing to create one of the five
new package names. Then register Trusted Publishers for those five names and
run the `Release` workflow; it verifies and skips the already-published crate
archives before publishing PyPI and creating the tags and GitHub release. This
bootstrap is not part of later releases.

## Zenodo integration boundary

All Zenodo network operations are isolated in
[`.github/scripts/zenodo.py`](../.github/scripts/zenodo.py). It currently uses
Zenodo's legacy deposition API, which Zenodo says
[will be deprecated](https://help.zenodo.org/docs/about/whats-changed/). Keep the
record manifest and workflow payload contract independent of that API so the
client can move to Zenodo's replacement API without changing release tags or
creating replacement DOI lineages.
