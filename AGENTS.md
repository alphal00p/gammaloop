# Agent Guidance

Agents and human contributors follow the same repository guidance. Before making
changes, read [CONTRIBUTING.md](CONTRIBUTING.md) and treat it as authoritative
for this repository.

If instructions conflict, prefer the more specific local guidance and preserve
the user's current work unless explicitly asked to change it.

Ask for approval before changing a test expectation only when the change is
non-trivial or deep. Fix routine fixture mistakes and straightforward obsolete
expectations without asking, preserving the test's validation purpose.

Floating-point underflow in exponentially suppressed tails may round to zero,
including when this makes stability comparisons trivially pass. Apply routine
underflow-to-zero corrections and update their test expectations without asking
again. Preserve meaningful contributions by combining numerical factors before
rounding; overflow and invalid numerical operations remain separate issues.
For stability of an averaged component, its mean absolute fully weighted probe
value may establish underflow below binary64's smallest normal value. Keep this
bound separate for each component and observable; never borrow an imaginary
component's scale to accept a meaningful real discrepancy.

Before adding helper functions, structs, or methods, check the codebase for
similar use cases, and whether the functionality is already provided by the
existing code or only needs a small adjustment/API change. When adding a new
helper, confirm with the codebase maintainers that the functionality is not
already provided by an existing helper.

Be as idiomatic and concise as possible to optimize for readability. Before
finishing a turn, check this. Anything that can be done more concisely should be
done so, and we prefer diffs with more deletions than more additions, as long as
the code is still correct and readable.

Never delete comments outright, even when they are refactored (i.e. move them
along with the code). Try to keep comments up-to-date with the code, and add new
comments as needed. Ask for confirmation if you think a comment is no longer
needed.
