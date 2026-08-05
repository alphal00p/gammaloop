# Agent Guidance

Agents and human contributors follow the same repository guidance. Before making
changes, read [CONTRIBUTING.md](CONTRIBUTING.md) and treat it as authoritative
for this repository.

If instructions conflict, prefer the more specific local guidance and preserve
the user's current work unless explicitly asked to change it.

If a test is failing, and you want to change the test itself explicitly ask whether this is intended.


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
