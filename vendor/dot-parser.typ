= DOT parser complete-document validation

`vendor/dot-parser` contains the crates.io source of `dot-parser` 0.5.1
(GPL-2.0-or-later), with its upstream package metadata and examples retained.

The local patch adds an end-of-input assertion to the `dotfile` grammar rule.
Upstream accepts a valid prefix of a document and silently discards a malformed
subsequent graph. Bulk FeynKit imports must reject such input rather than return
an incomplete list. Negative lookahead avoids introducing an extra parse-tree
node into the upstream AST conversion.

All workspace consumers, including Linnet and the standalone notebook host,
resolve this same source through the workspace dependency. Remove the local copy
once an upstream release provides equivalent complete-document validation.
