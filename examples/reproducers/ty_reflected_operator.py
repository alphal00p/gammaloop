"""Run with Python 3.11+ and `ty check` (reproduced with ty 0.0.84).

Runtime returns Child in both orders; ty infers Base for Base() * Child().
Related: https://github.com/astral-sh/ty/issues/2432
Dispatch policy: https://github.com/astral-sh/ruff/pull/26474
"""

from __future__ import annotations

from typing import assert_type


class Base:
    def __mul__(self, other: Base) -> Base:
        return Base()


class Child(Base):
    def __mul__(self, other: Base) -> Child:
        return Child()

    def __rmul__(self, other: Base) -> Child:
        return Child()


left = Base() * Child()
right = Child() * Base()
assert_type(left, Child)  # ty: inferred type is Base
assert_type(right, Child)  # Passes
assert type(left) is Child
assert type(right) is Child
print(type(left).__name__, type(right).__name__)  # Child Child
