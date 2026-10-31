"""Frozen strict pure grammar and whole-definition algorithms. No runtime wrapper."""

import re
import collections
from functools import cache

F = "spenso::{antisymmetric,real,spenso::tensor}::f"
SPELLINGS = {
    **{
        f"gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:{n.lower()}}}::{n}": f"gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}}::{n}"
        for n in ("Q", "K", "P")
    },
    "gammalooprs::{spenso::index,spenso::index-label:e}::edge": "gammalooprs::{spenso::index}::edge",
    "gammalooprs::{spenso::index,spenso::index-label:v}::vertex": "gammalooprs::{spenso::index}::vertex",
}
TOKEN = re.compile(r"(?:[^\W\d]\w*::)+(?:\{[^{}]*\}::)?[^\W\d]\w*|𝑖|[0-9]+|[()+*/^,\-]")

class Reject(Exception):
    def __init__(self, stage, detail):
        super().__init__(detail)
        self.stage = stage

def require(condition, stage, detail):
    if not condition:
        raise Reject(stage, detail)

@cache
def parse(text, historical=False):
    require(
        isinstance(text, str) and bool(text.strip()), "syntax", "nonempty expression"
    )
    pattern = re.compile(r"[0-9]+𝑖|" + TOKEN.pattern)
    tokens, position = [], 0
    while position < len(text):
        if text[position].isspace():
            position += 1
            continue
        match = pattern.match(text, position)
        require(
            match is not None,
            "syntax",
            f"unsupported syntax at {position}: {text[position : position + 30]}",
        )
        value = match.group()
        tokens.append(SPELLINGS.get(value, value) if historical else value)
        position = match.end()
    cursor = 0

    def expression(minimum=0):
        nonlocal cursor
        require(cursor < len(tokens), "syntax", "truncated expression")
        token = tokens[cursor]
        cursor += 1
        if token == "-":
            left = ("neg", expression(30))
        elif token == "(":
            left = expression()
            require(
                cursor < len(tokens) and tokens[cursor] == ")",
                "syntax",
                "missing closing parenthesis",
            )
            cursor += 1
        elif token.isdecimal():
            left = ("integer", token)
        elif token.endswith("𝑖") and token[:-1].isdecimal():
            left = ("*", ("integer", token[:-1]), ("symbol", "𝑖"))
        else:
            require(
                "::" in token or token == "𝑖", "syntax", "invalid primary " + token
            )
            if cursor < len(tokens) and tokens[cursor] == "(":
                cursor += 1
                args = []
                if cursor < len(tokens) and tokens[cursor] != ")":
                    while True:
                        args.append(expression())
                        if cursor >= len(tokens) or tokens[cursor] != ",":
                            break
                        cursor += 1
                require(
                    cursor < len(tokens) and tokens[cursor] == ")",
                    "syntax",
                    "invalid call termination",
                )
                cursor += 1
                left = ("call", token, tuple(args))
            else:
                left = ("symbol", token)
        while cursor < len(tokens):
            op = tokens[cursor]
            precedence = {"+": 10, "-": 10, "*": 20, "/": 20, "^": 40}.get(op, -1)
            if precedence < minimum:
                break
            cursor += 1
            right = expression(precedence if op == "^" else precedence + 1)
            left = (op, left, right)
        return left

    result = expression()
    require(cursor == len(tokens), "syntax", "trailing tokens")
    return result

@cache
def render(tree):
    kind = tree[0]
    if kind in ("integer", "symbol"):
        return tree[1]
    if kind == "call":
        return tree[1] + "(" + ",".join(render(x) for x in tree[2]) + ")"
    if kind == "neg":
        return "-(" + render(tree[1]) + ")"
    return "(" + render(tree[1]) + kind + render(tree[2]) + ")"

def definitions(tree, result=None):
    if result is None:
        result = {}
    if tree[0] in ("symbol", "call"):
        spelling = tree[1]
        if spelling != "𝑖":
            name = re.sub(r"::\{[^{}]*\}::", "::", spelling)
            require(
                name not in result or result[name] == spelling,
                "symbol-definition",
                f"inconsistent definition {name}",
            )
            result[name] = spelling
    children = (
        tree[2]
        if tree[0] == "call"
        else tree[1:]
        if tree[0] not in ("symbol", "integer")
        else ()
    )
    for child in children:
        definitions(child, result)
    return result

@cache
def cyclic(tree):
    if tree[0] == "call":
        args = tuple(cyclic(x) for x in tree[2])
        if tree[1] == F:
            require(len(args) == 3, "syntax", "f arity")
            args = min(args[i:] + args[:i] for i in range(3))
        return ("call", tree[1], args)
    if tree[0] in ("integer", "symbol"):
        return tree
    return (tree[0], *(cyclic(x) for x in tree[1:]))

def owners(payload):
    result = {}
    for i, vertex in enumerate(payload["vertices"]):
        result[f"vertex:{i}"] = (vertex["numerator"], vertex.get("num_interface"))
    for i, edge in enumerate(payload["edges"]):
        require(isinstance(edge, list) and len(edge) == 2, "edge-tuple", str(i))
        require(
            isinstance(edge[0], dict) and set(edge[0]) == {"source", "target"},
            "edge-tuple",
            str(i),
        )
        result[f"edge:{i}"] = (edge[1]["numerator"], edge[1].get("num_interface"))
    for name, interface in (
        ("numerator_prefactor", "numerator_prefactor_interface"),
        ("projector", "projector_interface"),
        ("numerator", "num_interface"),
    ):
        result[name] = (payload[name], payload.get(interface))
    return result

def permutation(left, right):
    require(
        len(left) == len(set(left)) and len(right) == len(set(right)),
        "port",
        "repeated port",
    )
    require(
        collections.Counter(left) == collections.Counter(right),
        "port",
        "typed port occurrence/dimension/variance mismatch",
    )
    return [right.index(x) for x in left]
