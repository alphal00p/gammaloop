"""Usage: python3 -I ggppss_catalogue.py <masters.jsonl> <out catalogue.json>
          [--data-version ggppss-5L-v1] [--terminal-map <map.json>] [--self-check]

Converts the GGPPSS five-loop masters of vakint-external-data
(p_integral/masters/5L-ggppss/masters.jsonl, arXiv:2104.08272, already converted to
Symbolica text by tools/ggppss_convert.py there) into a vakint p master catalogue
(schema vakint.p-master-catalogue.v1, src/p_integral/master_catalogue.rs):

- id "ggppss_int(k)", symbol "ggppss_int_k" (rendered vakint::PMaster(ggppss_int_k));
- leading_order from the source; every order up to the first order holding a
  vakint::pint_undetermined symbol is an exact coefficient; from that order on the
  master is unknown ("unknown_from"), so the finalizer's Oep guard refuses any result
  that needs it. A coefficient naming another master's undetermined constant counts as
  undetermined too (xd-ggppss finding (b));
- the normalisation is vakint's p-master one (G-scheme at Q^2 = 1, every loop divided
  by ep*G(1,1); the bubble is 1/ep), which is the GGPPSS normalisation (their eq. (3)).

--terminal-map merges the "terminals" rows of a WP3 terminal -> master map (rows
{node, key, status, combination: [{master, coefficient}], note}; masters by id
"ggppss_int(k)" or symbol). Stdlib only; the output is canonical (sorted keys, one
master per line).
"""
import json
import re
import sys

SCHEMA = "vakint.p-master-catalogue.v1"
VERSION = 1  # P_MASTER_CATALOGUE_VERSION
NORMALISATION = "g-scheme"  # P_MASTER_NORMALISATION
UNDETERMINED = "pint_undetermined"
ALLOWED = re.compile(r"^[0-9+\-*/^() ,.]*$")


def strip_constants(text):
    """The text without the constants vakint's catalogue accepts."""
    for pattern in (r"zeta\(\d+\)", r"MZV\([\d,]+\)", r"\bpi\b"):
        text = re.sub(pattern, "", text)
    return text


def convert_master(entry):
    k = entry["id"]
    leading = entry["leading_order"]
    orders = {int(order): text for order, text in entry["orders"].items()}
    last = entry["last_order"]
    if sorted(orders) != list(range(leading, last + 1)):
        raise ValueError(f"ggppss_int({k}): orders are not contiguous from {leading} to {last}")
    unknown_from = last + 1
    for order in range(leading, last + 1):
        if UNDETERMINED in orders[order]:
            unknown_from = order
            break
    if entry.get("undetermined_from") is not None and entry["undetermined_from"] != unknown_from:
        raise ValueError(f"ggppss_int({k}): undetermined_from {entry['undetermined_from']} != {unknown_from}")
    coefficients = {}
    for order in range(leading, unknown_from):
        text = orders[order].replace("vakint::", "")
        if not ALLOWED.match(strip_constants(text)):
            raise ValueError(f"ggppss_int({k}) epsilon^{order}: unexpected symbol in {text[:80]}")
        coefficients[str(order)] = {"exact": text}
    return {
        "id": f"ggppss_int({k})",
        "symbol": f"ggppss_int_{k}",
        "loops": entry["loops"],
        "graph": f"GGPPSS int[{k}]: {entry['lines']} lines, {entry['dots']} dots, {entry['components']} component(s)",
        "leading_order": leading,
        "known_through": unknown_from - 1,
        "unknown_from": unknown_from,
        "coefficients": coefficients,
    }


def convert(lines, data_version, terminals=()):
    masters = [convert_master(json.loads(line)) for line in lines if line.strip()]
    return {
        "schema": SCHEMA,
        "version": VERSION,
        "data_version": data_version,
        "normalisation": NORMALISATION,
        "description": "GGPPSS five-loop master p-integrals (arXiv:2104.08272), exact through the "
        "first undetermined order; values only, the terminal rows come from the WP3 map",
        "masters": masters,
        "terminals": list(terminals),
    }


def dump(catalogue):
    head = {key: value for key, value in catalogue.items() if key not in ("masters", "terminals")}
    out = json.dumps(head, sort_keys=True, ensure_ascii=False)[:-1]
    for key in ("masters", "terminals"):
        rows = catalogue[key]
        out += f', "{key}": [' + ("\n" if rows else "")
        out += ",\n".join(json.dumps(row, sort_keys=True, ensure_ascii=False) for row in rows)
        out += "]"
    return out + "}\n"


def self_check(catalogue, lines):
    """Counts of the conversion: 281 masters, 135 with an undetermined tail (xd-ggppss),
    and the histogram of the known depth (P5L_STATUS.json basis)."""
    masters = catalogue["masters"]
    sources = {json.loads(line)["id"]: json.loads(line) for line in lines if line.strip()}
    tails = sum(1 for m in masters
                if sources[int(m["id"][11:-1])]["undetermined_from"] is not None)
    histogram = {}
    for m in masters:
        histogram[m["known_through"]] = histogram.get(m["known_through"], 0) + 1
    assert len(masters) == 281, len(masters)
    assert tails == 135, tails
    return {"masters": len(masters), "undetermined_tails": tails,
            "known_through": dict(sorted(histogram.items()))}


def main(argv):
    if len(argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    source, target = argv[1], argv[2]
    data_version = "ggppss-5L-v1"
    terminals = []
    check = False
    rest = argv[3:]
    while rest:
        flag = rest.pop(0)
        if flag == "--data-version":
            data_version = rest.pop(0)
        elif flag == "--terminal-map":
            with open(rest.pop(0), encoding="utf-8") as handle:
                terminals = json.load(handle)["terminals"]
        elif flag == "--self-check":
            check = True
        else:
            print(f"unknown argument {flag}", file=sys.stderr)
            return 2
    with open(source, encoding="utf-8") as handle:
        lines = handle.readlines()
    catalogue = convert(lines, data_version, terminals)
    if check:
        print(json.dumps(self_check(catalogue, lines)))
    with open(target, "w", encoding="utf-8") as handle:
        handle.write(dump(catalogue))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
