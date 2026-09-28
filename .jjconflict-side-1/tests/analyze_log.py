# %%
import json
from datetime import datetime
from io import StringIO

import polars as pl
from symbolica import E, PrintMode, S

# 1) Lazy scan; discover all keys across the file
path = "./artifacts/feyn_gen_generation_test/logs/gammalog-feyngen.jsonl"

msg = "Initial Numerator"
# msg ="Network for canonization"


def _decode_filtered_lines(lines: list[str]) -> pl.DataFrame:
    if not lines:
        return pl.DataFrame()
    try:
        return pl.read_ndjson(StringIO("\n".join(lines)))
    except Exception:
        return pl.DataFrame([json.loads(s) for s in lines])


def df_by_message(path: str, message: str) -> pl.DataFrame:
    lf_idx = pl.scan_csv(path, has_header=False, new_columns=["raw"], separator="\n", infer_schema_length=0).with_columns(
        [
            pl.col("raw").str.json_path_match("$.timestamp").cast(pl.Utf8).alias("timestamp"),
            pl.col("raw").str.json_path_match("$.message").cast(pl.Utf8).alias("message"),
        ]
    )
    lines = lf_idx.filter(pl.col("message") == message).select("raw").collect()["raw"].to_list()
    return _decode_filtered_lines(lines)


def df_by_target(path: str, target: str) -> pl.DataFrame:
    # minimal index: keep raw + just the field we filter on
    lf_idx = pl.scan_csv(path, has_header=False, new_columns=["raw"], separator="\n", infer_schema_length=0).with_columns(
        [
            pl.col("raw").str.json_path_match("$.target").cast(pl.Utf8).alias("target"),
        ]
    )
    lines = lf_idx.filter(pl.col("target") == target).select("raw").collect()["raw"].to_list()
    return _decode_filtered_lines(lines)


def df_by_span(path: str, span: str, *, primary_only: bool = False) -> pl.DataFrame:
    """
    Filter by span name.
    - primary_only=False (default): match when span == $.span.name OR span appears in any $.spans[].name
    - primary_only=True: only match $.span.name
    """
    lf = pl.scan_csv(path, has_header=False, new_columns=["raw"], separator="\n", infer_schema_length=0)

    # Extract minimal indexes
    lf_idx = lf.with_columns(
        [
            pl.col("raw").str.json_path_match("$.span.name").cast(pl.Utf8).alias("span_name"),
            # We'll do a lightweight contains check on the JSON string for any spans[].name
            pl.col("raw").str.json_path_match("$.spans").alias("spans_json"),
        ]
    )

    if primary_only:
        pred = pl.col("span_name") == span
    else:
        # match primary OR any nested span entry; contains is fast and avoids global schema
        pred = (pl.col("span_name") == span) | (pl.col("spans_json").str.contains(f'"name":"{span}"', literal=True))

    lines = lf_idx.filter(pred).select("raw").collect()["raw"].to_list()
    return _decode_filtered_lines(lines)


# %%

num = df_by_message(path, "Initial Numerator")


samples = df_by_message(path, "Sample evaluation")

samples = df_by_message(path, "Sample evaluation")
# print(debug_dot_entries["debug_dot"][12])

# %%

num.filter(pl.col("diagram_id") == "11").glimpse()


# %%

for n in num.filter(pl.col("diagram_id") == "12").iter_rows(named=True):
    print(n["numerator"])
    try:
        print(E(n["numerator"]).coefficients_to_float(2).expand().format(show_namespaces=False, num_exp_as_superscript=False))
    except Exception as e:
        print(f"Error parsing numerator: {e}")
    print("-" * 80)


# %%
self = 12
other = 11

samples_filtered = samples.filter((pl.col("numerator_diagram_id") == str(self)) & (pl.col("denominator_diagram_id") == str(other)))


for num_row in samples_filtered.iter_rows(named=True):
    numerator_diagram_id = num_row["numerator_diagram_id"]
    denominator_diagram_id = num_row["denominator_diagram_id"]
    # ratio_value = num_row["ratio_value"]
    numerator_value = num_row["numerator"]
    denominator_value = num_row["denominator"]

    print(f"Diagram ID: {numerator_diagram_id} vs {denominator_diagram_id}")
    print("Numerator Value:")
    try:
        print(E(numerator_value).coefficients_to_float(2).expand().format(show_namespaces=False, num_exp_as_superscript=False))
    except Exception as e:
        print(f"Error parsing numerator_value: {e}")
        print(numerator_value)
    print("Denominator Value:")
    try:
        print(E(denominator_value).coefficients_to_float(2).expand().format(show_namespaces=False, num_exp_as_superscript=False))
    except Exception as e:
        print(f"Error parsing denominator_value: {e}")
        print(denominator_value)
    # print("Ratio Value:")
    # try:
    #     print(E(ratio_value).coefficients_to_float(2))
    # except Exception as e:
    #     print(f"Error parsing ratio_value: {e}")
    #     print(ratio_value)
    print("-" * 80)

# %%


# %%
