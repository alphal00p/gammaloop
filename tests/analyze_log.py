import polars as pl
from datetime import datetime
from symbolica import E,S,PrintMode

# Load the JSONL tracing data
df = pl.read_ndjson("./workspace/feyn_gen_generation_test/logs/gammalog-feyngen.jsonl")

df_with_spans = df.with_columns([
    pl.col("span").struct.field("self_id").cast(pl.Int32, strict=False).alias("self_id"),
    pl.col("span").struct.field("other_id").cast(pl.Int32, strict=False).alias("other_id"),
    pl.col("span").struct.field("name").alias("span_name"),
    pl.col("timestamp").str.to_datetime(format="%Y-%m-%dT%H:%M:%S%.fZ").alias("timestamp_parsed")
])

# Get comparison pairs
comparison_pairs = df_with_spans.filter(
    pl.col("message") == "Starting scalar rescaling comparison between diagrams"
).select(["self_id", "other_id", "timestamp_parsed"])

print(f"Total comparisons: {len(comparison_pairs)}")
print(comparison_pairs)
debug_dot_entries = df.filter(
    pl.col("span").is_not_null() & pl.col("span").struct.field("debug_dot").is_not_null()
).select([
    pl.col("span").struct.field("debug_dot").alias("debug_dot"),
    pl.col("timestamp"),
    pl.col("message")
])
print(debug_dot_entries["debug_dot"][12])



# self_diagram_id = %self.diagram_id,
# other_diagram_id = %other.diagram_id,
# ratio_value = ?rat.as_ref().map(|ra| ra.floatify(13).to_canonical_string()).unwrap_or("None".into()),
# numerator_value = %a.floatify(13).to_canonical_string(),
# denominator_value = %b.floatify(13).to_canonical_string(),
samples = df_with_spans.filter(
    pl.col("message")=="Detailed sample evaluation ratio information"
).select([pl.col("timestamp_parsed"),pl.col("self_diagram_id").map_elements(lambda x: x if x is None else int(x),return_dtype=pl.Int32).alias("self_diagram_id"),pl.col("other_diagram_id").map_elements(lambda x: x if x is None else int(x),return_dtype=pl.Int32).alias("other_diagram_id"),pl.col("ratio_value").str.strip_chars('"'),pl.col("numerator_value"),pl.col("denominator_value")])

samples_filtered = samples.filter(
    (pl.col("self_diagram_id")==6) & (pl.col("other_diagram_id")==5)
)


for num_row in samples_filtered.iter_rows(named=True):
    self_diagram_id = num_row["self_diagram_id"]
    other_diagram_id = num_row["other_diagram_id"]
    ratio_value = num_row["ratio_value"]
    numerator_value = num_row["numerator_value"]
    denominator_value = num_row["denominator_value"]

    print(f"Diagram ID: {self_diagram_id} vs {other_diagram_id}")
    print("Numerator Value:")
    try:
        print(E(numerator_value).format(show_namespaces=False,num_exp_as_superscript=False))
    except Exception as e:
        print(f"Error parsing numerator_value: {e}")
        print(numerator_value)
    print("Denominator Value:")
    try:
        print(E(denominator_value).format(show_namespaces=False,num_exp_as_superscript=False))
    except Exception as e:
        print(f"Error parsing denominator_value: {e}")
        print(denominator_value)
    print("Ratio Value:")
    try:
        print(E(ratio_value))
    except Exception as e:
        print(f"Error parsing ratio_value: {e}")
        print(ratio_value)
    print("-" * 80)
