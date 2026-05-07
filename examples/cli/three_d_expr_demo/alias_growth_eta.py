#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
import time
from pathlib import Path
from statistics import median


def read_last_row(path):
    if not path.exists():
        return None
    last = None
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            last = row
    return last


def process_elapsed_seconds(pid):
    out = subprocess.check_output(["ps", "-o", "etime=", "-p", str(pid)], text=True).strip()
    if not out:
        return None
    days = 0
    if "-" in out:
        day_text, out = out.split("-", 1)
        days = int(day_text)
    parts = [int(part) for part in out.split(":")]
    if len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:
        hours = 0
        minutes, seconds = parts
    else:
        hours = minutes = 0
        seconds = parts[0]
    return float((((days * 24) + hours) * 60 + minutes) * 60 + seconds)


def process_rss_kb(pid):
    out = subprocess.check_output(["ps", "-o", "rss=", "-p", str(pid)], text=True).strip()
    return int(out) if out else 0


def append_sample(samples_path, sample):
    samples_path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not samples_path.exists()
    with samples_path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sample))
        if new_file:
            writer.writeheader()
        writer.writerow(sample)


def read_samples(samples_path):
    if not samples_path.exists():
        return []
    with samples_path.open(newline="") as handle:
        return [
            {
                "timestamp": float(row["timestamp"]),
                "j": int(row["j"]),
                "elapsed": float(row["elapsed"]),
                "rss_kb": int(row["rss_kb"]),
            }
            for row in csv.DictReader(handle)
        ]


def fit_eta_from_samples(samples, current_j, target_j):
    intervals = []
    for left, right in zip(samples, samples[1:]):
        dj = right["j"] - left["j"]
        dt = right["timestamp"] - left["timestamp"]
        if dj > 0 and dt > 0:
            intervals.append(((left["j"] + right["j"]) / 2.0, dt / dj, dj))
    if len(intervals) < 8:
        return None

    xs = [math.log(max(mid, 1.0)) for mid, _sec, _w in intervals]
    ys = [math.log(max(sec, 1.0e-9)) for _mid, sec, _w in intervals]
    slopes = []
    for i in range(len(xs)):
        for k in range(i + 1, len(xs)):
            dx = xs[k] - xs[i]
            if abs(dx) > 1.0e-12:
                slopes.append((ys[k] - ys[i]) / dx)
    if not slopes:
        return None
    power = min(4.0, max(0.0, median(slopes)))
    log_a = median(y - power * x for x, y in zip(xs, ys))
    remaining = sum(math.exp(log_a) * (j**power) for j in range(current_j + 1, target_j + 1))
    return remaining, f"live fit dt/dj ~ j^{power:.2f}"


def fit_eta_from_atom_bytes(csv_path, current_j, target_j, elapsed):
    rows = []
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append((int(row["j"]), float(row["atom_bytes"])))
    if len(rows) < 8 or elapsed <= 0:
        return None

    xs = [math.log(j) for j, b in rows if j > 0 and b > 0]
    ys = [math.log(b) for j, b in rows if j > 0 and b > 0]
    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if denom <= 0:
        return None
    power = (n * sxy - sx * sy) / denom
    log_a = (sy - power * sx) / n

    done_weight = sum(math.exp(log_a) * (j**power) for j in range(1, current_j + 1))
    remaining_weight = sum(math.exp(log_a) * (j**power) for j in range(current_j + 1, target_j + 1))
    if done_weight <= 0:
        return None
    return elapsed * remaining_weight / done_weight, f"fallback atom_bytes ~ j^{power:.2f}"


def format_duration(seconds):
    seconds = max(0, int(seconds))
    hours, rem = divmod(seconds, 3600)
    minutes, seconds = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes}m"
    if minutes:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--pid", required=True, type=int)
    parser.add_argument("--target", type=int, default=686)
    args = parser.parse_args()

    row = read_last_row(args.csv)
    if row is None:
        print("ETA unavailable: no streamed rows yet")
        return

    j = int(row["j"])
    elapsed = process_elapsed_seconds(args.pid)
    rss_kb = process_rss_kb(args.pid)
    append_sample(
        args.samples,
        {
            "timestamp": f"{time.time():.6f}",
            "j": j,
            "elapsed": f"{elapsed:.6f}",
            "rss_kb": rss_kb,
        },
    )
    samples = read_samples(args.samples)

    eta = fit_eta_from_samples(samples, j, args.target)
    if eta is None and elapsed is not None:
        eta = fit_eta_from_atom_bytes(args.csv, j, args.target, elapsed)

    if eta is None:
        print(f"ETA unavailable: j={j}/{args.target}, rss={rss_kb / 1024 / 1024:.2f} GiB")
        return

    remaining, method = eta
    finish = time.strftime("%H:%M:%S", time.localtime(time.time() + remaining))
    print(
        f"ETA {format_duration(remaining)} (finish ~{finish}); "
        f"j={j}/{args.target}, rss={rss_kb / 1024 / 1024:.2f} GiB, {method}"
    )


if __name__ == "__main__":
    main()
