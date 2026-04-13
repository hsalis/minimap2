#!/usr/bin/env python3

import argparse
import json
import os
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path


BENCH_RE = re.compile(r"\[M::bench\] reads=(\d+) anchors=(\d+) avg_anchors_per_read=([0-9.]+) dp_invocations=(\d+)")
RSS_RE = re.compile(r"Peak RSS: ([0-9.]+) GB")


def run(cmd, cwd, env=None):
	start = time.perf_counter()
	proc = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True, check=True)
	elapsed = time.perf_counter() - start
	return elapsed, proc


def ensure_dataset(repo, dataset_dir, refs, reads):
	ref = dataset_dir / "references.fa"
	fq = dataset_dir / "reads.fq"
	if ref.exists() and fq.exists():
		return ref, fq
	dataset_dir.mkdir(parents=True, exist_ok=True)
	subprocess.run([
		sys.executable, str(repo / "bench" / "many_targets_generate.py"),
		"--output-dir", str(dataset_dir),
		"--references", str(refs),
		"--reads", str(reads),
	], cwd=repo, check=True)
	return ref, fq


def parse_metrics(stderr, read_len_total):
	bench = BENCH_RE.search(stderr)
	rss = RSS_RE.search(stderr)
	return {
		"reads": int(bench.group(1)) if bench else None,
		"anchors": int(bench.group(2)) if bench else None,
		"average_anchors_per_read": float(bench.group(3)) if bench else None,
		"dp_invocations": int(bench.group(4)) if bench else None,
		"peak_rss_gb": float(rss.group(1)) if rss else None,
		"bases_per_second": read_len_total,
	}


def total_query_bases(reads_path):
	total = 0
	with open(reads_path) as fh:
		for i, line in enumerate(fh):
			if i % 4 == 1:
				total += len(line.strip())
	return total


def summarize(label, samples, read_bases):
	walls = [s["wall_seconds"] for s in samples]
	cpu = [s.get("cpu_seconds") for s in samples if s.get("cpu_seconds") is not None]
	median_wall = statistics.median(walls)
	return {
		"label": label,
		"repeats": len(samples),
		"median_wall_seconds": median_wall,
		"median_reads_per_second": statistics.median(s["reads_per_second"] for s in samples),
		"median_bases_per_second": read_bases / median_wall if median_wall > 0 else None,
		"median_peak_rss_gb": statistics.median(s["peak_rss_gb"] for s in samples if s["peak_rss_gb"] is not None),
		"median_average_anchors_per_read": statistics.median(s["average_anchors_per_read"] for s in samples if s["average_anchors_per_read"] is not None),
		"median_dp_invocations": statistics.median(s["dp_invocations"] for s in samples if s["dp_invocations"] is not None),
		"samples": samples,
	}


def main():
	parser = argparse.ArgumentParser(description="Benchmark minimap2 many-target mode.")
	parser.add_argument("--minimap2", default="./minimap2")
	parser.add_argument("--dataset-dir", default="bench/generated/bench")
	parser.add_argument("--refs", type=int, default=10000)
	parser.add_argument("--reads", type=int, default=2000)
	parser.add_argument("--threads", type=int, default=2)
	parser.add_argument("--repeats", type=int, default=3)
	parser.add_argument("--output", default="bench/results/many_targets_bench.json")
	args = parser.parse_args()

	repo = Path(__file__).resolve().parents[1]
	minimap2 = (repo / args.minimap2).resolve() if not Path(args.minimap2).is_absolute() else Path(args.minimap2)
	dataset_dir = repo / args.dataset_dir
	output_path = repo / args.output
	output_path.parent.mkdir(parents=True, exist_ok=True)

	ref, reads = ensure_dataset(repo, dataset_dir, args.refs, args.reads)
	read_bases = total_query_bases(reads)
	index = dataset_dir / "references.mmi"
	sidecar = dataset_dir / "references.mmi.mts"

	subprocess.run([str(minimap2), "-x", "map-ont", "-d", str(index), str(ref)], cwd=repo, check=True)
	subprocess.run([str(minimap2), "-x", "map-ont", "-d", str(index), "--many-targets", "--write-many-targets-sidecar=yes",
	                "--many-targets-sidecar", str(sidecar), str(ref)], cwd=repo, check=True)

	env = os.environ.copy()
	env["MINIMAP2_BENCH_STATS"] = "1"
	base_samples = []
	opt_samples = []
	for label, extra, bucket in (
		("baseline", [], base_samples),
		("many_targets", ["--many-targets", "--many-targets-sidecar", str(sidecar)], opt_samples),
	):
		for _ in range(args.repeats):
			cmd = [str(minimap2), "-x", "map-ont", "-c", "-t", str(args.threads), *extra, str(index), str(reads)]
			wall, proc = run(cmd, repo, env=env)
			metrics = parse_metrics(proc.stderr, read_bases)
			reads_seen = metrics["reads"] or args.reads
			metrics.update({
				"wall_seconds": wall,
				"reads_per_second": reads_seen / wall if wall > 0 else None,
				"sidecar_used": bool(extra),
			})
			bucket.append(metrics)

	results = {
		"dataset": {
			"references": args.refs,
			"reads": args.reads,
			"reference_path": str(ref),
			"reads_path": str(reads),
		},
		"baseline": summarize("baseline", base_samples, read_bases),
		"many_targets": summarize("many_targets", opt_samples, read_bases),
	}
	results["speedup"] = (
		results["baseline"]["median_wall_seconds"] / results["many_targets"]["median_wall_seconds"]
		if results["many_targets"]["median_wall_seconds"] else None
	)
	output_path.write_text(json.dumps(results, indent=2) + "\n")
	print(json.dumps(results, indent=2))


if __name__ == "__main__":
	main()
