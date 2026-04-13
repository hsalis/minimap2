#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def run(cmd, cwd, env=None):
	return subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True, check=True)


def parse_paf(text):
	rows = []
	for line in text.splitlines():
		if not line or line.startswith("#"):
			continue
		fields = line.split("\t")
		tags = {}
		for tag in fields[12:]:
			parts = tag.split(":", 2)
			if len(parts) == 3:
				tags[parts[0]] = parts[2]
		rows.append((
			fields[0], fields[4], fields[5],
			int(fields[2]), int(fields[3]), int(fields[7]), int(fields[8]),
			tags.get("cg", ""), int(fields[11]), tags.get("tp", ""), tags.get("ts", ".")
		))
	return sorted(rows)


def parse_python(text):
	rows = []
	for line in text.splitlines():
		if not line:
			continue
		fields = line.split("\t")
		tags = {}
		for tag in fields[10:]:
			parts = tag.split(":", 2)
			if len(parts) == 3:
				tags[parts[0]] = parts[2]
		rows.append((
			fields[0], fields[3], fields[4],
			int(fields[1]), int(fields[2]), int(fields[6]), int(fields[7]),
			tags.get("cg", ""), int(fields[9]), tags.get("tp", ""), tags.get("ts", ".")
		))
	return sorted(rows)


def ensure_dataset(repo, dataset_dir, refs, reads):
	ref = dataset_dir / "references.fa"
	fq = dataset_dir / "reads.fq"
	if ref.exists() and fq.exists():
		return ref, fq
	dataset_dir.mkdir(parents=True, exist_ok=True)
	cmd = [
		sys.executable, str(repo / "bench" / "many_targets_generate.py"),
		"--output-dir", str(dataset_dir),
		"--references", str(refs),
		"--reads", str(reads),
	]
	run(cmd, repo)
	return ref, fq


def compare(name, lhs, rhs):
	if lhs != rhs:
		raise AssertionError(f"{name} regression mismatch: {len(lhs)} baseline rows vs {len(rhs)} optimized rows")


def cli_case(repo, minimap2, out_dir, ref, reads, threads):
	sidecar = out_dir / f"{Path(ref).name}.t{threads}.mts"
	index = out_dir / f"{Path(ref).name}.t{threads}.mmi"
	base_cmd = [str(minimap2), "-x", "map-ont", "-c", "-t", str(threads), str(ref), str(reads)]
	base = parse_paf(run(base_cmd, repo).stdout)
	run([str(minimap2), "-x", "map-ont", "-d", str(index), "--write-many-targets-sidecar=yes",
	     "--many-targets-sidecar", str(sidecar), str(ref)], repo)
	opt_cmd = [str(minimap2), "-x", "map-ont", "-c", "-t", str(threads), "--many-targets",
	           "--many-targets-sidecar", str(sidecar), str(index), str(reads)]
	opt = parse_paf(run(opt_cmd, repo).stdout)
	compare(f"cli:{Path(ref).name}:t{threads}", base, opt)


def python_case(repo, out_dir, ref, reads):
	env = os.environ.copy()
	env["PYTHONPATH"] = str(repo)
	sidecar = Path(str(ref) + ".mts")
	minimap2 = repo / "minimap2"
	run([str(minimap2), "-x", "map-ont", "-d", str(out_dir / "python.mmi"), "--write-many-targets-sidecar=yes",
	     "--many-targets-sidecar", str(sidecar), str(ref)], repo)
	# Run the Python comparisons in two explicit steps to keep the output simple.
	base_cmd = (
		"import mappy as mp\n"
		f"a = mp.Aligner({str(ref)!r}, preset='map-ont')\n"
		f"reads = {str(reads)!r}\n"
		"for name, seq, qual in mp.fastx_read(reads):\n"
		"  for h in a.map(seq):\n"
		"    ts='.' if h.trans_strand == 0 else '+' if h.trans_strand > 0 else '-'\n"
		"    print('\\t'.join(map(str, (name, '+' if h.strand > 0 else '-', h.ctg, h.q_st, h.q_en, h.r_st, h.r_en, h.cigar_str, h.mapq, 'P' if h.is_primary else 'S', ts))))\n"
	)
	opt_cmd = (
		"import mappy as mp\n"
		f"a = mp.Aligner({str(ref)!r}, preset='map-ont', many_targets=True, many_targets_sidecar={str(sidecar)!r})\n"
		f"reads = {str(reads)!r}\n"
		"for name, seq, qual in mp.fastx_read(reads):\n"
		"  for h in a.map(seq):\n"
		"    ts='.' if h.trans_strand == 0 else '+' if h.trans_strand > 0 else '-'\n"
		"    print('\\t'.join(map(str, (name, '+' if h.strand > 0 else '-', h.ctg, h.q_st, h.q_en, h.r_st, h.r_en, h.cigar_str, h.mapq, 'P' if h.is_primary else 'S', ts))))\n"
	)
	base = parse_python(run([sys.executable, "-c", base_cmd], repo, env=env).stdout)
	opt = parse_python(run([sys.executable, "-c", opt_cmd], repo, env=env).stdout)
	compare(f"python:{Path(ref).name}", base, opt)


def main():
	parser = argparse.ArgumentParser(description="Regression runner for minimap2 many-target mode.")
	parser.add_argument("--minimap2", default="./minimap2")
	parser.add_argument("--dataset-dir", default="bench/generated/regression")
	parser.add_argument("--refs", type=int, default=2048)
	parser.add_argument("--reads", type=int, default=256)
	parser.add_argument("--threads", nargs="+", type=int, default=[1, 2])
	parser.add_argument("--smoke-only", action="store_true")
	parser.add_argument("--python", action="store_true")
	args = parser.parse_args()

	repo = Path(__file__).resolve().parents[1]
	minimap2 = (repo / args.minimap2).resolve() if not Path(args.minimap2).is_absolute() else Path(args.minimap2)
	out_dir = repo / args.dataset_dir
	out_dir.mkdir(parents=True, exist_ok=True)

	smoke_ref = repo / "test" / "MT-human.fa"
	smoke_reads = repo / "test" / "MT-orang.fa"
	for threads in args.threads:
		cli_case(repo, minimap2, out_dir, smoke_ref, smoke_reads, threads)

	if not args.smoke_only:
		ref, reads = ensure_dataset(repo, out_dir, args.refs, args.reads)
		for threads in args.threads:
			cli_case(repo, minimap2, out_dir, ref, reads, threads)
		if args.python:
			python_case(repo, out_dir, ref, reads)

	print(json.dumps({
		"status": "ok",
		"smoke_only": args.smoke_only,
		"threads": args.threads,
		"python": args.python,
	}))


if __name__ == "__main__":
	main()
