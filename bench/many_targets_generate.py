#!/usr/bin/env python3

import argparse
import json
import random
from pathlib import Path

DNA = "ACGT"


def revcomp(seq):
	comp = str.maketrans("ACGT", "TGCA")
	return seq.translate(comp)[::-1]


def rand_seq(rng, length):
	return "".join(rng.choice(DNA) for _ in range(length))


def decorate_reference(rng, seq, idx):
	seq = list(seq)
	if idx % 11 == 0:
		pos = rng.randrange(0, max(1, len(seq) - 64))
		run = rng.choice("AT")
		seq[pos:pos + 32] = run * 32
	if idx % 13 == 0:
		pos = rng.randrange(0, max(1, len(seq) - 80))
		motif = "ACGT" * 10
		seq[pos:pos + len(motif)] = motif
	return "".join(seq)


def mutate_read(rng, seq, sub_rate, indel_rate):
	out = []
	i = 0
	while i < len(seq):
		if rng.random() < indel_rate:
			if rng.random() < 0.5 and len(seq) - i > 1:
				i += 1
				continue
			out.append(rng.choice(DNA))
		base = seq[i]
		if rng.random() < sub_rate:
			base = rng.choice([b for b in DNA if b != base])
		out.append(base)
		i += 1
	return "".join(out)


def write_fasta(path, records):
	with path.open("w") as fh:
		for name, seq in records:
			fh.write(f">{name}\n")
			for i in range(0, len(seq), 80):
				fh.write(seq[i:i + 80] + "\n")


def write_fastq(path, records):
	with path.open("w") as fh:
		for name, seq in records:
			fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def build_dataset(args):
	rng = random.Random(args.seed)
	refs = []
	for i in range(args.references):
		length = rng.randint(args.min_ref_len, args.max_ref_len)
		if i and i % 19 == 0:
			parent = refs[rng.randrange(max(1, len(refs) // 4))][1]
			length = min(length, len(parent))
			seq = mutate_read(rng, parent[:length], 0.01, 0.002)
		else:
			seq = rand_seq(rng, length)
		seq = decorate_reference(rng, seq, i)
		refs.append((f"ref{i:05d}", seq))

	reads = []
	for i in range(args.reads):
		mode = i % 8
		if mode == 7:
			length = rng.randint(args.min_read_len, args.max_read_len)
			seq = decorate_reference(rng, rand_seq(rng, length), i + 100000)
		else:
			ref_name, ref = refs[rng.randrange(len(refs))]
			read_len = min(len(ref), rng.randint(args.min_read_len, args.max_read_len))
			start = 0 if len(ref) == read_len else rng.randrange(0, len(ref) - read_len + 1)
			seq = ref[start:start + read_len]
			if mode in (1, 5):
				seq = revcomp(seq)
			if mode in (2, 6):
				clip = min(64, len(seq) // 8)
				seq = rand_seq(rng, clip) + seq[clip:]
			if mode in (3, 4, 5, 6):
				seq = mutate_read(rng, seq, args.sub_rate * 1.35, args.indel_rate * 1.25)
			else:
				seq = mutate_read(rng, seq, args.sub_rate, args.indel_rate)
			if len(seq) < args.min_read_len:
				seq += rand_seq(rng, args.min_read_len - len(seq))
		reads.append((f"read{i:05d}", seq))

	return refs, reads


def main():
	parser = argparse.ArgumentParser(description="Generate deterministic many-target minimap2 benchmark data.")
	parser.add_argument("--output-dir", required=True)
	parser.add_argument("--references", type=int, default=10000)
	parser.add_argument("--reads", type=int, default=2000)
	parser.add_argument("--min-ref-len", type=int, default=2000)
	parser.add_argument("--max-ref-len", type=int, default=8000)
	parser.add_argument("--min-read-len", type=int, default=2000)
	parser.add_argument("--max-read-len", type=int, default=8000)
	parser.add_argument("--sub-rate", type=float, default=0.08)
	parser.add_argument("--indel-rate", type=float, default=0.03)
	parser.add_argument("--seed", type=int, default=17)
	args = parser.parse_args()

	out_dir = Path(args.output_dir)
	out_dir.mkdir(parents=True, exist_ok=True)
	refs, reads = build_dataset(args)
	ref_path = out_dir / "references.fa"
	read_path = out_dir / "reads.fq"
	meta_path = out_dir / "dataset.json"
	write_fasta(ref_path, refs)
	write_fastq(read_path, reads)
	meta = {
		"references": len(refs),
		"reads": len(reads),
		"seed": args.seed,
		"reference_path": str(ref_path),
		"reads_path": str(read_path),
	}
	meta_path.write_text(json.dumps(meta, indent=2) + "\n")
	print(json.dumps(meta))


if __name__ == "__main__":
	main()
