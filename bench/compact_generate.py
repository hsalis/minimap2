#!/usr/bin/env python3
import argparse
import json
import random
from pathlib import Path

DNA = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def rand_dna(rng, length):
    return ''.join(rng.choice(DNA) for _ in range(length))


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def mutate_read(rng, seq, sub_rate=0.035, indel_rate=0.015):
    out = []
    i = 0
    while i < len(seq):
        r = rng.random()
        if r < indel_rate * 0.5:
            i += 1
            continue
        if r < indel_rate:
            out.append(rng.choice(DNA))
            continue
        base = seq[i]
        if rng.random() < sub_rate:
            choices = [c for c in DNA if c != base]
            out.append(rng.choice(choices))
        else:
            out.append(base)
        i += 1
    return ''.join(out) if out else seq[:1]


def generate_references(rng, n_refs, min_len, max_len):
    families = max(8, min(64, n_refs // 6 or 1))
    family_blocks = [rand_dna(rng, 320) for _ in range(families)]
    references = []
    for i in range(n_refs):
        family = i % families
        target_len = rng.randint(min_len, max_len)
        parts = []
        while sum(map(len, parts)) < target_len:
            if len(parts) % 3 == 1:
                motif = family_blocks[family]
                if (i + len(parts)) % 5 == 0:
                    motif = revcomp(motif)
                parts.append(motif)
            else:
                parts.append(rand_dna(rng, rng.randint(180, 420)))
        seq = ''.join(parts)[:target_len]
        references.append({"name": f"ref_{i:05d}", "sequence": seq})
    return references


def generate_reads(rng, references, n_reads, min_len, max_len):
    reads = []
    mapped_reads = max(1, int(n_reads * 0.85))
    for i in range(mapped_reads):
        ref = references[i % len(references)]
        ref_seq = ref['sequence']
        span = min(len(ref_seq), rng.randint(min_len, min(max_len, len(ref_seq))))
        start = 0 if len(ref_seq) == span else rng.randint(0, len(ref_seq) - span)
        seq = ref_seq[start:start + span]
        if i % 4 == 1:
            seq = revcomp(seq)
        seq = mutate_read(rng, seq)
        reads.append({"name": f"read_{i:05d}", "sequence": seq, "source": ref['name']})
    for i in range(mapped_reads, n_reads):
        seq = rand_dna(rng, rng.randint(min_len, max_len))
        reads.append({"name": f"read_{i:05d}", "sequence": seq, "source": None})
    return reads


def write_fasta(path, records):
    with path.open('w') as fh:
        for rec in records:
            fh.write(f">{rec['name']}\n{rec['sequence']}\n")


def main():
    ap = argparse.ArgumentParser(description='Generate a deterministic compact-mode benchmark dataset.')
    ap.add_argument('--output-dir', required=True)
    ap.add_argument('--refs', type=int, default=256)
    ap.add_argument('--reads', type=int, default=128)
    ap.add_argument('--seed', type=int, default=13)
    ap.add_argument('--min-ref-len', type=int, default=2000)
    ap.add_argument('--max-ref-len', type=int, default=8000)
    ap.add_argument('--min-read-len', type=int, default=2000)
    ap.add_argument('--max-read-len', type=int, default=6000)
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    refs = generate_references(rng, args.refs, args.min_ref_len, args.max_ref_len)
    reads = generate_reads(rng, refs, args.reads, args.min_read_len, args.max_read_len)
    write_fasta(outdir / 'references.fa', refs)
    write_fasta(outdir / 'reads.fa', reads)
    summary = {
        'seed': args.seed,
        'refs': len(refs),
        'reads': len(reads),
        'reference_lengths': {r['name']: len(r['sequence']) for r in refs},
        'read_lengths': {r['name']: len(r['sequence']) for r in reads},
    }
    (outdir / 'dataset.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
