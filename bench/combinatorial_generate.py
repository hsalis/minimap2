#!/usr/bin/env python3
import argparse
import csv
import itertools
import json
import math
import random
from collections import Counter
from pathlib import Path

DNA = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def rand_dna(rng, length):
    return ''.join(rng.choice(DNA) for _ in range(length))


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def mutate_consensus(rng, seq, sub_rate):
    out = []
    for base in seq:
        if rng.random() < sub_rate:
            choices = [c for c in DNA if c != base]
            out.append(rng.choice(choices))
        else:
            out.append(base)
    return ''.join(out)


def geometric_weights(n, ratio):
    weights = [ratio ** i for i in range(n)]
    total = sum(weights)
    return [w / total for w in weights]


def choose_variant(rng, weights):
    r = rng.random()
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r <= acc:
            return i
    return len(weights) - 1


def sample_clip_length(rng, clip_probability, clip_max):
    if rng.random() >= clip_probability:
        return 0
    x = int(round(-math.log(max(1e-9, 1.0 - rng.random())) * (clip_max / 5.0)))
    return max(0, min(clip_max, x))


def ont_mutate(rng, seq, sub_rate, ins_rate, del_rate, homopolymer_boost):
    out = []
    i = 0
    edit_distance = 0
    while i < len(seq):
        base = seq[i]
        prev_same = i > 0 and seq[i - 1] == base
        next_same = i + 1 < len(seq) and seq[i + 1] == base
        hp_factor = homopolymer_boost if (prev_same or next_same) else 1.0
        cur_del = min(0.95, del_rate * hp_factor)
        cur_ins = min(0.95, ins_rate * hp_factor)
        r = rng.random()
        if r < cur_del:
            edit_distance += 1
            i += 1
            continue
        if r < cur_del + cur_ins:
            ins_len = 1 + (1 if rng.random() < 0.12 else 0)
            for _ in range(ins_len):
                out.append(rng.choice(DNA))
                edit_distance += 1
            continue
        if rng.random() < sub_rate:
            choices = [c for c in DNA if c != base]
            out.append(rng.choice(choices))
            edit_distance += 1
        else:
            out.append(base)
        i += 1
    mutated = ''.join(out)
    if not mutated:
        mutated = seq[:1]
        edit_distance = max(1, edit_distance)
    return mutated, edit_distance


def build_blocks(rng, n_blocks, n_variants, min_block_len, max_block_len, variant_sub_rate):
    blocks = []
    for _ in range(n_blocks):
        block_len = rng.randint(min_block_len, max_block_len)
        consensus = rand_dna(rng, block_len)
        variants = []
        seen = set()
        while len(variants) < n_variants:
            seq = mutate_consensus(rng, consensus, variant_sub_rate)
            if seq in seen:
                continue
            seen.add(seq)
            variants.append(seq)
        blocks.append(variants)
    return blocks


def combo_name(combo):
    return ''.join(str(v + 1) for v in combo)


def combo_sequence(blocks, combo):
    return ''.join(blocks[i][combo[i]] for i in range(len(combo)))


def write_references(path, blocks, weights, n_reads, sampled_counts):
    reference_lengths = {}
    with path.open('w', encoding='utf-8') as handle:
        for combo in itertools.product(range(len(blocks[0])), repeat=len(blocks)):
            name = combo_name(combo)
            seq = combo_sequence(blocks, combo)
            handle.write(f'>{name}\n{seq}\n')
            reference_lengths[name] = len(seq)
    return reference_lengths


def write_block_fasta(path, blocks):
    with path.open('w', encoding='utf-8') as handle:
        for i, variants in enumerate(blocks, start=1):
            for j, seq in enumerate(variants, start=1):
                handle.write(f'>b_{i}_{j}\n{seq}\n')


def write_expected_counts(path, blocks, weights, n_reads, sampled_counts, ref_lengths):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['ref_name', 'probability', 'expected_count', 'sampled_count', 'ref_length'])
        for combo in itertools.product(range(len(blocks[0])), repeat=len(blocks)):
            name = combo_name(combo)
            prob = 1.0
            for i, variant_idx in enumerate(combo):
                prob *= weights[i][variant_idx]
            writer.writerow([
                name,
                f'{prob:.12f}',
                f'{n_reads * prob:.6f}',
                sampled_counts.get(name, 0),
                ref_lengths[name],
            ])


def generate_reads(path, truth_path, blocks, weights, n_reads, rng, sub_rate, ins_rate, del_rate, clip_probability, clip_max, homopolymer_boost):
    sampled_counts = Counter()
    truth_rows = []

    with path.open('w', encoding='utf-8') as fasta_handle:
        for idx in range(n_reads):
            combo = tuple(choose_variant(rng, weights[i]) for i in range(len(blocks)))
            ref_name = combo_name(combo)
            ref_seq = combo_sequence(blocks, combo)
            sampled_counts[ref_name] += 1

            strand = '-' if rng.random() < 0.5 else '+'
            template = revcomp(ref_seq) if strand == '-' else ref_seq
            mutated_core, edit_distance = ont_mutate(
                rng,
                template,
                sub_rate=sub_rate,
                ins_rate=ins_rate,
                del_rate=del_rate,
                homopolymer_boost=homopolymer_boost,
            )
            left_clip = sample_clip_length(rng, clip_probability, clip_max)
            right_clip = sample_clip_length(rng, clip_probability, clip_max)
            read_seq = rand_dna(rng, left_clip) + mutated_core + rand_dna(rng, right_clip)
            read_name = f'read_{idx:07d}'
            fasta_handle.write(f'>{read_name}\n{read_seq}\n')

            expected_error_rate = (edit_distance / len(mutated_core)) if mutated_core else 0.0
            truth_rows.append({
                'read_name': read_name,
                'true_ref': ref_name,
                'strand': strand,
                'expected_edit_distance': edit_distance,
                'expected_aligned_query_length': len(mutated_core),
                'expected_error_rate': expected_error_rate,
                'expected_correctness_percent': 100.0 * (1.0 - expected_error_rate),
                'left_soft_clip': left_clip,
                'right_soft_clip': right_clip,
                'read_length': len(read_seq),
            })

    with truth_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(truth_rows[0].keys()), delimiter='\t')
        writer.writeheader()
        writer.writerows(truth_rows)

    return sampled_counts, truth_rows


def main():
    ap = argparse.ArgumentParser(description='Generate a combinatorial DNA assembly benchmark with ONT-like reads.')
    ap.add_argument('--output-dir', required=True)
    ap.add_argument('--blocks', type=int, default=6)
    ap.add_argument('--variants', type=int, default=6)
    ap.add_argument('--reads', type=int, default=12000)
    ap.add_argument('--seed', type=int, default=23)
    ap.add_argument('--min-block-len', type=int, default=380)
    ap.add_argument('--max-block-len', type=int, default=520)
    ap.add_argument('--variant-sub-rate', type=float, default=0.10)
    ap.add_argument('--read-sub-rate', type=float, default=0.018)
    ap.add_argument('--read-ins-rate', type=float, default=0.012)
    ap.add_argument('--read-del-rate', type=float, default=0.020)
    ap.add_argument('--homopolymer-boost', type=float, default=1.6)
    ap.add_argument('--clip-probability', type=float, default=0.30)
    ap.add_argument('--clip-max', type=int, default=160)
    ap.add_argument('--weight-ratio', type=float, default=0.62)
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    blocks = build_blocks(
        rng,
        n_blocks=args.blocks,
        n_variants=args.variants,
        min_block_len=args.min_block_len,
        max_block_len=args.max_block_len,
        variant_sub_rate=args.variant_sub_rate,
    )
    weights = [geometric_weights(args.variants, args.weight_ratio) for _ in range(args.blocks)]

    reads_path = out_dir / 'reads.fa'
    truth_path = out_dir / 'read_truth.tsv'
    sampled_counts, truth_rows = generate_reads(
        reads_path,
        truth_path,
        blocks,
        weights,
        args.reads,
        rng,
        sub_rate=args.read_sub_rate,
        ins_rate=args.read_ins_rate,
        del_rate=args.read_del_rate,
        clip_probability=args.clip_probability,
        clip_max=args.clip_max,
        homopolymer_boost=args.homopolymer_boost,
    )

    write_block_fasta(out_dir / 'blocks.fa', blocks)
    ref_lengths = write_references(out_dir / 'references.fa', blocks, weights, args.reads, sampled_counts)
    write_expected_counts(out_dir / 'expected_reference_counts.tsv', blocks, weights, args.reads, sampled_counts, ref_lengths)

    summary = {
        'seed': args.seed,
        'n_blocks': args.blocks,
        'variants_per_block': args.variants,
        'n_references': args.variants ** args.blocks,
        'n_reads': args.reads,
        'weights': weights,
        'block_lengths': [len(blocks[i][0]) for i in range(args.blocks)],
        'reference_lengths': ref_lengths,
        'sampled_reference_counts': dict(sampled_counts),
        'mean_expected_error_rate': sum(row['expected_error_rate'] for row in truth_rows) / len(truth_rows),
        'mean_expected_correctness_percent': sum(row['expected_correctness_percent'] for row in truth_rows) / len(truth_rows),
        'mean_read_length': sum(row['read_length'] for row in truth_rows) / len(truth_rows),
    }
    (out_dir / 'dataset.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
