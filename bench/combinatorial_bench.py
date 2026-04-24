#!/usr/bin/env python3
import argparse
import csv
import json
import math
import re
import resource
import shlex
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import combinatorial_generate

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def run_timed(cmd, stdout_path=None):
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    if stdout_path is None:
        proc = subprocess.run(cmd, text=True, capture_output=True)
    else:
        with Path(stdout_path).open('w', encoding='utf-8') as handle:
            proc = subprocess.run(cmd, text=True, stdout=handle, stderr=subprocess.PIPE)
    wall = time.perf_counter() - start
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"STDOUT:\n{getattr(proc, 'stdout', '')}\nSTDERR:\n{proc.stderr}"
        )
    return {
        'command': cmd,
        'wall_sec': wall,
        'cpu_sec': (usage_after.ru_utime + usage_after.ru_stime) - (usage_before.ru_utime + usage_before.ru_stime),
        'peak_rss_kb': usage_after.ru_maxrss,
        'stderr': proc.stderr,
        'stdout_lines': 0 if stdout_path is not None else len([ln for ln in proc.stdout.splitlines() if ln.strip()]),
    }


def parse_tsv(path):
    with Path(path).open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def ensure_dataset(args):
    dataset_dir = Path(args.dataset_dir)
    dataset_file = dataset_dir / 'dataset.json'
    requested_shape = {
        'n_blocks': args.blocks,
        'variants_per_block': args.variants,
        'n_reads': args.reads,
        'seed': args.seed,
    }
    regenerate = args.force_generate or not dataset_file.exists()
    if dataset_file.exists() and not regenerate:
        existing = json.loads(dataset_file.read_text(encoding='utf-8'))
        existing_shape = {
            'n_blocks': existing.get('n_blocks'),
            'variants_per_block': existing.get('variants_per_block'),
            'n_reads': existing.get('n_reads'),
            'seed': existing.get('seed'),
        }
        if existing_shape != requested_shape:
            print(
                f"[bench] Reusing existing dataset at {dataset_dir} with shape {existing_shape}; "
                f"requested {requested_shape}. Pass --force-generate to rebuild.",
                file=sys.stderr,
                flush=True,
            )
    if regenerate:
        print(f"[bench] Generating dataset in {dataset_dir}", file=sys.stderr, flush=True)
        cmd = [
            sys.executable,
            str(SCRIPT_DIR / 'combinatorial_generate.py'),
            '--output-dir', str(dataset_dir),
            '--blocks', str(args.blocks),
            '--variants', str(args.variants),
            '--reads', str(args.reads),
            '--seed', str(args.seed),
            '--min-block-len', str(args.min_block_len),
            '--max-block-len', str(args.max_block_len),
            '--variant-sub-rate', str(args.variant_sub_rate),
            '--read-sub-rate', str(args.read_sub_rate),
            '--read-ins-rate', str(args.read_ins_rate),
            '--read-del-rate', str(args.read_del_rate),
            '--homopolymer-boost', str(args.homopolymer_boost),
            '--clip-probability', str(args.clip_probability),
            '--clip-max', str(args.clip_max),
            '--weight-ratio', str(args.weight_ratio),
        ]
        subprocess.run(cmd, check=True)
    else:
        print(f"[bench] Using existing dataset in {dataset_dir}", file=sys.stderr, flush=True)
    return dataset_dir, json.loads(dataset_file.read_text(encoding='utf-8'))


def cigar_query_aligned_length(cigar):
    if cigar == '*' or not cigar:
        return 0
    total = 0
    for length, op in CIGAR_RE.findall(cigar):
        length = int(length)
        if op in ('M', 'I', '=', 'X'):
            total += length
    return total


def parse_nm(optional_fields):
    for field in optional_fields:
        if field.startswith('NM:i:'):
            return int(field[5:])
    return None


def parse_sam_primary(sam_path):
    ref_counts = Counter()
    by_read = {}
    headers = {'sq_count': 0, 'pg_count': 0}
    with Path(sam_path).open('r', encoding='utf-8') as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith('@'):
                if line.startswith('@SQ'):
                    headers['sq_count'] += 1
                elif line.startswith('@PG'):
                    headers['pg_count'] += 1
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11:
                continue
            qname = fields[0]
            flag = int(fields[1])
            if flag & 0x100 or flag & 0x800:
                continue
            mapped = (flag & 0x4) == 0 and fields[2] != '*'
            if not mapped:
                by_read[qname] = {
                    'mapped': False,
                    'ref_name': None,
                    'strand': None,
                    'edit_distance': None,
                    'aligned_query_length': 0,
                    'error_rate': None,
                    'correctness_percent': None,
                    'mapq': int(fields[4]),
                    'cigar': fields[5],
                }
                continue
            ref_name = fields[2]
            cigar = fields[5]
            aligned_query_length = cigar_query_aligned_length(cigar)
            edit_distance = parse_nm(fields[11:])
            error_rate = (edit_distance / aligned_query_length) if (edit_distance is not None and aligned_query_length) else None
            correctness_percent = (100.0 * (1.0 - error_rate)) if error_rate is not None else None
            strand = '-' if (flag & 0x10) else '+'
            ref_counts[ref_name] += 1
            by_read[qname] = {
                'mapped': True,
                'ref_name': ref_name,
                'strand': strand,
                'edit_distance': edit_distance,
                'aligned_query_length': aligned_query_length,
                'error_rate': error_rate,
                'correctness_percent': correctness_percent,
                'mapq': int(fields[4]),
                'cigar': cigar,
            }
    return headers, ref_counts, by_read


def mean(values):
    return (sum(values) / len(values)) if values else None


def mae(expected, actual):
    return mean([abs(a - b) for a, b in zip(expected, actual)])


def rmse(expected, actual):
    return math.sqrt(mean([(a - b) ** 2 for a, b in zip(expected, actual)])) if expected else None


def pearson(xs, ys):
    if not xs or len(xs) != len(ys):
        return None
    mx = mean(xs)
    my = mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    deny = math.sqrt(sum((y - my) ** 2 for y in ys))
    if denx == 0.0 or deny == 0.0:
        return None
    return num / (denx * deny)


def write_count_comparison(path, expected_rows, actual_counts):
    sampled = []
    actual = []
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'ref_name',
            'probability',
            'expected_count',
            'sampled_count',
            'actual_primary_count',
            'delta_actual_vs_sampled',
            'delta_actual_vs_theoretical',
            'ref_length',
        ])
        for row in expected_rows:
            sampled_count = int(row['sampled_count'])
            actual_count = int(actual_counts.get(row['ref_name'], 0))
            expected_count = float(row['expected_count'])
            sampled.append(sampled_count)
            actual.append(actual_count)
            writer.writerow([
                row['ref_name'],
                row['probability'],
                row['expected_count'],
                sampled_count,
                actual_count,
                actual_count - sampled_count,
                f"{actual_count - expected_count:.6f}",
                row['ref_length'],
            ])
    return sampled, actual


def write_read_comparison(path, truth_rows, actual_by_read):
    error_expected = []
    error_actual = []
    correctness_expected = []
    correctness_actual = []
    mapping_matches = 0
    strand_matches = 0
    mapped_reads = 0
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow([
            'read_name',
            'true_ref',
            'actual_ref',
            'expected_strand',
            'actual_strand',
            'mapped',
            'correct_ref',
            'correct_ref_and_strand',
            'expected_edit_distance',
            'actual_edit_distance',
            'expected_aligned_query_length',
            'actual_aligned_query_length',
            'expected_error_rate',
            'actual_error_rate',
            'expected_correctness_percent',
            'actual_correctness_percent',
            'mapq',
            'cigar',
        ])
        for row in truth_rows:
            actual = actual_by_read.get(row['read_name'], {'mapped': False})
            mapped = bool(actual.get('mapped'))
            actual_ref = actual.get('ref_name')
            actual_strand = actual.get('strand')
            correct_ref = mapped and actual_ref == row['true_ref']
            correct_ref_and_strand = correct_ref and actual_strand == row['strand']
            if correct_ref:
                mapping_matches += 1
            if correct_ref_and_strand:
                strand_matches += 1
            if mapped:
                mapped_reads += 1
            expected_error_rate = float(row['expected_error_rate'])
            expected_correctness_percent = float(row['expected_correctness_percent'])
            actual_error_rate = actual.get('error_rate')
            actual_correctness_percent = actual.get('correctness_percent')
            if actual_error_rate is not None:
                error_expected.append(expected_error_rate)
                error_actual.append(actual_error_rate)
                correctness_expected.append(expected_correctness_percent)
                correctness_actual.append(actual_correctness_percent)
            writer.writerow([
                row['read_name'],
                row['true_ref'],
                actual_ref or '',
                row['strand'],
                actual_strand or '',
                int(mapped),
                int(correct_ref),
                int(correct_ref_and_strand),
                row['expected_edit_distance'],
                '' if actual.get('edit_distance') is None else actual['edit_distance'],
                row['expected_aligned_query_length'],
                actual.get('aligned_query_length', 0),
                f"{expected_error_rate:.8f}",
                '' if actual_error_rate is None else f"{actual_error_rate:.8f}",
                f"{expected_correctness_percent:.6f}",
                '' if actual_correctness_percent is None else f"{actual_correctness_percent:.6f}",
                '' if actual.get('mapq') is None else actual['mapq'],
                actual.get('cigar', ''),
            ])
    return {
        'mapped_reads': mapped_reads,
        'correct_ref_reads': mapping_matches,
        'correct_ref_and_strand_reads': strand_matches,
        'error_expected': error_expected,
        'error_actual': error_actual,
        'correctness_expected': correctness_expected,
        'correctness_actual': correctness_actual,
    }


def main():
    ap = argparse.ArgumentParser(description='Generate and benchmark combinatorial long-read mapping against many references.')
    ap.add_argument('--minimap2', required=True)
    ap.add_argument('--dataset-dir', default=str(SCRIPT_DIR / 'generated' / 'combinatorial-bench'))
    ap.add_argument('--output-dir', default=None)
    ap.add_argument('--threads', type=int, default=4)
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
    ap.add_argument('--force-generate', action='store_true')
    ap.add_argument('--compact-repeats', action='store_true')
    ap.add_argument('--compact-k', type=int, default=27)
    ap.add_argument('--compact-ratio', type=float, default=0.20)
    args = ap.parse_args()

    minimap2 = str(Path(args.minimap2).resolve())
    dataset_dir, dataset = ensure_dataset(args)
    output_dir = Path(args.output_dir) if args.output_dir else (dataset_dir / 'benchmark')
    output_dir.mkdir(parents=True, exist_ok=True)

    refs = dataset_dir / 'references.fa'
    reads = dataset_dir / 'reads.fa'
    truth_path = dataset_dir / 'read_truth.tsv'
    expected_counts_path = dataset_dir / 'expected_reference_counts.tsv'
    sam_path = output_dir / 'alignments.sam'
    index_path = output_dir / 'references.mmi'
    counts_compare_path = output_dir / 'reference_count_comparison.tsv'
    read_compare_path = output_dir / 'read_error_comparison.tsv'
    summary_path = output_dir / 'benchmark.json'

    compact_args = []
    if args.compact_repeats:
        compact_args = ['--compact-repeats', '--compact-k', str(args.compact_k), '--compact-ratio', str(args.compact_ratio)]

    print(f"[bench] Building index for {dataset['n_references']} references", file=sys.stderr, flush=True)
    build_metrics = run_timed([minimap2, '-x', 'map-ont'] + compact_args + ['-d', str(index_path), str(refs)])
    print(
        f"[bench] Mapping {dataset['n_reads']} reads with {args.threads} thread(s)",
        file=sys.stderr,
        flush=True,
    )
    map_metrics = run_timed([
        minimap2,
        '-a',
        '-x', 'map-ont',
        '--secondary=no',
        '-t', str(args.threads),
        str(index_path),
        str(reads),
    ], stdout_path=sam_path)

    print("[bench] Parsing alignments and writing comparisons", file=sys.stderr, flush=True)
    expected_rows = parse_tsv(expected_counts_path)
    truth_rows = parse_tsv(truth_path)
    sam_headers, actual_counts, actual_by_read = parse_sam_primary(sam_path)
    sampled_counts, actual_primary_counts = write_count_comparison(counts_compare_path, expected_rows, actual_counts)
    read_metrics = write_read_comparison(read_compare_path, truth_rows, actual_by_read)

    mapped_reads = read_metrics['mapped_reads']
    total_reads = len(truth_rows)
    summary = {
        'dataset_dir': str(dataset_dir),
        'output_dir': str(output_dir),
        'minimap2': minimap2,
        'threads': args.threads,
        'compact_repeats': args.compact_repeats,
        'compact_k': args.compact_k if args.compact_repeats else None,
        'compact_ratio': args.compact_ratio if args.compact_repeats else None,
        'dataset': {
            'seed': dataset['seed'],
            'n_blocks': dataset['n_blocks'],
            'variants_per_block': dataset['variants_per_block'],
            'n_references': dataset['n_references'],
            'n_reads': dataset['n_reads'],
            'mean_expected_error_rate': dataset['mean_expected_error_rate'],
            'mean_expected_correctness_percent': dataset.get('mean_expected_correctness_percent'),
            'mean_read_length': dataset['mean_read_length'],
        },
        'build': build_metrics,
        'map': {
            **map_metrics,
            'reads_per_sec': (total_reads / map_metrics['wall_sec']) if map_metrics['wall_sec'] else None,
            'mapped_reads_per_sec': (mapped_reads / map_metrics['wall_sec']) if map_metrics['wall_sec'] else None,
            'sam_path': str(sam_path),
            'sam_headers': sam_headers,
        },
        'mapping_accuracy': {
            'total_reads': total_reads,
            'mapped_primary_reads': mapped_reads,
            'unmapped_reads': total_reads - mapped_reads,
            'correct_ref_reads': read_metrics['correct_ref_reads'],
            'correct_ref_and_strand_reads': read_metrics['correct_ref_and_strand_reads'],
            'mapped_rate': (mapped_reads / total_reads) if total_reads else None,
            'correct_ref_rate': (read_metrics['correct_ref_reads'] / total_reads) if total_reads else None,
            'correct_ref_and_strand_rate': (read_metrics['correct_ref_and_strand_reads'] / total_reads) if total_reads else None,
        },
        'reference_count_comparison': {
            'total_expected_sampled_reads': sum(sampled_counts),
            'total_actual_primary_reads': sum(actual_primary_counts),
            'refs_with_expected_reads': sum(1 for x in sampled_counts if x > 0),
            'refs_with_actual_primary_reads': sum(1 for x in actual_primary_counts if x > 0),
            'mae_vs_sampled_count': mae(sampled_counts, actual_primary_counts),
            'rmse_vs_sampled_count': rmse(sampled_counts, actual_primary_counts),
            'pearson_vs_sampled_count': pearson(sampled_counts, actual_primary_counts),
            'comparison_tsv': str(counts_compare_path),
        },
        'read_correctness_comparison': {
            'reads_with_actual_error_rate': len(read_metrics['error_actual']),
            'mean_expected_error_rate': mean(read_metrics['error_expected']),
            'mean_actual_error_rate': mean(read_metrics['error_actual']),
            'mae_error_rate': mae(read_metrics['error_expected'], read_metrics['error_actual']),
            'rmse_error_rate': rmse(read_metrics['error_expected'], read_metrics['error_actual']),
            'mean_expected_correctness_percent': mean(read_metrics['correctness_expected']),
            'mean_actual_correctness_percent': mean(read_metrics['correctness_actual']),
            'mae_correctness_percent': mae(read_metrics['correctness_expected'], read_metrics['correctness_actual']),
            'rmse_correctness_percent': rmse(read_metrics['correctness_expected'], read_metrics['correctness_actual']),
            'comparison_tsv': str(read_compare_path),
        },
    }

    summary_path.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(f"[bench] Wrote summary to {summary_path}", file=sys.stderr, flush=True)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
