#!/usr/bin/env python3
import argparse
import json
import resource
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import compact_generate


def run_timed(cmd):
    rss_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    cpu_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime + resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime
    start = time.perf_counter()
    proc = subprocess.run(cmd, text=True, capture_output=True)
    wall = time.perf_counter() - start
    cpu_after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime + resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime
    rss_after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(shlex.quote(c) for c in cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return {
        'command': cmd,
        'wall_sec': wall,
        'cpu_sec': cpu_after - cpu_before,
        'peak_rss_kb': rss_after,
        'stdout_lines': len([ln for ln in proc.stdout.splitlines() if ln.strip()]),
        'stderr': proc.stderr,
    }, proc


def parse_lengths(fasta):
    total = 0
    n = 0
    seq = []
    for line in Path(fasta).read_text().splitlines():
        if line.startswith('>'):
            if seq:
                total += len(''.join(seq))
                n += 1
                seq = []
        else:
            seq.append(line.strip())
    if seq:
        total += len(''.join(seq))
        n += 1
    return n, total


def ensure_dataset(dataset_dir, refs, reads):
    dataset_dir = Path(dataset_dir)
    if not (dataset_dir / 'dataset.json').exists():
        subprocess.run([sys.executable, str(SCRIPT_DIR / 'compact_generate.py'), '--output-dir', str(dataset_dir), '--refs', str(refs), '--reads', str(reads)], check=True)
    return dataset_dir


def main():
    ap = argparse.ArgumentParser(description='Benchmark compactified-reference minimap2 mode.')
    ap.add_argument('--minimap2', required=True)
    ap.add_argument('--dataset-dir', default=str(SCRIPT_DIR / 'generated' / 'compact-bench'))
    ap.add_argument('--refs', type=int, default=512)
    ap.add_argument('--reads', type=int, default=256)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--repeats', type=int, default=3)
    ap.add_argument('--output', default=str(SCRIPT_DIR / 'results' / 'compact-bench.json'))
    args = ap.parse_args()

    minimap2 = str(Path(args.minimap2).resolve())
    dataset_dir = ensure_dataset(args.dataset_dir, args.refs, args.reads)
    refs = dataset_dir / 'references.fa'
    reads = dataset_dir / 'reads.fa'
    n_reads, total_bases = parse_lengths(reads)

    results = {
        'dataset_dir': str(dataset_dir),
        'threads': args.threads,
        'repeats': args.repeats,
        'read_count': n_reads,
        'read_bases': total_bases,
        'runs': [],
    }

    with tempfile.TemporaryDirectory(prefix='compact-bench-') as tempdir:
        tempdir = Path(tempdir)
        for mode_name, compact_args in [
            ('baseline', []),
            ('compact', ['--compact-repeats', '--compact-k', '27', '--compact-ratio', '0.20']),
        ]:
            index_path = tempdir / f'{mode_name}.mmi'
            build_metrics, _ = run_timed([minimap2, '-x', 'map-ont'] + compact_args + ['-d', str(index_path), str(refs)])
            mode_runs = []
            for _ in range(args.repeats):
                map_metrics, proc = run_timed([minimap2, '-c', '-x', 'map-ont', '-t', str(args.threads), str(index_path), str(reads)])
                map_metrics['reads_per_sec'] = n_reads / map_metrics['wall_sec'] if map_metrics['wall_sec'] else None
                map_metrics['bases_per_sec'] = total_bases / map_metrics['wall_sec'] if map_metrics['wall_sec'] else None
                map_metrics['record_count'] = len([ln for ln in proc.stdout.splitlines() if ln.strip()])
                mode_runs.append(map_metrics)
            results['runs'].append({
                'mode': mode_name,
                'build': build_metrics,
                'maps': mode_runs,
            })

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
