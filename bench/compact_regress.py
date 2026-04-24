#!/usr/bin/env python3
import argparse
import json
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import compact_generate


def run(cmd, cwd=None):
    proc = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(shlex.quote(c) for c in cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc


def parse_fasta_lengths(path):
    lengths = {}
    name = None
    seq = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                lengths[name] = len(''.join(seq))
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        lengths[name] = len(''.join(seq))
    return lengths


def parse_paf(text):
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split('	')
        tags = tuple(sorted(fields[12:]))
        rows.append((fields[0], fields[4], fields[5], int(fields[1]), int(fields[2]), int(fields[3]), int(fields[6]), int(fields[7]), int(fields[8]), int(fields[11]), tags))
    return rows


def validate_paf(rows, ref_lengths):
    for row in rows:
        tname = row[2]
        tlen = row[6]
        ts = row[7]
        te = row[8]
        if tname not in ref_lengths:
            raise AssertionError(f'unknown target {tname}')
        if ref_lengths[tname] != tlen:
            raise AssertionError(f'target length mismatch for {tname}: paf={tlen} ref={ref_lengths[tname]}')
        if not (0 <= ts <= te <= tlen):
            raise AssertionError(f'invalid target interval for {tname}: {ts}-{te}/{tlen}')


def maybe_run_mappy_smoke(refs, reads):
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))
    try:
        import mappy as mp
    except Exception as exc:
        print(f'[skip] python mappy smoke unavailable: {exc}')
        return
    first_read = None
    name = None
    seq = []
    for line in Path(reads).read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                first_read = ''.join(seq)
                break
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if first_read is None and seq:
        first_read = ''.join(seq)
    if not first_read:
        raise AssertionError('python smoke could not read a query sequence')
    try:
        aln = mp.Aligner(str(refs), preset='map-ont', compact_repeats=True, compact_k=27, compact_ratio=0.20)
    except TypeError as exc:
        print(f'[skip] python mappy compact smoke unavailable: {exc}')
        return
    if not aln:
        raise AssertionError('python smoke could not build compact aligner')
    list(aln.map(first_read))
    print('[ok] python mappy compact smoke succeeded')


def ensure_dataset(dataset_dir, refs, reads):
    dataset_dir = Path(dataset_dir)
    dataset_json = dataset_dir / 'dataset.json'
    if not dataset_json.exists():
        compact_generate.main.__globals__['__name__'] = '__main__'
        subprocess.run([sys.executable, str(SCRIPT_DIR / 'compact_generate.py'), '--output-dir', str(dataset_dir), '--refs', str(refs), '--reads', str(reads)], check=True)
    return dataset_dir


def compare_mapping_modes(minimap2, refs, reads, threads, compact_args, compare_label, tempdir):
    normal_fasta = parse_paf(run([minimap2, '-c', '-x', 'map-ont', '-t', str(threads)] + compact_args + [str(refs), str(reads)]).stdout)
    out_index = Path(tempdir) / f'{compare_label}.mmi'
    run([minimap2, '-x', 'map-ont'] + compact_args + ['-d', str(out_index), str(refs)])
    indexed = parse_paf(run([minimap2, '-c', '-x', 'map-ont', '-t', str(threads), str(out_index), str(reads)]).stdout)
    if normal_fasta != indexed:
        raise AssertionError(f'{compare_label} FASTA and prebuilt index results differ at t={threads}')
    return normal_fasta


def main():
    ap = argparse.ArgumentParser(description='Regression runner for compactified-reference minimap2 mode.')
    ap.add_argument('--minimap2', required=True)
    ap.add_argument('--dataset-dir', default=str(SCRIPT_DIR / 'generated' / 'compact-regression'))
    ap.add_argument('--refs', type=int, default=128)
    ap.add_argument('--reads', type=int, default=48)
    ap.add_argument('--threads', nargs='+', type=int, default=[1, 2])
    ap.add_argument('--smoke-only', action='store_true')
    args = ap.parse_args()

    minimap2 = str(Path(args.minimap2).resolve())
    smoke_cases = [
        ('t2', Path(__file__).resolve().parents[1] / 'test' / 't2.fa', Path(__file__).resolve().parents[1] / 'test' / 'q2.fa'),
        ('x3s', Path(__file__).resolve().parents[1] / 'test' / 'x3s-ref.fa', Path(__file__).resolve().parents[1] / 'test' / 'x3s-qry.fa'),
    ]
    cases = list(smoke_cases)
    if not args.smoke_only:
        dataset_dir = ensure_dataset(args.dataset_dir, args.refs, args.reads)
        cases.append(('synthetic', dataset_dir / 'references.fa', dataset_dir / 'reads.fa'))

    with tempfile.TemporaryDirectory(prefix='compact-regress-') as tempdir:
        for label, refs, reads in cases:
            ref_lengths = parse_fasta_lengths(refs)
            for threads in args.threads:
                baseline_rows = compare_mapping_modes(minimap2, refs, reads, threads, [], f'{label}.baseline.t{threads}', tempdir)
                validate_paf(baseline_rows, ref_lengths)
                compact_rows = compare_mapping_modes(minimap2, refs, reads, threads, ['--compact-repeats', '--compact-k', '27', '--compact-ratio', '0.20'], f'{label}.compact.t{threads}', tempdir)
                validate_paf(compact_rows, ref_lengths)
                print(f'[ok] {label} t={threads}: baseline/index and compact/index agree internally; compact emitted {len(compact_rows)} records')
        smoke_refs = smoke_cases[-1][1]
        smoke_reads = smoke_cases[-1][2]
        maybe_run_mappy_smoke(smoke_refs, smoke_reads)


if __name__ == '__main__':
    main()
