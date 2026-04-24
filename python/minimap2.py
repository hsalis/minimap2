#!/usr/bin/env python

import argparse
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import mappy as mp

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=XB])")
DEFAULT_OUTPUT_FILES = (
    "coverage_percent.npy",
    "high_quality_mappings.npy",
    "matching_percent.npy",
    "mean_depth.npy",
    "offsets.npy",
    "oligomers.tsv",
    "ref_lengths.npy",
    "ref_names.npy",
    "summary.tsv",
    "total_mappings.npy",
)


def _bool_flag(value):
    if value is None:
        return True
    value = value.strip().lower()
    if value in ("yes", "y", "true", "t", "1"):
        return True
    if value in ("no", "n", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("expected yes or no")


def _count_reads(query_path):
    count = 0
    for _name, _seq, _qual in mp.fastx_read(query_path):
        count += 1
    return count


def _parse_tags(tag_fields):
    tags = {}
    for field in tag_fields:
        parts = field.split(":", 2)
        if len(parts) != 3:
            continue
        key, tag_type, raw_value = parts
        if tag_type in ("i", "I"):
            try:
                value = int(raw_value)
            except ValueError:
                value = raw_value
        elif tag_type == "f":
            try:
                value = float(raw_value)
            except ValueError:
                value = raw_value
        else:
            value = raw_value
        tags[key] = value
    return tags


def _cigar_metrics(cigar, ref_start):
    if not cigar or cigar == "*":
        return {
            "intervals": [],
            "aligned_ref_bases": 0,
            "insertions": 0,
            "deletions": 0,
            "m_bases": 0,
            "eq_bases": 0,
            "x_bases": 0,
            "ref_end": ref_start,
        }

    intervals = []
    segment_start = None
    ref_pos = ref_start
    aligned_ref_bases = 0
    insertions = 0
    deletions = 0
    m_bases = 0
    eq_bases = 0
    x_bases = 0

    for length_s, op in CIGAR_RE.findall(cigar):
        length = int(length_s)
        if op in ("M", "=", "X", "D"):
            if segment_start is None:
                segment_start = ref_pos
            if op == "M":
                m_bases += length
            elif op == "=":
                eq_bases += length
            elif op == "X":
                x_bases += length
            elif op == "D":
                deletions += length
            aligned_ref_bases += length
            ref_pos += length
        elif op == "N":
            if segment_start is not None and ref_pos > segment_start:
                intervals.append((segment_start, ref_pos))
                segment_start = None
            ref_pos += length
        elif op == "I":
            insertions += length
        else:
            continue

    if segment_start is not None and ref_pos > segment_start:
        intervals.append((segment_start, ref_pos))

    return {
        "intervals": intervals,
        "aligned_ref_bases": aligned_ref_bases,
        "insertions": insertions,
        "deletions": deletions,
        "m_bases": m_bases,
        "eq_bases": eq_bases,
        "x_bases": x_bases,
        "ref_end": ref_pos,
    }


def _matched_bases_from_nm(m_bases, eq_bases, x_bases, insertions, deletions, nm, fallback=None):
    if fallback is not None:
        return max(0, int(fallback))
    if nm is None:
        return max(0, eq_bases + m_bases - x_bases)
    mismatches_in_m = max(0, int(nm) - int(insertions) - int(deletions) - int(x_bases))
    return max(0, int(eq_bases) + int(m_bases) - mismatches_in_m)


def _load_reference_metadata(ref_path, aligner, fmt, raw_headers=None):
    ref_names = []
    ref_lengths = []

    if raw_headers:
        ref_names = [name for name, _len in raw_headers]
        ref_lengths = [int(length) for _name, length in raw_headers]
    elif ref_path and str(ref_path).lower().endswith(".mmi"):
        ref_names = list(aligner.seq_names)
        ref_lengths = [len(aligner.seq(name)) for name in ref_names]
    else:
        for name, seq, _qual in mp.fastx_read(ref_path):
            ref_names.append(name)
            ref_lengths.append(len(seq))

    ref_index = {name: i for i, name in enumerate(ref_names)}
    return ref_names, ref_lengths, ref_index


def _parse_legacy_line(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    tags = _parse_tags(fields[12:])
    return {
        "qname": fields[0],
        "qlen": int(fields[1]),
        "qstart": int(fields[2]),
        "qend": int(fields[3]),
        "strand": fields[4],
        "rname": fields[5],
        "rlen": int(fields[6]),
        "rstart": int(fields[7]),
        "rend": int(fields[8]),
        "mlen": int(fields[9]),
        "blen": int(fields[10]),
        "mapq": int(fields[11]),
        "tags": tags,
        "cigar": tags.get("cg"),
        "nm": tags.get("NM"),
        "divergence": tags.get("de"),
        "is_unmapped": False,
        "is_secondary": tags.get("tp") == "S",
        "is_supplementary": False,
    }


def _parse_paf_line(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    tags = _parse_tags(fields[12:])
    return {
        "qname": fields[0],
        "qlen": int(fields[1]),
        "qstart": int(fields[2]),
        "qend": int(fields[3]),
        "strand": fields[4],
        "rname": fields[5],
        "rlen": int(fields[6]),
        "rstart": int(fields[7]),
        "rend": int(fields[8]),
        "mlen": int(fields[9]),
        "blen": int(fields[10]),
        "mapq": int(fields[11]),
        "tags": tags,
        "cigar": tags.get("cg"),
        "nm": tags.get("NM"),
        "divergence": tags.get("de"),
        "is_unmapped": False,
        "is_secondary": tags.get("tp") == "S",
        "is_supplementary": False,
    }


def _parse_sam_line(line, ref_lengths_by_name):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        return None
    flag = int(fields[1])
    tags = _parse_tags(fields[11:])
    rname = fields[2]
    cigar = fields[5]
    rstart = int(fields[3]) - 1 if fields[3] != "0" else 0
    metrics = _cigar_metrics(cigar, rstart)
    rlen = int(ref_lengths_by_name.get(rname, 0))
    nm = tags.get("NM")
    matched = _matched_bases_from_nm(
        metrics["m_bases"],
        metrics["eq_bases"],
        metrics["x_bases"],
        metrics["insertions"],
        metrics["deletions"],
        nm,
    )
    return {
        "qname": fields[0],
        "qlen": len(fields[9]) if fields[9] != "*" else 0,
        "qstart": 0,
        "qend": len(fields[9]) if fields[9] != "*" else 0,
        "strand": "-" if (flag & 0x10) else "+",
        "rname": rname,
        "rlen": rlen,
        "rstart": rstart,
        "rend": metrics["ref_end"],
        "mlen": matched,
        "blen": metrics["aligned_ref_bases"],
        "mapq": int(fields[4]),
        "tags": tags,
        "cigar": cigar,
        "nm": nm,
        "divergence": tags.get("de"),
        "is_unmapped": bool(flag & 0x4),
        "is_secondary": bool(flag & 0x100),
        "is_supplementary": bool(flag & 0x800),
    }


def _parse_mapping_file(raw_path, fmt, ref_names, ref_lengths, ref_index, nm_threshold=20, mapq_threshold=10, max_divergence=None, primary_only=True):
    total_mappings = np.zeros(len(ref_names), dtype=np.int64)
    high_quality_mappings = np.zeros(len(ref_names), dtype=np.int64)
    matched_sums = np.zeros(len(ref_names), dtype=np.float64)
    aligned_sums = np.zeros(len(ref_names), dtype=np.float64)
    coverage_events = defaultdict(list)
    oligomer_reads = defaultdict(list)
    sam_headers = []
    ref_lengths_by_name = {name: int(length) for name, length in zip(ref_names, ref_lengths)}

    with open(raw_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            if fmt == "sam" and line.startswith("@"):
                if line.startswith("@SQ"):
                    fields = line.rstrip("\n").split("\t")
                    sq = {}
                    for field in fields[1:]:
                        if ":" in field:
                            key, value = field.split(":", 1)
                            sq[key] = value
                    if "SN" in sq and "LN" in sq:
                        sam_headers.append((sq["SN"], int(sq["LN"])))
                        ref_lengths_by_name[sq["SN"]] = int(sq["LN"])
                continue

            if fmt == "legacy":
                rec = _parse_legacy_line(line)
            elif fmt == "paf":
                rec = _parse_paf_line(line)
            elif fmt == "sam":
                rec = _parse_sam_line(line, ref_lengths_by_name)
            else:
                raise ValueError(f"unsupported output format: {fmt}")

            if rec is None or rec["is_unmapped"]:
                continue
            if rec["rname"] not in ref_index:
                continue
            if primary_only and rec["is_secondary"]:
                continue

            ref_id = ref_index[rec["rname"]]
            total_mappings[ref_id] += 1

            high_quality = True
            if rec["nm"] is not None and int(rec["nm"]) > nm_threshold:
                high_quality = False
            if int(rec["mapq"]) < mapq_threshold:
                high_quality = False
            if max_divergence is not None and rec["divergence"] is not None and float(rec["divergence"]) > max_divergence:
                high_quality = False

            if high_quality:
                high_quality_mappings[ref_id] += 1
                metrics = _cigar_metrics(rec["cigar"], rec["rstart"])
                if metrics["aligned_ref_bases"] == 0 and rec["rend"] > rec["rstart"]:
                    metrics["intervals"] = [(rec["rstart"], rec["rend"])]
                    metrics["aligned_ref_bases"] = rec["rend"] - rec["rstart"]
                matched_sums[ref_id] += float(rec["mlen"])
                aligned_sums[ref_id] += float(metrics["aligned_ref_bases"])
                for start, end in metrics["intervals"]:
                    if end <= start:
                        continue
                    coverage_events[ref_id].append((int(start), 1))
                    coverage_events[ref_id].append((int(end), -1))

            if not rec["is_secondary"]:
                oligomer_reads[rec["qname"]].append(rec["rname"])

    oligomer_counts = Counter()
    for refs in oligomer_reads.values():
        if not refs:
            continue
        counts = Counter(refs)
        for rname, count in counts.items():
            if count >= 2:
                oligomer_counts[(rname, count)] += 1
        unique_refs = tuple(sorted(set(refs)))
        if len(unique_refs) > 1:
            oligomer_counts[("chimeric_multi_ref", unique_refs)] += 1

    return {
        "total_mappings": total_mappings,
        "high_quality_mappings": high_quality_mappings,
        "matched_sums": matched_sums,
        "aligned_sums": aligned_sums,
        "coverage_events": coverage_events,
        "oligomer_counts": oligomer_counts,
        "sam_headers": sam_headers,
    }


def _coverage_from_events(events, ref_len):
    if ref_len <= 0:
        return 0, 0.0
    if not events:
        return 0, 0.0
    position_deltas = defaultdict(int)
    for pos, delta in events:
        if pos < 0:
            pos = 0
        if pos > ref_len:
            pos = ref_len
        position_deltas[pos] += delta

    covered = 0
    area = 0.0
    current_depth = 0
    prev = 0
    for pos in sorted(position_deltas):
        if pos > prev and current_depth > 0:
            span = pos - prev
            covered += span
            area += current_depth * span
        current_depth += position_deltas[pos]
        prev = pos
    return covered, area


def _write_analysis_outputs(out_dir, ref_names, ref_lengths, metrics):
    out_dir.mkdir(parents=True, exist_ok=True)
    offsets = np.zeros(len(ref_names) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(np.array(ref_lengths, dtype=np.int64))
    total_mappings = metrics["total_mappings"]
    high_quality_mappings = metrics["high_quality_mappings"]
    matched_sums = metrics["matched_sums"]
    aligned_sums = metrics["aligned_sums"]

    matching_percent = np.full(len(ref_names), np.nan, dtype=np.float64)
    coverage_percent = np.zeros(len(ref_names), dtype=np.float64)
    mean_depth = np.zeros(len(ref_names), dtype=np.float64)

    for i, ref_len in enumerate(ref_lengths):
        covered, area = _coverage_from_events(metrics["coverage_events"].get(i, []), int(ref_len))
        if aligned_sums[i] > 0:
            matching_percent[i] = 100.0 * matched_sums[i] / aligned_sums[i]
        coverage_percent[i] = 100.0 * covered / float(ref_len) if ref_len else 0.0
        mean_depth[i] = area / float(ref_len) if ref_len else 0.0

    np.save(out_dir / "ref_names.npy", np.array(ref_names, dtype=object))
    np.save(out_dir / "ref_lengths.npy", np.array(ref_lengths, dtype=np.int64))
    np.save(out_dir / "offsets.npy", offsets)
    np.save(out_dir / "total_mappings.npy", total_mappings)
    np.save(out_dir / "high_quality_mappings.npy", high_quality_mappings)
    np.save(out_dir / "matching_percent.npy", matching_percent)
    np.save(out_dir / "coverage_percent.npy", coverage_percent)
    np.save(out_dir / "mean_depth.npy", mean_depth)

    with open(out_dir / "summary.tsv", "w", encoding="utf-8") as handle:
        handle.write(
            "ref_name\tref_length\ttotal_mappings\t"
            "high_quality_mappings\tmatching_percent\t"
            "coverage_percent\tmean_depth\n"
        )
        for i, name in enumerate(ref_names):
            match_value = matching_percent[i]
            match_str = f"{match_value:.4f}" if np.isfinite(match_value) else "nan"
            handle.write(
                f"{name}\t{int(ref_lengths[i])}\t{int(total_mappings[i])}\t"
                f"{int(high_quality_mappings[i])}\t{match_str}\t"
                f"{coverage_percent[i]:.4f}\t{mean_depth[i]:.4f}\n"
            )

    with open(out_dir / "oligomers.tsv", "w", encoding="utf-8") as handle:
        handle.write("event\tcount\n")
        for event, count in metrics["oligomer_counts"].items():
            handle.write(f"{event}\t{count}\n")

    return {
        "ref_names": ref_names,
        "ref_lengths": ref_lengths,
        "offsets": offsets,
        "total_mappings": total_mappings,
        "high_quality_mappings": high_quality_mappings,
        "matching_percent": matching_percent,
        "coverage_percent": coverage_percent,
        "mean_depth": mean_depth,
        "oligomer_counts": metrics["oligomer_counts"],
    }


def analyze_mappings(raw_path, fmt, ref_path, aligner, out_dir, nm_threshold=20, mapq_threshold=10, max_divergence=None, primary_only=True):
    provisional_ref_names, provisional_ref_lengths, provisional_ref_index = _load_reference_metadata(ref_path, aligner, fmt)
    parsed = _parse_mapping_file(
        raw_path=raw_path,
        fmt=fmt,
        ref_names=provisional_ref_names,
        ref_lengths=provisional_ref_lengths,
        ref_index=provisional_ref_index,
        nm_threshold=nm_threshold,
        mapq_threshold=mapq_threshold,
        max_divergence=max_divergence,
        primary_only=primary_only,
    )
    if fmt == "sam" and parsed["sam_headers"]:
        ref_names, ref_lengths, ref_index = _load_reference_metadata(ref_path, aligner, fmt, raw_headers=parsed["sam_headers"])
        parsed = _parse_mapping_file(
            raw_path=raw_path,
            fmt=fmt,
            ref_names=ref_names,
            ref_lengths=ref_lengths,
            ref_index=ref_index,
            nm_threshold=nm_threshold,
            mapq_threshold=mapq_threshold,
            max_divergence=max_divergence,
            primary_only=primary_only,
        )
    else:
        ref_names, ref_lengths = provisional_ref_names, provisional_ref_lengths
    return _write_analysis_outputs(Path(out_dir), ref_names, ref_lengths, parsed)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="minimap2.py",
        description="Map query sequences with the Python mappy interface using shared-index multithreaded batch mapping, then parse the mappings into summary arrays and tables.")
    parser.add_argument("ref", help="reference FASTA/FASTQ or prebuilt .mmi index")
    parser.add_argument("query", help="query FASTA/FASTQ (use '-' for stdin when supported by minimap2)")
    parser.add_argument("-x", dest="preset", help="preset: sr, map-pb, map-ont, asm5, asm10 or splice")
    parser.add_argument("-n", dest="min_cnt", type=int, help="minimum number of minimizers")
    parser.add_argument("-m", dest="min_chain_score", type=int, help="minimum chaining score")
    parser.add_argument("-k", dest="k", type=int, help="k-mer length")
    parser.add_argument("-w", dest="w", type=int, help="minimizer window length")
    parser.add_argument("-r", dest="bw", type=int, help="band width")
    parser.add_argument("-t", "--threads", dest="n_threads", type=int, default=3, help="mapping worker threads sharing one loaded index [3]")
    parser.add_argument("-o", "--output", dest="output", default=".", help="output directory for parsed .npy/.tsv files [.]")
    parser.add_argument("--raw-output", dest="raw_output", default=None, help="optional path to retain the raw mapping output stream")
    parser.add_argument("--output-format", choices=("legacy", "paf", "sam"), default="legacy", help="raw mapping format used before parsing [legacy]")
    parser.add_argument("-c", dest="out_cs", action="store_true", help="output the cs tag in raw mappings when supported")
    parser.add_argument("-M", dest="out_MD", action="store_true", help="output the MD tag in raw mappings when supported")
    parser.add_argument("--nm-threshold", dest="nm_threshold", type=int, default=20, help="maximum NM for a high-quality mapping [20]")
    parser.add_argument("--mapq-threshold", dest="mapq_threshold", type=int, default=10, help="minimum MAPQ for a high-quality mapping [10]")
    parser.add_argument("--max-divergence", dest="max_divergence", type=float, default=None, help="optional maximum divergence for a high-quality mapping")
    parser.add_argument("--primary-only", nargs="?", const=True, default=True, type=_bool_flag, help="count only non-secondary alignments in the summary metrics [yes]")
    parser.add_argument(
        "--compact-repeats",
        nargs="?",
        const=True,
        default=False,
        type=_bool_flag,
        help="build a compactified reference index; accepts optional yes/no and is ignored when loading a prebuilt .mmi")
    parser.add_argument("--compact-k", dest="compact_k", type=int, default=27, help="k-mer size for compactification [27]")
    parser.add_argument("--compact-ratio", dest="compact_ratio", type=float, default=0.20, help="remove k-mers with counts > int(ratio*n_seq)+1 during compactification [0.20]")
    parser.add_argument(
        "--verbose",
        nargs="?",
        const=True,
        default=True,
        type=_bool_flag,
        help="console logging only: yes prints per-progress batch updates, no prints summary progress updates [yes]")
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = build_parser().parse_args(argv[1:])
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    aligner = mp.Aligner(
        args.ref,
        preset=args.preset,
        min_cnt=args.min_cnt,
        min_chain_score=args.min_chain_score,
        k=args.k,
        w=args.w,
        bw=args.bw,
        n_threads=args.n_threads,
        compact_repeats=args.compact_repeats,
        compact_k=args.compact_k,
        compact_ratio=args.compact_ratio)
    if not aligner:
        raise Exception("ERROR: failed to load/build index file '{}'".format(args.ref))

    total_reads = 0
    progress_every = 1000
    if args.query != '-':
        total_reads = _count_reads(args.query)
        progress_every = 1000 if total_reads >= 10000 else max(1, total_reads // 10 if total_reads > 0 else 1)
        print(
            f"[mappy] starting threaded mapping for {total_reads} reads with {args.n_threads} threads and update interval {progress_every}",
            file=sys.stderr,
            flush=True,
        )
    else:
        print(
            f"[mappy] starting threaded mapping from stdin with {args.n_threads} threads; ETA unavailable without a pre-count",
            file=sys.stderr,
            flush=True,
        )

    remove_raw_after = False
    if args.raw_output:
        raw_path = Path(args.raw_output)
        raw_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        fd, temp_name = tempfile.mkstemp(prefix=f"minimap2_{args.output_format}_", suffix=f".{args.output_format}", dir=str(output_dir))
        os.close(fd)
        raw_path = Path(temp_name)
        remove_raw_after = True

    try:
        aligner.map_file(
            args.query,
            output_path=str(raw_path),
            output_format=args.output_format,
            n_threads=args.n_threads,
            cs=args.out_cs,
            MD=args.out_MD,
            verbose=args.verbose,
            total_reads=total_reads,
            progress_every=progress_every)

        analyze_mappings(
            raw_path=str(raw_path),
            fmt=args.output_format,
            ref_path=args.ref,
            aligner=aligner,
            out_dir=output_dir,
            nm_threshold=args.nm_threshold,
            mapq_threshold=args.mapq_threshold,
            max_divergence=args.max_divergence,
            primary_only=args.primary_only,
        )
    finally:
        if remove_raw_after and raw_path.exists():
            raw_path.unlink()

    print(
        f"[mappy] wrote parsed outputs to {output_dir}: {', '.join(DEFAULT_OUTPUT_FILES)}",
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main(sys.argv)
