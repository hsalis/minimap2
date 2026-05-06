import pysam
import numpy as np
from collections import Counter
from pathlib import Path

def phred_weight(q):
    if q is None:
        return 1.0
    return 1.0 - 10 ** (-q / 10.0)

def get_tag_or_none(read, tag):
    try:
        return read.get_tag(tag)
    except KeyError:
        return None

def parse_sa_tag(sa_tag):
    if not sa_tag:
        return []

    out = []
    for item in sa_tag.rstrip(";").split(";"):
        if not item:
            continue
        fields = item.split(",")
        if len(fields) >= 6:
            out.append({
                "rname": fields[0],
                "pos": int(fields[1]),
                "strand": fields[2],
                "cigar": fields[3],
                "mapq": int(fields[4]),
                "nm": int(fields[5]),
            })
    return out


def process_alignment_correctness_with_reference(
    read,
    reference_fasta,
    correct_vec,
    aligned_vec,
    ref_offset,
    use_quality=True,
):
    """
    Uses the reference FASTA sequence to determine per-position correctness.

    Rules:
      M / = : compare query base to reference base
      X     : mismatch, incorrect
      D     : deletion, incorrect
      I     : insertion, neutral
      S/H   : ignored
      N     : reference skip, neutral
    """

    ref_name = read.reference_name
    ref_start = read.reference_start  # 0-based
    ref_end = read.reference_end      # 0-based exclusive

    ref_seq = reference_fasta.fetch(ref_name, ref_start, ref_end).upper()
    query_seq = read.query_sequence.upper()
    quals = read.query_qualities

    qpos = 0
    rpos = ref_start

    for op, length in read.cigartuples:
        # 0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X

        if op in (0, 7, 8):  # M, =, X
            for _ in range(length):
                idx = ref_offset + rpos
                local_ref_pos = rpos - ref_start

                ref_base = ref_seq[local_ref_pos]
                query_base = query_seq[qpos]

                aligned_vec[idx] += 1

                if op == 8:
                    # explicit mismatch
                    pass
                elif query_base == ref_base:
                    q = quals[qpos] if quals is not None else None
                    correct_vec[idx] += phred_weight(q) if use_quality else 1.0

                qpos += 1
                rpos += 1

        elif op == 1:  # insertion in read
            qpos += length

        elif op == 2:  # deletion from read relative to reference
            for _ in range(length):
                idx = ref_offset + rpos
                aligned_vec[idx] += 1
                rpos += 1

        elif op == 3:  # skipped reference region
            rpos += length

        elif op == 4:  # soft clip
            qpos += length

        elif op in (5, 6):  # hard clip or padding
            pass


def analyze_sam(
    sam_path,
    ref_fasta_path,
    out_dir,
    nm_threshold=20,
    mapq_threshold=10,
    max_divergence=None,
    use_quality=True,
    primary_only=True,
):
    sam_path = Path(sam_path)
    ref_fasta_path = Path(ref_fasta_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sam = pysam.AlignmentFile(str(sam_path), "r")
    reference_fasta = pysam.FastaFile(str(ref_fasta_path))

    ref_names = list(sam.references)
    ref_lengths = list(sam.lengths)

    n_refs = len(ref_names)

    offsets = np.zeros(n_refs + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(ref_lengths)
    total_ref_nt = int(offsets[-1])

    total_mappings = np.zeros(n_refs, dtype=np.int64)
    high_quality_mappings = np.zeros(n_refs, dtype=np.int64)

    correct_vec = np.memmap(
        out_dir / "correct_weighted_nt.dat",
        dtype=np.float32,
        mode="w+",
        shape=(total_ref_nt,),
    )

    aligned_vec = np.memmap(
        out_dir / "aligned_nt.dat",
        dtype=np.uint32,
        mode="w+",
        shape=(total_ref_nt,),
    )

    oligomer_counts = Counter()
    split_alignment_counts = Counter()

    for read in sam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        if primary_only and (read.is_secondary or read.is_supplementary):
            continue

        ref_id = read.reference_id
        ref_name = read.reference_name
        ref_offset = offsets[ref_id]

        total_mappings[ref_id] += 1

        nm = get_tag_or_none(read, "NM")
        divergence = get_tag_or_none(read, "de")

        high_quality = True

        if nm is not None and nm > nm_threshold:
            high_quality = False

        if read.mapping_quality < mapq_threshold:
            high_quality = False

        if max_divergence is not None:
            if divergence is not None and divergence > max_divergence:
                high_quality = False

        if high_quality:
            high_quality_mappings[ref_id] += 1

            process_alignment_correctness_with_reference(
                read=read,
                reference_fasta=reference_fasta,
                correct_vec=correct_vec,
                aligned_vec=aligned_vec,
                ref_offset=ref_offset,
                use_quality=use_quality,
            )

        sa_tag = get_tag_or_none(read, "SA")
        sa_entries = parse_sa_tag(sa_tag)

        if sa_entries:
            split_alignment_counts[ref_name] += 1

            all_refs = [ref_name] + [x["rname"] for x in sa_entries]
            c = Counter(all_refs)

            for r, count in c.items():
                if count >= 2:
                    oligomer_counts[(r, count)] += 1

            if len(set(all_refs)) > 1:
                oligomer_counts[("chimeric_multi_ref", tuple(all_refs))] += 1

    sam.close()
    reference_fasta.close()

    correct_vec.flush()
    aligned_vec.flush()

    matching_percent = np.zeros(n_refs, dtype=np.float64)
    coverage_percent = np.zeros(n_refs, dtype=np.float64)
    mean_depth = np.zeros(n_refs, dtype=np.float64)

    for i in range(n_refs):
        start = offsets[i]
        end = offsets[i + 1]

        ref_correct = correct_vec[start:end]
        ref_aligned = aligned_vec[start:end]

        correct_sum = float(ref_correct.sum())
        aligned_sum = int(ref_aligned.sum())

        covered_positions = int(np.count_nonzero(ref_aligned))
        ref_len = ref_lengths[i]

        if aligned_sum > 0:
            matching_percent[i] = 100.0 * correct_sum / aligned_sum
        else:
            matching_percent[i] = np.nan

        coverage_percent[i] = 100.0 * covered_positions / ref_len
        mean_depth[i] = aligned_sum / ref_len

    np.save(out_dir / "ref_names.npy", np.array(ref_names, dtype=object))
    np.save(out_dir / "ref_lengths.npy", np.array(ref_lengths, dtype=np.int64))
    np.save(out_dir / "offsets.npy", offsets)
    np.save(out_dir / "total_mappings.npy", total_mappings)
    np.save(out_dir / "high_quality_mappings.npy", high_quality_mappings)
    np.save(out_dir / "matching_percent.npy", matching_percent)
    np.save(out_dir / "coverage_percent.npy", coverage_percent)
    np.save(out_dir / "mean_depth.npy", mean_depth)

    with open(out_dir / "summary.tsv", "w") as f:
        f.write(
            "ref_name\tref_length\ttotal_mappings\t"
            "high_quality_mappings\tmatching_percent\t"
            "coverage_percent\tmean_depth\n"
        )

        for i, name in enumerate(ref_names):
            f.write(
                f"{name}\t"
                f"{ref_lengths[i]}\t"
                f"{total_mappings[i]}\t"
                f"{high_quality_mappings[i]}\t"
                f"{matching_percent[i]:.4f}\t"
                f"{coverage_percent[i]:.4f}\t"
                f"{mean_depth[i]:.4f}\n"
            )

    with open(out_dir / "oligomers.tsv", "w") as f:
        f.write("event\tcount\n")
        for event, count in oligomer_counts.items():
            f.write(f"{event}\t{count}\n")

    return {
        "ref_names": ref_names,
        "ref_lengths": ref_lengths,
        "total_mappings": total_mappings,
        "high_quality_mappings": high_quality_mappings,
        "matching_percent": matching_percent,
        "coverage_percent": coverage_percent,
        "mean_depth": mean_depth,
        "oligomer_counts": oligomer_counts,
    }


if __name__ == "__main__":
    results = analyze_sam(
        sam_path="results.sam",
        ref_fasta_path="compiled_references.fa",
        out_dir="sam_metrics",
        nm_threshold=20,
        mapq_threshold=10,
        max_divergence=0.03,
        use_quality=True,
        primary_only=True,
    )


