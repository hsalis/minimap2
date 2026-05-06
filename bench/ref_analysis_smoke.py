#!/usr/bin/env python3
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def mutate(seq, pos, base):
    return seq[:pos] + base + seq[pos + 1:]


def insert(seq, pos, ins):
    return seq[:pos] + ins + seq[pos:]


def delete(seq, pos):
    return seq[:pos] + seq[pos + 1:]


def write_fasta(path, records):
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def write_fastq(path, records):
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def read_fasta(path):
    out = {}
    name = None
    seq = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(seq)
            name = line[1:]
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        out[name] = "".join(seq)
    return out


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return {row["ref_name"]: row for row in csv.DictReader(handle, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description="Smoke-test per-reference consensus analysis.")
    ap.add_argument("--minimap2", required=True)
    args = ap.parse_args()

    refs = {
        "ref_mismatch": "ACGTTGCATGACCTGATCGTAGGCTAACGTTAGCGTACCTGATGGCATTCGATCGTACGGTACCTAGTCGAT",
        "ref_insertion": "TTGACCAAGTCCGATGATCGTACCGGATTCGATGGCATACCGTAGCTAGGCTTACGATCGGATCCGTAGTACC",
        "ref_deletion": "GGCATGTATCGGATCACGTAGCTAGTCCGATGGCATACGTTAGCCGATTCGATGGTACCATGCTAGTCGATGC",
        "ref_unmapped": "CATGGTACGATTCGAGTCCATGATCGTAGCTACCGTATGGCATCGATGACTGATCGTACCGATAGTCCGATTA",
    }

    reads = [
        ("mm_1", refs["ref_mismatch"]),
        ("mm_2", mutate(refs["ref_mismatch"], 20, "C")),
        ("mm_3", mutate(refs["ref_mismatch"], 20, "C")),
        ("ins_1", refs["ref_insertion"]),
        ("ins_2", insert(refs["ref_insertion"], 30, "G")),
        ("ins_3", insert(refs["ref_insertion"], 30, "G")),
        ("del_1", delete(refs["ref_deletion"], 41)),
        ("del_2", delete(refs["ref_deletion"], 41)),
        ("del_3", delete(refs["ref_deletion"], 41)),
    ]

    expected_consensus = {
        "ref_mismatch": mutate(refs["ref_mismatch"], 20, "C"),
        "ref_insertion": insert(refs["ref_insertion"], 30, "G"),
        "ref_deletion": delete(refs["ref_deletion"], 41),
    }
    expected_cigar = {
        "ref_mismatch": "20=1X51=",
        "ref_insertion": "30=1I43=",
        "ref_deletion": "41=1D31=",
    }

    with tempfile.TemporaryDirectory(prefix="mm-ref-analysis-") as tmpdir:
        tmpdir = Path(tmpdir)
        ref_path = tmpdir / "refs.fa"
        read_path = tmpdir / "reads.fq"
        prefix = tmpdir / "out"
        write_fasta(ref_path, refs.items())
        write_fastq(read_path, reads)

        cmd = [
            str(Path(args.minimap2).resolve()),
            "-a",
            "-x",
            "sr",
            "-k",
            "5",
            "-w",
            "3",
            "-n",
            "1",
            "-m",
            "5",
            "--secondary=no",
            "--ref-analysis-prefix",
            str(prefix),
            str(ref_path),
            str(read_path),
        ]
        proc = subprocess.run(cmd, text=True, capture_output=True)
        if proc.returncode != 0:
            raise SystemExit(f"minimap2 failed\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

        consensus = read_fasta(prefix.with_suffix(".consensus.fa"))
        stats = read_tsv(prefix.with_suffix(".reference_stats.tsv"))
        cigars = read_tsv(prefix.with_suffix(".consensus_cigar.tsv"))

        for ref_name, seq in expected_consensus.items():
            if consensus.get(ref_name) != seq:
                raise SystemExit(f"unexpected consensus for {ref_name}: {consensus.get(ref_name)!r} != {seq!r}")
            if cigars[ref_name]["cigar"] != expected_cigar[ref_name]:
                raise SystemExit(f"unexpected cigar for {ref_name}: {cigars[ref_name]['cigar']} != {expected_cigar[ref_name]}")
            if stats[ref_name]["primary_reads_total"] != "3":
                raise SystemExit(f"unexpected primary count for {ref_name}: {stats[ref_name]['primary_reads_total']}")
            if stats[ref_name]["primary_reads_used"] != "3":
                raise SystemExit(f"unexpected used primary count for {ref_name}: {stats[ref_name]['primary_reads_used']}")

        if "ref_unmapped" in consensus:
            raise SystemExit("unmapped reference should not be emitted to consensus FASTA")
        if stats["ref_unmapped"]["primary_reads_total"] != "0":
            raise SystemExit("unmapped reference should have zero primary reads")

        print("ref-analysis smoke test passed")


if __name__ == "__main__":
    main()
