#!/usr/bin/env python

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import zipfile
from pathlib import Path

import mappy as mp

import minimap2

VALID_DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
GENBANK_SUFFIXES = {".gb", ".gbk", ".genbank"}
FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".ffn"}
RAW_OUTPUT_SUFFIX = {
    "legacy": "legacy.tsv",
    "paf": "paf",
    "sam": "sam",
}


def _normalize_reference_name(name):
    text = str(name).strip()
    if not text:
        raise ValueError("reference name is empty")
    text = re.sub(r"\s+", "_", text)
    return text


def _normalize_sequence(sequence, source_id):
    seq = "".join(str(sequence).split()).upper().replace("U", "T")
    if not seq:
        raise ValueError(f"{source_id}: sequence is empty")
    if not VALID_DNA_RE.fullmatch(seq):
        raise ValueError(f"{source_id}: sequence contains unsupported characters")
    return seq


def _resolve_output_path(path_value, output_dir, default_name):
    if path_value:
        path = Path(path_value)
        if not path.is_absolute():
            path = output_dir / path
        return path
    return output_dir / default_name


def _load_seqio():
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise SystemExit("Biopython is required to parse GenBank reference inputs") from exc
    return SeqIO


def _detect_reference_format(path, user_format):
    if user_format != "auto":
        return user_format
    suffix = path.suffix.lower()
    if suffix in FASTA_SUFFIXES:
        return "fasta"
    if suffix == ".csv":
        return "csv"
    if suffix in GENBANK_SUFFIXES:
        return "genbank"
    if suffix == ".zip":
        return "genbank-zip"
    if suffix in (".jsonl", ".json"):
        return "gsc-jsonl"
    raise ValueError(f"could not auto-detect reference format for {path}")


def _iter_fasta_records(path):
    SeqIO = _load_seqio()
    with open(path, "r", encoding="utf-8") as handle:
        records = SeqIO.parse(handle, "fasta")
        for record in records:
            name = record.description or record.id
            source_id = name
            seq = str(record.seq)
            yield {
                "name": _normalize_reference_name(name),
                "sequence": _normalize_sequence(seq, source_id),
                "source_format": "fasta",
                "source_id": source_id,
                "has_annotations": False,
                "metadata": {},
            }


def _iter_csv_records(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV reference file is missing a header row")
        missing = {"name", "sequence"} - set(reader.fieldnames)
        if missing:
            raise ValueError(f"CSV reference file is missing required columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, start=2):
            source_id = row.get("name") or f"csv_row_{line_number}"
            yield {
                "name": _normalize_reference_name(row["name"]),
                "sequence": _normalize_sequence(row["sequence"], source_id),
                "source_format": "csv",
                "source_id": source_id,
                "has_annotations": False,
                "metadata": {},
            }


def _iter_genbank_records(path, zip_member=None):
    SeqIO = _load_seqio()
    if zip_member is None:
        with open(path, "r", encoding="utf-8") as handle:
            records = SeqIO.parse(handle, "genbank")
            for record in records:
                record_id = record.id or record.name
                yield {
                    "name": _normalize_reference_name(record_id),
                    "sequence": _normalize_sequence(str(record.seq), record_id),
                    "source_format": "genbank",
                    "source_id": record_id,
                    "has_annotations": True,
                    "metadata": {
                        "record_id": record.id,
                        "record_name": record.name,
                        "description": record.description,
                    },
                }
        return

    member_name = Path(zip_member).stem
    with zipfile.ZipFile(path, "r") as zf:
        with zf.open(zip_member, "r") as raw_handle:
            with io.TextIOWrapper(raw_handle, encoding="utf-8") as handle:
                records = SeqIO.parse(handle, "genbank")
                for record in records:
                    record_id = record.id or record.name
                    source_id = f"{zip_member}:{record_id}"
                    yield {
                        "name": _normalize_reference_name(f"{member_name}__{record_id}"),
                        "sequence": _normalize_sequence(str(record.seq), source_id),
                        "source_format": "genbank-zip",
                        "source_id": source_id,
                        "has_annotations": True,
                        "metadata": {
                            "record_id": record.id,
                            "record_name": record.name,
                            "description": record.description,
                            "zip_member": zip_member,
                        },
                    }


def _iter_genbank_zip_records(path):
    with zipfile.ZipFile(path, "r") as zf:
        members = [member for member in zf.namelist() if not member.endswith("/")]
    for member in members:
        member_suffix = Path(member).suffix.lower()
        if member_suffix not in GENBANK_SUFFIXES:
            raise ValueError(f"unsupported archive member in GenBank ZIP: {member}")
        yield from _iter_genbank_records(path, zip_member=member)


def _iter_gsc_jsonl_records(path):
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            data = json.loads(line)
            genetic_system = data.get("genetic_system") or {}
            parts = genetic_system.get("parts")
            if not isinstance(parts, list) or not parts:
                raise ValueError(f"systems.jsonl line {line_number}: missing genetic_system.parts")
            sequences = []
            for part_index, part in enumerate(parts, start=1):
                if "sequence" not in part:
                    raise ValueError(f"systems.jsonl line {line_number}: part {part_index} is missing sequence")
                sequences.append(str(part["sequence"]))
            system_info = data.get("system_info") or {}
            name = system_info.get("name") or genetic_system.get("id") or f"system_{line_number}"
            source_id = genetic_system.get("id") or name
            yield {
                "name": _normalize_reference_name(name),
                "sequence": _normalize_sequence("".join(sequences), source_id),
                "source_format": "gsc-jsonl",
                "source_id": source_id,
                "has_annotations": True,
                "metadata": {
                    "system_name": system_info.get("name"),
                    "genetic_system_id": genetic_system.get("id"),
                    "run_id": data.get("run_id"),
                    "part_count": len(parts),
                },
            }


def _collect_references(path, reference_format):
    if reference_format == "fasta":
        records = list(_iter_fasta_records(path))
    elif reference_format == "csv":
        records = list(_iter_csv_records(path))
    elif reference_format == "genbank":
        records = list(_iter_genbank_records(path))
    elif reference_format == "genbank-zip":
        records = list(_iter_genbank_zip_records(path))
    elif reference_format == "gsc-jsonl":
        records = list(_iter_gsc_jsonl_records(path))
    else:
        raise ValueError(f"unsupported reference format: {reference_format}")

    seen = set()
    for record in records:
        if record["name"] in seen:
            raise ValueError(f"duplicate normalized reference name: {record['name']}")
        seen.add(record["name"])
    if not records:
        raise ValueError("no reference records were found")
    return records


def _write_compiled_fasta(records, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    with open(path, "w", encoding="utf-8") as handle:
        for record in records:
            header = f">{record['name']}\n"
            handle.write(header)
            digest.update(header.encode("utf-8"))
            sequence = record["sequence"]
            for start in range(0, len(sequence), 80):
                chunk = sequence[start:start + 80]
                handle.write(chunk)
                handle.write("\n")
                digest.update(chunk.encode("utf-8"))
                digest.update(b"\n")
    return digest.hexdigest()


def _write_manifest(records, prefix, source_path, reference_format, compiled_fasta, compiled_index):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = Path(f"{prefix}.tsv")
    json_path = Path(f"{prefix}.json")

    rows = []
    for record in records:
        row = {
            "name": record["name"],
            "source_format": record["source_format"],
            "source_path": str(source_path),
            "source_id": record["source_id"],
            "sequence_length": len(record["sequence"]),
            "has_annotations": record["has_annotations"],
        }
        row.update(record["metadata"])
        rows.append(row)

    columns = [
        "name",
        "source_format",
        "source_path",
        "source_id",
        "sequence_length",
        "has_annotations",
        "record_id",
        "record_name",
        "description",
        "zip_member",
        "system_name",
        "genetic_system_id",
        "run_id",
        "part_count",
    ]

    with open(tsv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out_row = dict(row)
            out_row["has_annotations"] = "yes" if row["has_annotations"] else "no"
            writer.writerow(out_row)

    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "source_path": str(source_path),
                "reference_format": reference_format,
                "compiled_fasta": str(compiled_fasta),
                "compiled_index": str(compiled_index),
                "records": rows,
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
    return tsv_path, json_path


def _index_metadata_path(compiled_index):
    return compiled_index.with_suffix(compiled_index.suffix + ".meta.json")


def _expected_index_metadata(args, compiled_fasta, fasta_sha256):
    return {
        "compiled_fasta": str(compiled_fasta),
        "compiled_fasta_sha256": fasta_sha256,
        "preset": args.preset,
        "k": args.k,
        "w": args.w,
        "compact_repeats": bool(args.compact_repeats),
        "compact_k": args.compact_k,
        "compact_ratio": args.compact_ratio,
    }


def _load_index_metadata(path):
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_index_metadata(path, metadata):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _build_or_reuse_index(compiled_fasta, compiled_index, args, fasta_sha256):
    compiled_index.parent.mkdir(parents=True, exist_ok=True)
    meta_path = _index_metadata_path(compiled_index)
    expected = _expected_index_metadata(args, compiled_fasta, fasta_sha256)
    existing = _load_index_metadata(meta_path)
    if compiled_index.exists() and existing == expected:
        print(f"[AutoAlign] reusing compiled index {compiled_index}", file=sys.stderr, flush=True)
        return

    print(f"[AutoAlign] building compiled index {compiled_index}", file=sys.stderr, flush=True)
    aligner = mp.Aligner(
        str(compiled_fasta),
        preset=args.preset,
        min_cnt=args.min_cnt,
        min_chain_score=args.min_chain_score,
        k=args.k,
        w=args.w,
        bw=args.bw,
        n_threads=args.n_threads,
        fn_idx_out=str(compiled_index),
        compact_repeats=args.compact_repeats,
        compact_k=args.compact_k,
        compact_ratio=args.compact_ratio,
    )
    if not aligner:
        raise RuntimeError(f"failed to build index for {compiled_fasta}")
    _write_index_metadata(meta_path, expected)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="AutoAlign.py",
        description="Normalize references from FASTA, CSV, GenBank, GenBank ZIP, or GSC JSONL, then map long reads using the Python minimap2 interface.",
    )
    parser.add_argument("references", help="reference input file (FASTA, CSV, GenBank, GenBank ZIP, or GSC systems.jsonl)")
    parser.add_argument("query", help="query FASTA/FASTQ (use '-' for stdin when supported by minimap2)")
    parser.add_argument(
        "--reference-format",
        choices=("auto", "fasta", "csv", "genbank", "genbank-zip", "gsc-jsonl"),
        default="auto",
        help="reference input format [auto]",
    )
    parser.add_argument("--compiled-ref-fasta", dest="compiled_ref_fasta", default=None, help="optional path for the generated normalized reference FASTA")
    parser.add_argument("--compiled-index", dest="compiled_index", default=None, help="optional path for the generated/reused .mmi index")
    parser.add_argument("--manifest-prefix", dest="manifest_prefix", default=None, help="optional prefix for the reference manifest TSV/JSON outputs")
    minimap2.add_mapping_arguments(parser)
    return parser


def run_with_args(args):
    source_path = Path(args.references)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    compiled_fasta = _resolve_output_path(args.compiled_ref_fasta, output_dir, "references.compiled.fa")
    compiled_index = _resolve_output_path(args.compiled_index, output_dir, "references.compiled.mmi")
    manifest_prefix = _resolve_output_path(args.manifest_prefix, output_dir, "references.manifest")

    if args.raw_output is None:
        raw_name = f"mappings.{RAW_OUTPUT_SUFFIX[args.output_format]}"
        args.raw_output = str(output_dir / raw_name)

    reference_format = _detect_reference_format(source_path, args.reference_format)
    print(f"[AutoAlign] parsing references from {source_path} as {reference_format}", file=sys.stderr, flush=True)
    records = _collect_references(source_path, reference_format)
    fasta_sha256 = _write_compiled_fasta(records, compiled_fasta)
    manifest_tsv, manifest_json = _write_manifest(records, manifest_prefix, source_path, reference_format, compiled_fasta, compiled_index)
    print(
        f"[AutoAlign] wrote {len(records)} normalized references to {compiled_fasta} and manifests to {manifest_tsv} / {manifest_json}",
        file=sys.stderr,
        flush=True,
    )

    _build_or_reuse_index(compiled_fasta, compiled_index, args, fasta_sha256)

    args.ref = str(compiled_index)
    return minimap2.run_with_args(args)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = build_parser().parse_args(argv[1:])
    return run_with_args(args)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
