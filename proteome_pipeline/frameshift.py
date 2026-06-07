from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
import subprocess

from .fasta import FastaRecord, write_fasta
from .genetics import reverse_complement, translate_until_stop
from .gtf import load_id_set, load_transcripts
from .vcf import iter_frameshift_variants

GenomeFetcher = Callable[[str, int, int], str]


def samtools_fetcher(genome_path: str | Path) -> GenomeFetcher:
    genome = str(genome_path)

    def fetch(chrom: str, start: int, end: int) -> str:
        region = f"{chrom}:{start}-{end}"
        result = subprocess.check_output(["samtools", "faidx", genome, region], text=True)
        return "".join(result.splitlines()[1:]).upper()

    return fetch


def generate_frameshift_records(
    vcf_path: str | Path,
    gtf_path: str | Path,
    allowed_ids_path: str | Path,
    fetch_seq: GenomeFetcher,
    min_aa: int = 6,
) -> list[FastaRecord]:
    allowed = load_id_set(allowed_ids_path, strip_version=True)
    transcripts = load_transcripts(gtf_path, allowed)
    by_chrom: dict[str, list] = {}
    for tx in transcripts.values():
        by_chrom.setdefault(tx.chrom, []).append(tx)

    seen_by_gene: dict[str, set[str]] = {}
    records: list[FastaRecord] = []

    for variant in iter_frameshift_variants(vcf_path):
        for tx in by_chrom.get(variant.chrom, []):
            exons = sorted(tx.exons)
            if not any(start <= variant.pos <= end for start, end in exons):
                continue

            wt_cdna_parts: list[str] = []
            mut_cdna_parts: list[str] = []
            for start, end in exons:
                seq = fetch_seq(tx.chrom, start, end)
                wt_cdna_parts.append(seq)
                if start <= variant.pos <= end:
                    offset = variant.pos - start
                    mut_cdna_parts.append(seq[:offset] + variant.alt + seq[offset + len(variant.ref):])
                else:
                    mut_cdna_parts.append(seq)

            wt_cdna = "".join(wt_cdna_parts)
            mut_cdna = "".join(mut_cdna_parts)
            if tx.strand == "-":
                wt_cdna = reverse_complement(wt_cdna)
                mut_cdna = reverse_complement(mut_cdna)

            wt_protein = translate_until_stop(wt_cdna)
            mut_protein = translate_until_stop(mut_cdna, min_aa=min_aa)
            if not mut_protein:
                continue

            seen = seen_by_gene.setdefault(tx.gene_id, set())
            if mut_protein in seen:
                continue
            seen.add(mut_protein)

            shift = len(variant.alt) - len(variant.ref)
            kind = "DEL" if shift < 0 else "INS"
            header = (
                f">FRAMESHIFT_{tx.gene_id}_{tx.transcript_id} "
                f"type={kind} shift={shift}bp WT={len(wt_protein)} "
                f"MUT={len(mut_protein)} diff={len(mut_protein) - len(wt_protein)}"
            )
            records.append(FastaRecord(header, mut_protein))
    return records


def generate_frameshift_fasta(
    vcf_path: str | Path,
    gtf_path: str | Path,
    genome_path: str | Path,
    allowed_ids_path: str | Path,
    output_path: str | Path,
    min_aa: int = 6,
) -> int:
    records = generate_frameshift_records(
        vcf_path,
        gtf_path,
        allowed_ids_path,
        samtools_fetcher(genome_path),
        min_aa=min_aa,
    )
    return write_fasta(records, output_path)
