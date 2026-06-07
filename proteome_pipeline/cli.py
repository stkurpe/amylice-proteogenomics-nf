from __future__ import annotations

import argparse
from pathlib import Path

from .fasta import clean_stop_codons, combine_fastas, read_fasta, write_fasta
from .frameshift import generate_frameshift_fasta
from .gtf import clean_transcript_ids
from .nonsense import write_nonsense_report
from .verify import write_verification_report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="proteome-pipeline")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("clean-ids")
    p.add_argument("--gtf", required=True)
    p.add_argument("--abundance")
    p.add_argument("--out", required=True)
    p.add_argument("--min-tpm", type=float, default=1.0)
    p.add_argument("--min-cds-bp", type=int, default=10)

    p = sub.add_parser("clean-proteins")
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)

    p = sub.add_parser("frameshifts")
    p.add_argument("--vcf", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--ids", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--min-aa", type=int, default=6)

    p = sub.add_parser("combine")
    p.add_argument("--out", required=True)
    p.add_argument("fastas", nargs="+")

    p = sub.add_parser("nonsense-report")
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)

    p = sub.add_parser("verify")
    p.add_argument("--proteo-dir", required=True)
    p.add_argument("--out", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "clean-ids":
        ids = clean_transcript_ids(args.gtf, args.abundance, args.out, args.min_tpm, args.min_cds_bp)
        print(f"Wrote {len(ids)} clean transcript ids to {args.out}")
    elif args.command == "clean-proteins":
        records = clean_stop_codons(read_fasta(args.input))
        print(f"Wrote {write_fasta(records, args.out)} cleaned proteins to {args.out}")
    elif args.command == "frameshifts":
        count = generate_frameshift_fasta(args.vcf, args.gtf, args.genome, args.ids, args.out, args.min_aa)
        print(f"Wrote {count} unique frameshift proteins to {args.out}")
    elif args.command == "combine":
        records = combine_fastas([Path(p) for p in args.fastas])
        print(f"Wrote {write_fasta(records, args.out)} combined proteins to {args.out}")
    elif args.command == "nonsense-report":
        count = write_nonsense_report(args.input, args.out)
        print(f"Wrote {count} nonsense candidates to {args.out}")
    elif args.command == "verify":
        ok = write_verification_report(args.proteo_dir, args.out)
        print(f"Verification {'OK' if ok else 'requires attention'}: {args.out}")
        return 0 if ok else 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
