#!/usr/bin/env python3
"""Upload one complete gene list and run GO BP FAT plus KEGG via DAVID-WS."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

from functions.gene_id_conversion_functions import (
    DavidApiError,
    prepare_analysis_ids,
    read_analysis_ids,
)
from functions.enrichment_webservice_functions import (
    MAX_WEBSERVICE_GENES,
    DavidWebServiceClient,
    write_chart_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "gene list 전체를 DAVID Web Service에 한 번 업로드하고 "
            "GO BP FAT와 KEGG 결과를 각각 TSV로 저장합니다."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="gene ID 목록",
    )
    parser.add_argument(
        "--output",
        dest="output_prefix",
        required=True,
        type=Path,
        help="출력 prefix(.go_bp_fat.tsv/.kegg_pathway.tsv 자동 추가)",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("DAVID_EMAIL"),
        help="등록한 DAVID 이메일(또는 DAVID_EMAIL 환경변수)",
    )
    parser.add_argument(
        "--id-type",
        choices=("auto", "ensembl", "entrez", "symbol"),
        default="auto",
    )
    parser.add_argument(
        "--species",
        required=True,
        help='DAVID에서 선택할 species의 학명(예: "Homo sapiens")',
    )
    parser.add_argument("--keep-version", action="store_true")
    parser.add_argument(
        "--ease",
        type=float,
        default=1.0,
        help="chart report 최대 EASE score/P-value(기본값: 1.0)",
    )
    parser.add_argument(
        "--count", type=int, default=1, help="최소 gene count(기본값: 1)"
    )
    parser.add_argument("--list-name", default="uploaded_gene_list")
    parser.add_argument("--timeout", type=float, default=600.0)
    args = parser.parse_args()
    args.species = args.species.strip()
    if not args.species:
        parser.error("--species에는 비어 있지 않은 종명을 입력해야 합니다.")
    if not args.email:
        parser.error("--email 또는 DAVID_EMAIL 환경변수가 필요합니다.")
    if not (0.0 <= args.ease <= 1.0):
        parser.error("--ease는 0과 1 사이여야 합니다.")
    if args.count < 1:
        parser.error("--count는 1 이상이어야 합니다.")
    return args


def output_paths(prefix: Path) -> tuple[Path, Path]:
    base = re.sub(r"\.tsv$", "", str(prefix), flags=re.I)
    return Path(f"{base}.go_bp_fat.tsv"), Path(f"{base}.kegg_pathway.tsv")


def main() -> int:
    args = parse_args()
    try:
        inputs, detected_type = read_analysis_ids(
            args.input,
            args.id_type,
            args.keep_version,
            enforce_total_limit=False,
        )
        if len(inputs) > MAX_WEBSERVICE_GENES:
            raise DavidApiError(
                f"DAVID Web Service 권장 범위({MAX_WEBSERVICE_GENES} genes)를 "
                f"초과했습니다: {len(inputs)}개"
            )
        mygene_species = args.species.casefold().replace(" ", "_")
        ids, submitted_type, mapped_count = prepare_analysis_ids(
            inputs,
            detected_type,
            species=mygene_species,
            timeout=args.timeout,
            max_genes=MAX_WEBSERVICE_GENES,
        )
        if detected_type == "GENE_SYMBOL":
            print(
                f"gene symbols: {mapped_count}/{len(inputs)} IDs를 "
                f"Entrez로 변환했습니다(species={args.species}).",
                file=sys.stderr,
            )

        client = DavidWebServiceClient(args.timeout)
        client.authenticate(args.email)
        print(f"DAVID 인증 완료: {args.email}", file=sys.stderr)
        mapped_fraction = client.add_list(ids, submitted_type, args.list_name)
        print(
            f"목록 업로드 완료: {len(ids)} IDs, "
            f"DAVID mapping={mapped_fraction * 100:.2f}%",
            file=sys.stderr,
        )
        selected_species = client.select_species(args.species)
        print(f"species 선택: {selected_species}", file=sys.stderr)

        go_path, kegg_path = output_paths(args.output_prefix)
        client.set_categories("GOTERM_BP_FAT")
        go_records = client.get_chart_report(args.ease, args.count)
        write_chart_tsv(go_path, go_records)
        print(f"GO BP FAT 완료: {go_path} ({len(go_records)} records)", file=sys.stderr)

        client.set_categories("KEGG_PATHWAY")
        kegg_records = client.get_chart_report(args.ease, args.count)
        write_chart_tsv(kegg_path, kegg_records)
        print(f"KEGG 완료: {kegg_path} ({len(kegg_records)} records)", file=sys.stderr)
    except DavidApiError as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1

    print(f"완료: {go_path}, {kegg_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
