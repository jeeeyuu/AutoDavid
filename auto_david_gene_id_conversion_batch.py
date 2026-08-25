#!/usr/bin/env python3
"""Run URL-length-aware batch conversion for DAVID Ensembl gene IDs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from functions.gene_id_conversion_functions import (
    DavidApiClient,
    DavidApiError,
    GeneMapping,
    fetch_entrez_mappings,
    parse_gene_report,
    read_ensembl_ids,
    split_ids_for_api,
    unique_query_ids,
    write_conversion_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Ensembl gene ID 목록을 DAVID URL 길이에 맞게 자동 분할하여 "
            "gene ID conversion 결과를 하나의 TSV로 합칩니다."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Ensembl ID 목록(공백/쉼표 구분 가능)",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="통합 출력 TSV",
    )
    parser.add_argument(
        "--keep-version",
        action="store_true",
        help="Ensembl ID 끝의 버전(.숫자)을 제거하지 않음",
    )
    parser.add_argument(
        "--species",
        required=True,
        help='매핑할 species의 학명(예: "Homo sapiens")',
    )
    parser.add_argument("--timeout", type=float, default=120.0, help="HTTP timeout 초")
    args = parser.parse_args()
    args.species = args.species.strip()
    if not args.species:
        parser.error("--species에는 비어 있지 않은 종명을 입력해야 합니다.")
    return args


def main() -> int:
    args = parse_args()
    try:
        inputs = read_ensembl_ids(
            args.input,
            args.keep_version,
            enforce_total_limit=False,
        )
        chunks = split_ids_for_api(unique_query_ids(inputs), tool="geneReport")
        client = DavidApiClient(args.timeout)
        mappings: list[GeneMapping] = []

        print(
            f"{len(inputs)} unique IDs를 URL 길이에 따라 "
            f"{len(chunks)}개 요청으로 나눕니다.",
            file=sys.stderr,
        )
        for index, chunk in enumerate(chunks, start=1):
            print(
                f"[{index}/{len(chunks)}] {len(chunk)} IDs 요청 중...",
                file=sys.stderr,
            )
            page = client.run(chunk, tool="geneReport", allow_empty=True)
            if page:
                mappings.extend(parse_gene_report(page))

        print("Entrez Gene ID batch 변환 중...", file=sys.stderr)
        entrez_mappings = fetch_entrez_mappings(
            unique_query_ids(inputs),
            timeout=args.timeout,
        )
        write_conversion_tsv(
            args.output,
            inputs,
            mappings,
            expected_species=args.species,
            entrez_by_source=entrez_mappings,
        )
    except DavidApiError as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1

    selected_mappings = [
        mapping
        for mapping in mappings
        if mapping.species.casefold() == args.species.casefold()
    ]
    mapped_sources = {mapping.source_id for mapping in selected_mappings}
    entrez_sources = {gene_id for gene_id, values in entrez_mappings.items() if values}
    print(
        f"완료: {args.output} "
        f"({len(mapped_sources)}/{len(inputs)} unique IDs mapped, "
        f"{len(chunks)} DAVID API requests, "
        f"{len(entrez_sources)}/{len(inputs)} Entrez IDs, "
        f"species={args.species})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
