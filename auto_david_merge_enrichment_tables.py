#!/usr/bin/env python3
"""Merge DAVID enrichment TSV/CSV files side-by-side by Category + Term."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


# 특정 파일에 term이 없을 때 채울 값(R enrichment merge 로직과 동일)
MISSING_COUNT = "0"
MISSING_PVALUE = "1"
MISSING_GENES = "#N/A"


SCRIPT_DIR = Path(__file__).resolve().parent
REQUIRED_COLUMNS = ("Category", "Term", "Count", "PValue", "Genes")


class MergeError(RuntimeError):
    """Raised for malformed input or ambiguous enrichment records."""


@dataclass(frozen=True)
class EnrichmentRecord:
    category: str
    term: str
    count: str
    pvalue: str
    genes: str


def resolve_path(value: str | Path, base_dir: Path = SCRIPT_DIR) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else base_dir / path


def canonical_column(value: str) -> str:
    """Normalize harmless header spelling differences such as P-Value."""
    return re.sub(r"[^a-z0-9]", "", value.casefold())


def detect_delimiter(path: Path, sample: str) -> str:
    suffix = path.suffix.casefold()
    if suffix == ".csv":
        return ","
    if suffix in {".tsv", ".txt"}:
        return "\t"
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,").delimiter
    except csv.Error as exc:
        raise MergeError(
            f"구분자를 판별할 수 없습니다(.tsv 또는 .csv 권장): {path}"
        ) from exc


def read_enrichment_file(path: Path) -> list[EnrichmentRecord]:
    if not path.is_file():
        raise MergeError(f"입력 파일이 없습니다: {path}")
    try:
        text = path.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError) as exc:
        raise MergeError(f"입력 파일을 읽을 수 없습니다: {path}: {exc}") from exc
    if not text.strip():
        raise MergeError(f"입력 파일이 비어 있습니다: {path}")

    delimiter = detect_delimiter(path, text[:8192])
    rows = list(csv.reader(text.splitlines(), delimiter=delimiter))
    header_index = next(
        (
            index
            for index, row in enumerate(rows)
            if row and canonical_column(row[0]) == "category"
        ),
        None,
    )
    if header_index is None:
        raise MergeError(f"'Category' header row가 없습니다: {path.name}")

    header = rows[header_index]
    positions = {canonical_column(name): index for index, name in enumerate(header)}
    missing = [
        name
        for name in REQUIRED_COLUMNS
        if canonical_column(name) not in positions
    ]
    if missing:
        raise MergeError(
            f"{path.name}에 필요한 열이 없습니다: {', '.join(missing)}"
        )

    def value(row: list[str], column: str) -> str:
        index = positions[canonical_column(column)]
        return row[index].strip() if index < len(row) else ""

    records: list[EnrichmentRecord] = []
    seen: set[tuple[str, str]] = set()
    for row_number, row in enumerate(rows[header_index + 1 :], start=header_index + 2):
        if not row or not any(cell.strip() for cell in row):
            continue
        category = value(row, "Category")
        term = value(row, "Term")
        if not category and not term:
            continue
        if not category or not term:
            raise MergeError(
                f"{path.name}:{row_number} Category 또는 Term이 비어 있습니다."
            )
        key = (category, term)
        if key in seen:
            raise MergeError(
                f"{path.name}:{row_number}에 중복 Category+Term이 있습니다: "
                f"{category} / {term}"
            )
        seen.add(key)
        records.append(
            EnrichmentRecord(
                category=category,
                term=term,
                count=value(row, "Count"),
                pvalue=value(row, "PValue"),
                genes=value(row, "Genes"),
            )
        )
    return records


def merge_files(input_paths: list[Path], output_path: Path) -> tuple[int, int]:
    if not input_paths:
        raise MergeError("입력 파일 목록이 비어 있습니다.")
    duplicate_paths = [
        str(path)
        for index, path in enumerate(input_paths)
        if path in input_paths[:index]
    ]
    if duplicate_paths:
        raise MergeError(f"입력 목록에 중복 파일이 있습니다: {duplicate_paths[0]}")

    tables = [read_enrichment_file(path) for path in input_paths]
    lookups = [
        {(record.category, record.term): record for record in table}
        for table in tables
    ]

    ordered_keys: list[tuple[str, str]] = []
    seen_keys: set[tuple[str, str]] = set()
    for table in tables:
        for record in table:
            key = (record.category, record.term)
            if key not in seen_keys:
                seen_keys.add(key)
                ordered_keys.append(key)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            file_header = ["", ""]
            column_header = ["Category", "Term"]
            for index, path in enumerate(input_paths, start=1):
                file_header.extend([path.name, "", ""])
                column_header.extend(
                    [f"Count{index}", f"PValue{index}", f"Genes{index}"]
                )
            writer.writerow(file_header)
            writer.writerow(column_header)

            for category, term in ordered_keys:
                output_row = [category, term]
                for lookup in lookups:
                    record = lookup.get((category, term))
                    if record is None:
                        output_row.extend(
                            [MISSING_COUNT, MISSING_PVALUE, MISSING_GENES]
                        )
                    else:
                        output_row.extend(
                            [
                                record.count or MISSING_COUNT,
                                record.pvalue or MISSING_PVALUE,
                                record.genes or MISSING_GENES,
                            ]
                        )
                writer.writerow(output_row)
    except OSError as exc:
        raise MergeError(f"출력 파일을 쓸 수 없습니다: {output_path}: {exc}") from exc
    return len(ordered_keys), len(input_paths)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "DAVID enrichment TSV/CSV를 Category+Term union 기준으로 "
            "입력한 파일 순서대로 가로 병합합니다."
        )
    )
    parser.add_argument(
        "--input",
        dest="inputs",
        nargs="+",
        required=True,
        type=Path,
        help="입력 TSV/CSV 목록(입력한 순서 유지)",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="출력 TSV",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_paths = [resolve_path(path, Path.cwd()) for path in args.inputs]
    output_path = resolve_path(args.output, Path.cwd())
    try:
        term_count, file_count = merge_files(input_paths, output_path)
    except MergeError as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1
    print(f"완료: {output_path} ({term_count} terms x {file_count} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
