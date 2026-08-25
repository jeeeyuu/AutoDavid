#!/usr/bin/env bash

# AutoDavid Python 스크립트 실행 예시입니다.
# 필요한 작업의 변수값을 수정한 뒤, 해당 python 명령의 주석(#)을 제거하여
# 터미널에 붙여 넣거나 이 파일에서 실행하세요.

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="python3"


# -----------------------------------------------------------------------------
# 1. Gene ID conversion
# -----------------------------------------------------------------------------
CONVERSION_INPUT="example.txt"
CONVERSION_OUTPUT="example_conversion.tsv"
CONVERSION_SPECIES="Homo sapiens"

# "$PYTHON_BIN" "$SCRIPT_DIR/auto_david_gene_id_conversion_batch.py" \
#   --input "$CONVERSION_INPUT" \
#   --output "$CONVERSION_OUTPUT" \
#   --species "$CONVERSION_SPECIES"


# -----------------------------------------------------------------------------
# 2. GO BP FAT + KEGG enrichment
# -----------------------------------------------------------------------------
ENRICHMENT_INPUT="example.txt"
ENRICHMENT_OUTPUT="example_result"
DAVID_EMAIL="your_registered_email@example.org"
ID_TYPE="auto"                  # auto, ensembl, entrez, symbol
ENRICHMENT_SPECIES="Homo sapiens"
EASE_CUTOFF="1.0"
MINIMUM_GENE_COUNT="1"

# "$PYTHON_BIN" "$SCRIPT_DIR/auto_david_webservice_enrichment.py" \
#   --input "$ENRICHMENT_INPUT" \
#   --output "$ENRICHMENT_OUTPUT" \
#   --email "$DAVID_EMAIL" \
#   --id-type "$ID_TYPE" \
#   --species "$ENRICHMENT_SPECIES" \
#   --ease "$EASE_CUTOFF" \
#   --count "$MINIMUM_GENE_COUNT"


# -----------------------------------------------------------------------------
# 3. Enrichment 결과 병합
# -----------------------------------------------------------------------------
MERGE_INPUTS=(
  "group1.go_bp_fat.tsv"
  "group2.go_bp_fat.tsv"
  "group3.go_bp_fat.csv"
)
MERGE_OUTPUT="merged_go_bp_fat.tsv"

# "$PYTHON_BIN" "$SCRIPT_DIR/auto_david_merge_enrichment_tables.py" \
#   --input "${MERGE_INPUTS[@]}" \
#   --output "$MERGE_OUTPUT"
