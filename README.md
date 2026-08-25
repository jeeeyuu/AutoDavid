# AutoDavid

**AutoDavid**는 DAVID ID 변환, 대용량 enrichment, 여러 enrichment 결과 비교 병합을 위한 독립 프로그램입니다. Python 3 표준 라이브러리만 사용하므로 별도의 `pip install`은 필요하지 않습니다.

## 이 패키지를 사용하는 이유

DAVID linking API(`api.jsp`)는 ID를 GET URL에 넣으므로 전체 URL 약 2,048자와 light-duty 목록 약 400 genes 제한을 받습니다. Ensembl ID는 길어서 실제로는 400개보다 적은 수에서도 URL 제한에 걸릴 수 있습니다.

ID conversion은 각 gene의 매핑이 독립적이므로 목록을 제한에 맞게 나눠 요청한 뒤 다시 합칠 수 있습니다. `auto_david_gene_id_conversion_batch.py`가 이 과정을 자동화합니다.

GO/KEGG enrichment는 목록을 나누면 gene 수와 배경이 달라져 p-value와 fold enrichment가 바뀝니다. 따라서 `auto_david_webservice_enrichment.py`는 URL 제한이 없는 DAVID SOAP Web Service로 전체 목록을 한 번에 업로드하고 GO BP FAT와 KEGG를 계산합니다.

`auto_david_merge_enrichment_tables.py`는 여러 결과의 `Category + Term` 합집합을 기준으로 필요한 열만 가로 병합합니다.

## 파일 구성

- `auto_david_gene_id_conversion_batch.py`: Ensembl → DAVID ID/Entrez ID/gene 정보 변환
- `auto_david_webservice_enrichment.py`: 전체 gene list의 GO BP FAT 및 KEGG enrichment
- `auto_david_merge_enrichment_tables.py`: 여러 enrichment TSV/CSV 가로 병합
- `auto_david_command_examples.sh`: 입력 변수를 설정하고 Python 명령을 작성하는 Bash 예시
- `example.txt`: 테스트용 Ensembl ID 목록

`functions/`는 위 Python 실행 파일들이 공통으로 사용하는 내부 함수 모음입니다. 직접 실행하지 않으며, 패키지를 옮길 때 함께 유지해야 합니다.

- `functions/gene_id_conversion_functions.py`: 입력 처리, DAVID linking API, species 필터, MyGene.info 조회
- `functions/enrichment_webservice_functions.py`: DAVID SOAP 인증, 목록 업로드, chart 조회, 결과 저장
- `functions/__init__.py`: 내부 Python package 정의

## 요구사항

- Python 3.10 이상 권장
- 인터넷 연결
- ID conversion: 로그인 불필요
- Web Service enrichment: [DAVID Web Service에 등록한 이메일](https://davidbioinformatics.nih.gov/webservice/register.htm) 필요
- 추가 Python package 불필요

## Species 입력

Gene ID conversion과 GO/KEGG enrichment는 `--species` 입력이 필수입니다. 기본 species를 임의로 적용하지 않으므로 분석 대상의 학명을 명시해야 합니다.

대표적인 모델 생물 species 예시:

- 사람: `Homo sapiens`
- 생쥐: `Mus musculus`
- 랫드: `Rattus norvegicus`
- 초파리: `Drosophila melanogaster`
- 애기장대: `Arabidopsis thaliana`

종명은 따옴표로 묶어 입력하는 것을 권장합니다. 결과표 병합은 이미 생성된 파일을 단순히 조합하며 DAVID에 gene을 제출하지 않으므로 species를 입력하지 않습니다.

스크립트는 자신의 위치를 기준으로 `functions/`를 찾습니다. 패키지는 원하는 위치에 둘 수 있고, 입력·출력 상대 경로는 명령을 실행한 현재 작업 디렉터리를 기준으로 해석됩니다.

아래 예시에서는 패키지 위치를 변수로 지정합니다.

```bash
AUTO_DAVID_DIR="/path/to/260825_auto_david_general"
```

모든 Python 실행 파일은 순서에 의존하는 위치 인자를 사용하지 않습니다. `--input`과 `--output`을 명시하므로 다른 옵션과 순서를 바꿔도 됩니다.

## Bash 명령 예시

`auto_david_command_examples.sh`에는 세 작업의 입력·출력 및 주요 옵션 변수가 별도 영역으로 정리되어 있습니다. 필요한 값을 수정하고 해당 Python 명령의 주석을 제거해 사용하면 됩니다. 별도의 중간 실행 계층이나 subcommand 단계는 없습니다.

## 1. Gene ID conversion batch

Ensembl gene ID 목록을 DAVID URL 길이와 요청당 gene 수 제한에 맞춰 자동 분할합니다. DAVID 규정에 맞춰 API 요청 사이에 최소 10초 대기하고, 결과를 입력 순서대로 하나의 TSV로 합칩니다.

```bash
python "${AUTO_DAVID_DIR}/auto_david_gene_id_conversion_batch.py" \
  --input genes.txt \
  --output gene_conversion.tsv \
  --species "Homo sapiens"
```

옵션 순서는 자유롭게 바꿀 수 있습니다.

```bash
python "${AUTO_DAVID_DIR}/auto_david_gene_id_conversion_batch.py" \
  --species "Homo sapiens" \
  --output gene_conversion.tsv \
  --input genes.txt
```

주요 출력 열:

```text
Original_ID
Query_ID
DAVID_ID
ENTREZ_GENE_ID
Gene_Symbol
Gene_Name
Species
Status
```

`DAVID_ID`는 DAVID 내부 gene cluster ID이며 Entrez ID가 아닙니다. `ENTREZ_GENE_ID`는 동일한 Ensembl ID를 MyGene.info batch API로 조회한 결과입니다.

`--species`는 필수이며 지정한 species의 매핑만 결과에 사용합니다. `ENSG00000123456.11`처럼 Ensembl ID 끝에 붙은 `.숫자` version suffix는 기본적으로 제거한 뒤 DAVID와 MyGene.info에 전달합니다. Conversion 결과에서는 원래 ID를 `Original_ID`, version을 제거한 ID를 `Query_ID`에 기록합니다. Version을 제거하지 않으려면 `--keep-version`을 사용합니다.

## 2. DAVID Web Service enrichment

먼저 [DAVID Web Service 등록 페이지](https://davidbioinformatics.nih.gov/webservice/register.htm)에서 이메일을 등록해야 합니다.

```bash
python "${AUTO_DAVID_DIR}/auto_david_webservice_enrichment.py" \
  --input genes.txt \
  --output result \
  --email "your_registered_email@example.org" \
  --id-type auto \
  --species "Homo sapiens" \
  --ease 1.0 \
  --count 1
```

`--output`은 파일명이 아니라 출력 prefix입니다. 위 명령은 다음 두 파일을 만듭니다.

```text
result.go_bp_fat.tsv
result.kegg_pathway.tsv
```

등록 이메일은 환경변수로 전달할 수도 있습니다.

```bash
export DAVID_EMAIL="your_registered_email@example.org"
python "${AUTO_DAVID_DIR}/auto_david_webservice_enrichment.py" \
  --input genes.txt --output result --species "Homo sapiens"
```

`--id-type`은 `auto`, `ensembl`, `entrez`, `symbol`을 지원합니다. 기본값 `auto`는 전체 목록이 `ENS...G숫자`이면 Ensembl, 모두 숫자이면 Entrez, 그 외에는 gene symbol로 판별합니다. 한 입력 파일에는 한 종류의 ID만 넣으세요.

Enrichment에서도 `ENSG00000123456.11`처럼 끝에 version이 붙은 Ensembl ID를 자동으로 인식하고 `.11` 부분을 제거한 뒤 전체 목록을 DAVID Web Service에 업로드합니다. 이 처리는 `--id-type auto`와 `--id-type ensembl` 모두에 적용되며, version을 유지하려는 경우에만 `--keep-version`을 사용합니다.

DAVID Web Service는 gene symbol 입력 유형을 직접 제공하지 않으므로 symbol 목록은 species 기준 Entrez ID로 먼저 변환합니다. Ensembl과 Entrez 목록은 분할하지 않고 SOAP POST body로 한 번에 업로드합니다.

기본 chart 조건은 최대 EASE score/P-value `1.0`, 최소 gene count `1`입니다. 따라서 EASE와 count 기준으로 결과를 추가로 거르지 않고 DAVID가 반환할 수 있는 항목을 모두 받습니다. 필요한 경우 `--ease 0.05` 또는 `--count 2`처럼 명시적으로 제한할 수 있습니다. SOAP 응답의 `afdr`는 백분율이므로 출력 `FDR` 열은 0–1 범위로 변환합니다.

## 3. 여러 enrichment 결과 병합

이 기능은 연구실에서 사용하던 원본 MATLAB 코드 `DAVIDtable_20170418.m`의 DAVID table 병합 로직을 참고했습니다. 원본 코드는 하나의 template 파일에 있는 `Term`을 기준으로 다른 파일의 `Count`, `PValue`, `Genes`를 VLOOKUP처럼 붙였지만, AutoDavid에서는 별도 template 없이 입력한 모든 파일의 `Category + Term` 합집합을 사용하도록 변경했습니다. 특정 파일에 해당 항목이 없으면 원본 로직과 동일하게 `Count=0`, `PValue=1`, `Genes=#N/A`를 입력합니다.

`--input` 뒤에 파일을 원하는 순서대로 나열합니다. TSV와 CSV를 섞을 수 있으며 2개, 3개, 20개 등 파일 수 제한 없이 사용할 수 있습니다.

```bash
python "${AUTO_DAVID_DIR}/auto_david_merge_enrichment_tables.py" \
  --input group1.tsv group2.csv group3.tsv \
  --output merged.tsv
```

파일은 입력 순서대로 `Count1/PValue1/Genes1`, `Count2/PValue2/Genes2` 형태로 쌓입니다.

```text
            group1.tsv                          group2.csv
Category    Term    Count1 PValue1 Genes1        Count2 PValue2 Genes2
```

- 병합 key: `Category + Term`
- 행 순서: 첫 파일 순서부터 시작하여 다음 파일에서 처음 발견된 term을 차례로 추가
- 별도 base/union 파일: 사용하지 않음
- 해당 term이 특정 파일에 없음: `Count=0`, `PValue=1`, `Genes=#N/A`
- term은 있지만 개별 값이 비어 있음: 해당 열에도 위 기본값 적용
- 입력 파일 내 중복 `Category + Term`: 오류로 중단
- 필수 열: `Category`, `Term`, `Count`, `PValue`, `Genes`

## 상세 옵션 확인

```bash
python "${AUTO_DAVID_DIR}/auto_david_gene_id_conversion_batch.py" --help
python "${AUTO_DAVID_DIR}/auto_david_webservice_enrichment.py" --help
python "${AUTO_DAVID_DIR}/auto_david_merge_enrichment_tables.py" --help
```

## 권장 전체 흐름

```text
gene ID 목록
    ├─ ID 확인/변환 → auto_david_gene_id_conversion_batch.py
    └─ GO/KEGG 분석 → auto_david_webservice_enrichment.py
                              ↓
                       조건별 GO/KEGG TSV
                              ↓
                 auto_david_merge_enrichment_tables.py
                              ↓
                       비교용 merged TSV
```

## DAVID 사용 제한과 주의사항

- Linking API: 전체 URL 약 2,048자 제한
- Linking API: light-duty gene list 약 400 genes 제한
- Linking API: 요청 사이 최소 10초, 사용자/컴퓨터당 하루 200 hits
- Web Service: 일반적으로 최대 3,000 genes 범위 권장
- Web Service: 사용자/컴퓨터당 하루 200 jobs
- DAVID는 부적절한 대량 반복 호출을 예고 없이 제한할 수 있음
- 3,000개를 넘는 목록은 DAVID 팀 문의 또는 DAVID Knowledgebase 기반 로컬 분석 권장

공식 문서:

- [DAVID linking API](https://davidbioinformatics.nih.gov/content.jsp?file=DAVID_API.html)
- [DAVID Web Service](https://davidbioinformatics.nih.gov/content.jsp?file=WS.html)
- [DAVID Web Service operations](https://davidbioinformatics.nih.gov/webservice/services/listServices)
