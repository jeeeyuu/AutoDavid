#!/usr/bin/env python3
"""Shared client utilities for DAVID's public linking API."""

from __future__ import annotations

import csv
import fcntl
import html
import json
import re
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from html.parser import HTMLParser
from http.cookiejar import CookieJar
from pathlib import Path
from typing import Iterable, Iterator


BASE_URL = "https://davidbioinformatics.nih.gov/"
API_URL = urllib.parse.urljoin(BASE_URL, "api.jsp")
MAX_GENES = 400
MAX_URL_LENGTH = 2048
MIN_REQUEST_INTERVAL = 10.0
MYGENE_GENE_URL = "https://mygene.info/v3/gene"
MYGENE_QUERY_URL = "https://mygene.info/v3/query"
MYGENE_MAX_IDS = 1000
DAVID_ID_TYPES = {
    "ensembl": "ENSEMBL_GENE_ID",
    "entrez": "ENTREZ_GENE_ID",
    "symbol": "GENE_SYMBOL",
}


class DavidApiError(RuntimeError):
    """Raised when DAVID rejects a request or returns an unexpected page."""


@dataclass(frozen=True)
class InputId:
    original: str
    query: str


def read_ensembl_ids(
    path: Path,
    keep_version: bool = False,
    enforce_total_limit: bool = True,
) -> list[InputId]:
    """Read comma/whitespace-separated IDs, preserving the original values."""
    try:
        text = path.read_text(encoding="utf-8-sig")
    except OSError as exc:
        raise DavidApiError(f"입력 파일을 읽을 수 없습니다: {path}: {exc}") from exc

    tokens = [token for token in re.split(r"[\s,]+", text.strip()) if token]
    if not tokens:
        raise DavidApiError("입력 파일에 Ensembl ID가 없습니다.")

    records: list[InputId] = []
    seen: set[str] = set()
    for token in tokens:
        query = token if keep_version else re.sub(r"\.\d+$", "", token)
        if query in seen:
            continue
        seen.add(query)
        records.append(InputId(original=token, query=query))

    if enforce_total_limit and len(records) > MAX_GENES:
        raise DavidApiError(
            f"단일 DAVID 요청 제한({MAX_GENES} genes)을 초과했습니다: "
            f"중복 제거 후 {len(records)}개"
        )
    return records


def read_analysis_ids(
    path: Path,
    id_type: str = "auto",
    keep_version: bool = False,
    enforce_total_limit: bool = True,
) -> tuple[list[InputId], str]:
    """Read one uniform ID list and return records plus DAVID's ID type name."""
    try:
        text = path.read_text(encoding="utf-8-sig")
    except OSError as exc:
        raise DavidApiError(f"입력 파일을 읽을 수 없습니다: {path}: {exc}") from exc
    tokens = [token for token in re.split(r"[\s,]+", text.strip()) if token]
    if not tokens:
        raise DavidApiError("입력 파일에 gene ID가 없습니다.")

    if id_type == "auto":
        looks_ensembl = [
            bool(re.fullmatch(r"ENS[A-Z0-9]*G\d+(?:\.\d+)?", token, re.I))
            for token in tokens
        ]
        looks_entrez = [bool(re.fullmatch(r"\d+", token)) for token in tokens]
        if all(looks_ensembl):
            id_type = "ensembl"
        elif all(looks_entrez):
            id_type = "entrez"
        elif not any(looks_ensembl) and not any(looks_entrez):
            id_type = "symbol"
        else:
            raise DavidApiError(
                "입력 목록에 서로 다른 ID 유형이 섞여 있어 자동 판별할 수 없습니다. "
                "Ensembl, Entrez 또는 gene symbol 중 한 유형만 넣거나 "
                "--id-type을 명시하세요."
            )
    if id_type not in DAVID_ID_TYPES:
        raise DavidApiError(f"지원하지 않는 ID 유형입니다: {id_type}")

    records: list[InputId] = []
    seen: set[str] = set()
    for token in tokens:
        query = token
        if id_type == "ensembl" and not keep_version:
            query = re.sub(r"\.\d+$", "", token)
        if query in seen:
            continue
        seen.add(query)
        records.append(InputId(original=token, query=query))
    if enforce_total_limit and len(records) > MAX_GENES:
        raise DavidApiError(
            f"단일 DAVID 요청 제한({MAX_GENES} genes)을 초과했습니다: "
            f"중복 제거 후 {len(records)}개"
        )
    return records, DAVID_ID_TYPES[id_type]


def unique_query_ids(records: Iterable[InputId]) -> list[str]:
    return [record.query for record in records]


def build_api_url(
    ids: list[str],
    tool: str,
    annotation: str = "",
    david_id_type: str = "ENSEMBL_GENE_ID",
) -> str:
    params = {
        "type": david_id_type,
        "ids": ",".join(ids),
        "tool": tool,
    }
    if annotation:
        params["annot"] = annotation
    return f"{API_URL}?{urllib.parse.urlencode(params)}"


def split_ids_for_api(
    ids: list[str], tool: str, annotation: str = ""
) -> list[list[str]]:
    """Greedily create the largest chunks whose encoded URLs fit DAVID's limit."""
    chunks: list[list[str]] = []
    current: list[str] = []
    for gene_id in ids:
        candidate = [*current, gene_id]
        if (
            len(candidate) <= MAX_GENES
            and len(build_api_url(candidate, tool, annotation)) <= MAX_URL_LENGTH
        ):
            current = candidate
            continue
        if not current:
            raise DavidApiError(
                f"ID 하나만으로 DAVID URL 제한({MAX_URL_LENGTH}자)을 초과합니다: "
                f"{gene_id}"
            )
        chunks.append(current)
        current = [gene_id]
        if len(build_api_url(current, tool, annotation)) > MAX_URL_LENGTH:
            raise DavidApiError(
                f"ID 하나만으로 DAVID URL 제한({MAX_URL_LENGTH}자)을 초과합니다: "
                f"{gene_id}"
            )
    if current:
        chunks.append(current)
    return chunks


@contextmanager
def _global_rate_limit() -> Iterator[None]:
    """Serialize API hits and maintain DAVID's required 10-second interval."""
    lock_path = Path(tempfile.gettempdir()) / "david_api_request.lock"
    stamp_path = Path(tempfile.gettempdir()) / "david_api_last_request.txt"
    with lock_path.open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            try:
                last_request = float(stamp_path.read_text(encoding="ascii"))
            except (OSError, ValueError):
                last_request = 0.0
            delay = MIN_REQUEST_INTERVAL - (time.time() - last_request)
            if delay > 0:
                print(
                    f"DAVID 호출 간격 준수를 위해 {delay:.1f}초 대기합니다...",
                    file=sys.stderr,
                )
                time.sleep(delay)
            yield
        finally:
            # A failed HTTP attempt is still an API hit for rate-limit purposes.
            stamp_path.write_text(str(time.time()), encoding="ascii")
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


class DavidApiClient:
    def __init__(self, timeout: float = 120.0) -> None:
        self.timeout = timeout
        cookie_jar = CookieJar()
        self.opener = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(cookie_jar)
        )
        self.opener.addheaders = [("User-Agent", "david-api-scripts/1.0")]

    def _open(self, request: str | urllib.request.Request) -> bytes:
        try:
            with self.opener.open(request, timeout=self.timeout) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            raise DavidApiError(f"DAVID HTTP 오류: {exc.code} {exc.reason}") from exc
        except urllib.error.URLError as exc:
            raise DavidApiError(f"DAVID 연결 오류: {exc.reason}") from exc

    def run(
        self,
        ids: list[str],
        tool: str,
        annotation: str = "",
        allow_empty: bool = False,
        david_id_type: str = "ENSEMBL_GENE_ID",
    ) -> str:
        if david_id_type not in DAVID_ID_TYPES.values():
            raise DavidApiError(f"지원하지 않는 DAVID ID 유형입니다: {david_id_type}")
        api_url = build_api_url(ids, tool, annotation, david_id_type)
        if len(api_url) > MAX_URL_LENGTH:
            raise DavidApiError(
                f"요청 URL이 DAVID 제한({MAX_URL_LENGTH}자)을 초과합니다: "
                f"{len(api_url)}자. ID 수를 줄여 다시 실행하세요."
            )

        with _global_rate_limit():
            landing_html = self._open(api_url).decode("utf-8", errors="replace")

        rowids = _extract_js_value(landing_html, "rowids")
        annot = _extract_js_value(landing_html, "annot", required=False)
        action = _extract_form_action(landing_html)
        expected_action = {
            "geneReport": "geneReport.jsp",
            "chartReport": "chartReport.jsp",
        }.get(tool)
        if expected_action is None or action != expected_action:
            raise DavidApiError(
                f"DAVID가 예상하지 못한 결과 페이지를 반환했습니다: {action or '없음'}"
            )
        if not rowids and allow_empty:
            return ""
        if not rowids:
            raise DavidApiError("DAVID에 매핑된 ID가 없습니다.")

        post_data = urllib.parse.urlencode(
            {"rowids": rowids, "annot": annot}
        ).encode("ascii")
        request = urllib.request.Request(
            urllib.parse.urljoin(BASE_URL, action), data=post_data, method="POST"
        )
        return self._open(request).decode("utf-8", errors="replace")

    def download_report(self, report_html: str) -> bytes:
        match = re.search(
            r'href=["\']([^"\']*data/download/[^"\']+)["\']', report_html
        )
        if not match:
            plain = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", report_html))
            detail = html.unescape(plain).strip()[:300]
            raise DavidApiError(
                "DAVID 결과에서 다운로드 파일을 찾지 못했습니다. "
                f"응답 요약: {detail}"
            )
        return self._open(urllib.parse.urljoin(BASE_URL, html.unescape(match.group(1))))


def _extract_js_value(page: str, name: str, required: bool = True) -> str:
    match = re.search(rf"\b{name}\.value\s*=\s*[\"']([^\"']*)[\"']", page)
    if match:
        return html.unescape(match.group(1))
    if required:
        raise DavidApiError(f"DAVID 응답에서 {name!r} 값을 찾지 못했습니다.")
    return ""


def _extract_form_action(page: str) -> str:
    match = re.search(r"apiForm\.action\s*=\s*[\"']([^\"']+)[\"']", page)
    return html.unescape(match.group(1)) if match else ""


@dataclass(frozen=True)
class GeneMapping:
    source_id: str
    david_id: str
    gene_name: str
    gene_symbol: str
    species: str


class _GeneReportParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.in_data_table = False
        self.table_depth = 0
        self.in_cell = False
        self.cell_text: list[str] = []
        self.cells: list[str] = []
        self.row_david_id = ""
        self.rows: list[GeneMapping] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag == "table":
            if self.in_data_table:
                self.table_depth += 1
            elif "dataTable" in (attributes.get("class") or "").split():
                self.in_data_table = True
                self.table_depth = 1
        if not self.in_data_table:
            return
        if tag == "tr":
            self.cells = []
            self.row_david_id = ""
        elif tag in {"td", "th"}:
            self.in_cell = True
            self.cell_text = []
        elif tag == "a" and self.in_cell:
            href = attributes.get("href") or ""
            match = re.search(r"geneReportFull\.jsp\?rowids=([^&]+)", href)
            if match:
                self.row_david_id = urllib.parse.unquote(match.group(1))

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.cell_text.append(data)

    def handle_endtag(self, tag: str) -> None:
        if not self.in_data_table:
            return
        if tag in {"td", "th"} and self.in_cell:
            self.cells.append(" ".join("".join(self.cell_text).split()))
            self.in_cell = False
        elif tag == "tr" and len(self.cells) >= 4 and self.row_david_id:
            full_name = self.cells[1]
            symbol_match = re.search(r"\(([^()]*)\)\s*$", full_name)
            symbol = symbol_match.group(1) if symbol_match else ""
            gene_name = full_name[: symbol_match.start()].strip() if symbol_match else full_name
            self.rows.append(
                GeneMapping(
                    source_id=self.cells[0],
                    david_id=self.row_david_id,
                    gene_name=gene_name,
                    gene_symbol=symbol,
                    species=self.cells[3],
                )
            )
        elif tag == "table":
            self.table_depth -= 1
            if self.table_depth == 0:
                self.in_data_table = False


def parse_gene_report(page: str) -> list[GeneMapping]:
    parser = _GeneReportParser()
    parser.feed(page)
    return parser.rows


def fetch_entrez_mappings(
    ids: list[str], timeout: float = 120.0
) -> dict[str, list[str]]:
    """Map Ensembl gene IDs to explicit Entrez IDs with MyGene.info batch API."""
    mappings: dict[str, list[str]] = {}
    for start in range(0, len(ids), MYGENE_MAX_IDS):
        chunk = ids[start : start + MYGENE_MAX_IDS]
        request = urllib.request.Request(
            MYGENE_GENE_URL,
            data=urllib.parse.urlencode(
                {"ids": ",".join(chunk), "fields": "entrezgene,taxid"}
            ).encode("ascii"),
            headers={
                "Accept": "application/json",
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": "david-api-scripts/1.0",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            raise DavidApiError(
                f"MyGene.info Entrez 변환 HTTP 오류: {exc.code} {exc.reason}"
            ) from exc
        except urllib.error.URLError as exc:
            raise DavidApiError(
                f"MyGene.info Entrez 변환 연결 오류: {exc.reason}"
            ) from exc
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise DavidApiError(
                "MyGene.info Entrez 변환 응답을 해석하지 못했습니다."
            ) from exc

        if not isinstance(payload, list):
            raise DavidApiError("MyGene.info가 예상하지 못한 응답을 반환했습니다.")
        for item in payload:
            if not isinstance(item, dict) or item.get("notfound"):
                continue
            query = str(item.get("query", ""))
            entrez = item.get("entrezgene")
            if query and entrez is not None:
                value = str(entrez)
                values = mappings.setdefault(query, [])
                if value not in values:
                    values.append(value)
    return mappings


def fetch_symbol_entrez_mappings(
    symbols: list[str], species: str, timeout: float = 120.0
) -> dict[str, list[str]]:
    """Resolve gene symbols to Entrez IDs for one species with MyGene.info."""
    mappings: dict[str, list[str]] = {}
    for start in range(0, len(symbols), MYGENE_MAX_IDS):
        chunk = symbols[start : start + MYGENE_MAX_IDS]
        request = urllib.request.Request(
            MYGENE_QUERY_URL,
            data=urllib.parse.urlencode(
                {
                    "q": ",".join(chunk),
                    "scopes": "symbol",
                    "fields": "entrezgene,taxid",
                    "species": species,
                }
            ).encode("utf-8"),
            headers={
                "Accept": "application/json",
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": "david-api-scripts/1.0",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            raise DavidApiError(
                f"MyGene.info symbol 변환 HTTP 오류: {exc.code} {exc.reason}"
            ) from exc
        except urllib.error.URLError as exc:
            raise DavidApiError(
                f"MyGene.info symbol 변환 연결 오류: {exc.reason}"
            ) from exc
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise DavidApiError(
                "MyGene.info symbol 변환 응답을 해석하지 못했습니다."
            ) from exc

        if not isinstance(payload, list):
            raise DavidApiError("MyGene.info가 예상하지 못한 응답을 반환했습니다.")
        for item in payload:
            if not isinstance(item, dict) or item.get("notfound"):
                continue
            query = str(item.get("query", ""))
            entrez = item.get("entrezgene")
            if query and entrez is not None:
                value = str(entrez)
                values = mappings.setdefault(query, [])
                if value not in values:
                    values.append(value)
    return mappings


def prepare_analysis_ids(
    records: list[InputId],
    david_id_type: str,
    species: str,
    timeout: float = 120.0,
    max_genes: int | None = MAX_GENES,
) -> tuple[list[str], str, int]:
    """Convert symbol input to Entrez when needed and return DAVID-ready IDs."""
    query_ids = unique_query_ids(records)
    if david_id_type != "GENE_SYMBOL":
        return query_ids, david_id_type, len(query_ids)

    mappings = fetch_symbol_entrez_mappings(query_ids, species, timeout)
    resolved: list[str] = []
    seen: set[str] = set()
    for symbol in query_ids:
        for entrez_id in mappings.get(symbol, []):
            if entrez_id not in seen:
                seen.add(entrez_id)
                resolved.append(entrez_id)
    if not resolved:
        raise DavidApiError(
            f"species={species!r}에서 Entrez ID로 변환된 gene symbol이 없습니다."
        )
    if max_genes is not None and len(resolved) > max_genes:
        raise DavidApiError(
            f"symbol 변환 후 gene 제한({max_genes} genes)을 "
            f"초과했습니다: {len(resolved)}개 Entrez IDs"
        )
    return resolved, "ENTREZ_GENE_ID", len(mappings)


def write_bytes(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.write_bytes(content)
    except OSError as exc:
        raise DavidApiError(f"결과 파일을 쓸 수 없습니다: {path}: {exc}") from exc


def write_conversion_tsv(
    path: Path,
    inputs: list[InputId],
    mappings: list[GeneMapping],
    expected_species: str | None = None,
    entrez_by_source: dict[str, list[str]] | None = None,
) -> None:
    by_source: dict[str, list[GeneMapping]] = {}
    for mapping in mappings:
        by_source.setdefault(mapping.source_id, []).append(mapping)

    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "Original_ID",
                    "Query_ID",
                    "DAVID_ID",
                    "ENTREZ_GENE_ID",
                    "Gene_Symbol",
                    "Gene_Name",
                    "Species",
                    "Status",
                ]
            )
            for record in inputs:
                entrez_ids = (entrez_by_source or {}).get(record.query, [])
                entrez_text = ",".join(entrez_ids)
                all_hits = by_source.get(record.query, [])
                hits = all_hits
                if expected_species is not None:
                    hits = [
                        hit
                        for hit in all_hits
                        if hit.species.casefold() == expected_species.casefold()
                    ]
                if not hits:
                    if all_hits and expected_species is not None:
                        observed = ",".join(
                            sorted({hit.species for hit in all_hits if hit.species})
                        )
                        status = f"species_mismatch:{observed or 'unknown'}"
                    else:
                        status = "unmapped"
                    writer.writerow(
                        [
                            record.original,
                            record.query,
                            "",
                            entrez_text,
                            "",
                            "",
                            expected_species or "",
                            "entrez_only"
                            if status == "unmapped" and entrez_ids
                            else status,
                        ]
                    )
                    continue
                status = "mapped" if len(hits) == 1 else "ambiguous"
                for hit in hits:
                    writer.writerow(
                        [
                            record.original,
                            record.query,
                            hit.david_id,
                            entrez_text,
                            hit.gene_symbol,
                            hit.gene_name,
                            hit.species,
                            status,
                        ]
                    )
    except OSError as exc:
        raise DavidApiError(f"결과 파일을 쓸 수 없습니다: {path}: {exc}") from exc
