#!/usr/bin/env python3
"""Minimal stateful DAVID SOAP Web Service client using only stdlib."""

from __future__ import annotations

import csv
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from http.cookiejar import CookieJar
from pathlib import Path

from .gene_id_conversion_functions import DavidApiError


SOAP_NS = "http://schemas.xmlsoap.org/soap/envelope/"
SERVICE_NS = "http://service.session.sample"
SOAP_ENDPOINT = (
    "https://davidbioinformatics.nih.gov/webservice/services/"
    "DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/"
)
MAX_WEBSERVICE_GENES = 3000

CHART_COLUMNS = [
    ("Category", "categoryName"),
    ("Term", "termName"),
    ("Count", "listHits"),
    ("%", "percent"),
    ("PValue", "ease"),
    ("Genes", "geneIds"),
    ("List Total", "listTotals"),
    ("Pop Hits", "popHits"),
    ("Pop Total", "popTotals"),
    ("Fold Enrichment", "foldEnrichment"),
    ("Bonferroni", "bonferroni"),
    ("Benjamini", "benjamini"),
    ("FDR", "afdr"),
]


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


class DavidWebServiceClient:
    def __init__(self, timeout: float = 600.0) -> None:
        self.timeout = timeout
        self.opener = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(CookieJar())
        )

    def call(self, operation: str, *arguments: object) -> ET.Element:
        envelope = ET.Element(ET.QName(SOAP_NS, "Envelope"))
        body = ET.SubElement(envelope, ET.QName(SOAP_NS, "Body"))
        method = ET.SubElement(body, ET.QName(SERVICE_NS, operation))
        for index, argument in enumerate(arguments):
            element = ET.SubElement(method, ET.QName(SERVICE_NS, f"args{index}"))
            element.text = str(argument)
        payload = ET.tostring(envelope, encoding="utf-8", xml_declaration=True)
        request = urllib.request.Request(
            SOAP_ENDPOINT,
            data=payload,
            headers={
                "Content-Type": "text/xml; charset=utf-8",
                "SOAPAction": f'"urn:{operation}"',
                "User-Agent": "david-webservice-scripts/1.0",
            },
            method="POST",
        )
        try:
            with self.opener.open(request, timeout=self.timeout) as response:
                raw = response.read()
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")[:500]
            raise DavidApiError(
                f"DAVID Web Service HTTP 오류: {exc.code} {exc.reason}: {detail}"
            ) from exc
        except urllib.error.URLError as exc:
            raise DavidApiError(f"DAVID Web Service 연결 오류: {exc.reason}") from exc
        try:
            root = ET.fromstring(raw)
        except ET.ParseError as exc:
            raise DavidApiError("DAVID SOAP 응답 XML을 해석하지 못했습니다.") from exc

        fault = root.find(f".//{{{SOAP_NS}}}Fault")
        if fault is not None:
            message = next(
                (
                    child.text or ""
                    for child in fault.iter()
                    if _local_name(child.tag) in {"faultstring", "Exception"}
                    and child.text
                ),
                "알 수 없는 SOAP fault",
            )
            raise DavidApiError(f"DAVID Web Service 오류: {message}")
        return root

    @staticmethod
    def return_elements(root: ET.Element) -> list[ET.Element]:
        return [element for element in root.iter() if _local_name(element.tag) == "return"]

    def authenticate(self, email: str) -> None:
        values = self.return_elements(self.call("authenticate", email))
        authenticated = values and (values[0].text or "").strip().lower() == "true"
        if not authenticated:
            raise DavidApiError(
                f"DAVID Web Service 인증에 실패했습니다: {email}. 등록 상태를 확인하세요."
            )

    def add_list(self, ids: list[str], id_type: str, list_name: str) -> float:
        values = self.return_elements(
            self.call("addList", ",".join(ids), id_type, list_name, 0)
        )
        try:
            return float(values[0].text or "0")
        except (IndexError, ValueError) as exc:
            raise DavidApiError("DAVID addList 결과를 해석하지 못했습니다.") from exc

    def get_species(self) -> list[str]:
        return [
            (element.text or "").strip()
            for element in self.return_elements(self.call("getSpecies"))
            if (element.text or "").strip()
        ]

    def select_species(self, species_name: str) -> str:
        species = self.get_species()
        wanted = species_name.casefold().replace("_", " ").strip()
        for index, label in enumerate(species):
            name = label.rsplit("(", 1)[0].strip()
            if name.casefold() == wanted:
                self.call("setCurrentSpecies", str(index))
                return label
        available = ", ".join(species) if species else "없음"
        raise DavidApiError(
            f"업로드 목록에서 species {species_name!r}를 찾지 못했습니다. "
            f"DAVID 후보: {available}"
        )

    def set_categories(self, category: str) -> None:
        values = self.return_elements(self.call("setCategories", category))
        accepted = ",".join((value.text or "") for value in values)
        if category not in accepted:
            raise DavidApiError(
                f"DAVID가 annotation category를 승인하지 않았습니다: {category}"
            )

    def get_chart_report(self, ease: float, count: int) -> list[dict[str, str]]:
        returns = self.return_elements(self.call("getChartReport", ease, count))
        records: list[dict[str, str]] = []
        for record in returns:
            values = {
                _local_name(child.tag): (child.text or "")
                for child in list(record)
            }
            if values:
                # SOAP's afdr is a percentage; api.jsp TSV uses a 0-1 fraction.
                if values.get("afdr"):
                    try:
                        values["afdr"] = str(float(values["afdr"]) / 100.0)
                    except ValueError:
                        pass
                records.append(values)
        return records


def write_chart_tsv(path: Path, records: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([header for header, _ in CHART_COLUMNS])
            for record in records:
                writer.writerow([record.get(field, "") for _, field in CHART_COLUMNS])
    except OSError as exc:
        raise DavidApiError(f"결과 파일을 쓸 수 없습니다: {path}: {exc}") from exc
