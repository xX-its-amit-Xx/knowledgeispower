"""Shared helpers for the Phase B bulk-atlas reference panels (procan, gtex_proteome, cptac). Not a representation.

* ``gene_accessions``: gene id (Ensembl ``ENSG...[.v]`` or symbol) -> canonical UniProt accessions via HGNC
  (``ids.resolve_symbol`` -> ``ids.uniprot_for_symbol``); ambiguous resolutions keep every candidate.
* ``build_index``: accession -> panel row. When several rows claim one accession, the row with the most observed
  samples wins (ties: the earlier row), so each accession points to exactly one row.
* ``read_xlsx_rows``: dependency-free streaming reader for one worksheet of an .xlsx (zip + iterparse), used to
  parse large supplementary tables without openpyxl and without loading the sheet into memory.
* ``corr_graph`` / ``block``: the same construction as ``reps.ccle`` (pairwise-complete Pearson kNN graph; SVD
  loadings of the standardised panel), parameterised by panel.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import networkx as nx
import numpy as np

from ..features import _finish
from ..ids import resolve_symbol, uniprot_for_symbol
from . import _refpanel as rp

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def gene_accessions(gene: str) -> list[str]:
    syms, _ = resolve_symbol(str(gene).split(".")[0] if str(gene).upper().startswith("ENSG") else str(gene))
    out: dict[str, None] = {}
    for s in syms:
        for a in uniprot_for_symbol(s):
            out.setdefault(a, None)
    return list(out)


def build_index(acc_lists: list[list[str]], nobs: np.ndarray) -> dict[str, int]:
    best: dict[str, int] = {}
    for i, accs in enumerate(acc_lists):
        for a in accs:
            j = best.get(a)
            if j is None or nobs[i] > nobs[j]:
                best[a] = i
    return best


def _col(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_xlsx_rows(path: Path, sheet_name: str):
    """Yield each row of worksheet `sheet_name` as {column index: str value} (shared strings resolved)."""
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rid = next(s.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
                   for s in wb.iter(_NS + "sheet") if s.get("name") == sheet_name)
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in rels if r.get("Id") == rid)
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        ss: list[str] = []
        with z.open("xl/sharedStrings.xml") as fh:
            for _, el in ET.iterparse(fh):
                if el.tag == _NS + "si":
                    ss.append("".join(t.text or "" for t in el.iter(_NS + "t")))
                    el.clear()
        with z.open(target) as fh:
            for _, el in ET.iterparse(fh):
                if el.tag != _NS + "row":
                    continue
                row: dict[int, str] = {}
                for c in el.iter(_NS + "c"):
                    v = c.find(_NS + "v")
                    if v is None or v.text is None:
                        continue
                    row[_col(c.get("r"))] = ss[int(v.text)] if c.get("t") == "s" else v.text
                yield row
                el.clear()


def corr_graph(var_names, X: np.ndarray, index: dict[str, int], top: int, min_shared: int, prov: dict,
               unit: str) -> nx.Graph:
    rows = rp.align(var_names, index)
    idx = np.flatnonzero(rows >= 0)
    r, n = rp.pairwise_corr(X[rows[idx]])
    g = rp.knn_graph(var_names, rows, r, top, valid=n >= min_shared)
    g.graph["provenance"] = {**prov, "construction": f"top-{top} pairwise-complete Pearson, >= {min_shared} {unit}"}
    return g


def block(var_names, index: dict[str, int], L: np.ndarray, ok: np.ndarray) -> np.ndarray:
    rows = rp.align(var_names, index)
    has = np.zeros(len(rows), bool)
    has[rows >= 0] = ok[rows[rows >= 0]]
    M = np.zeros((len(rows), L.shape[1]))
    M[has] = L[rows[has]]
    return _finish(M, has)
