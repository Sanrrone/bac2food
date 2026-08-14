"""Detect raw food-DB source file format from magic bytes + extension.

Returns a structured FileInfo describing how to open the file. Used by every
v2 ingester so it never silently misreads (e.g. semicolon-CSV as comma-CSV,
or .xls as .xlsx).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class FileInfo:
    path: Path
    fmt: str          # "xlsx" | "xls" | "csv" | "tsv" | "pdf" | "unknown"
    encoding: Optional[str] = None
    delimiter: Optional[str] = None
    note: str = ""


_XLSX_MAGIC = b"PK\x03\x04"
_XLS_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
_PDF_MAGIC = b"%PDF"


def _sniff_text(path: Path) -> tuple[str, str]:
    """Return (encoding, delimiter) for a text file."""
    raw = path.read_bytes()[: 256 * 1024]
    for enc in ("utf-8", "utf-8-sig", "latin-1", "cp1252"):
        try:
            sample = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    else:
        sample = raw.decode("utf-8", errors="replace")
        enc = "utf-8(repl)"
    try:
        dialect = csv.Sniffer().sniff(sample[:8192], delimiters=",;\t|")
        delim = dialect.delimiter
    except csv.Error:
        # Fallback: count occurrences in first line
        first = sample.split("\n", 1)[0]
        delim = max(",;\t|", key=lambda d: first.count(d))
    return enc, delim


def detect(path: str | Path) -> FileInfo:
    p = Path(path)
    if not p.exists():
        return FileInfo(p, "missing", note="file does not exist")
    head = p.open("rb").read(16)
    ext = p.suffix.lower()

    if head.startswith(_XLSX_MAGIC):
        return FileInfo(p, "xlsx", note=f"OOXML zip magic; ext={ext}")
    if head.startswith(_XLS_MAGIC):
        return FileInfo(p, "xls", note=f"legacy OLE magic; ext={ext}")
    if head.startswith(_PDF_MAGIC):
        return FileInfo(p, "pdf", note="PDF — must be converted to xlsx/csv before ingestion")

    if ext in (".csv", ".tsv", ".txt"):
        enc, delim = _sniff_text(p)
        fmt = "tsv" if delim == "\t" else "csv"
        return FileInfo(p, fmt, encoding=enc, delimiter=delim,
                        note=f"text; ext={ext}; sniffed delim={delim!r}")

    if ext in (".xlsx", ".xlsm"):
        return FileInfo(p, "xlsx", note="ext-only; magic missing")
    if ext == ".xls":
        return FileInfo(p, "xls", note="ext-only; magic missing")

    return FileInfo(p, "unknown", note=f"unrecognized; head={head[:8]!r} ext={ext}")


def describe(info: FileInfo) -> str:
    parts = [f"[format] {info.path.name} → {info.fmt}"]
    if info.encoding:
        parts.append(f"enc={info.encoding}")
    if info.delimiter:
        parts.append(f"delim={info.delimiter!r}")
    if info.note:
        parts.append(info.note)
    return " ".join(parts)
