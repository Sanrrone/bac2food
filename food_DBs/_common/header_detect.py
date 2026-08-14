"""Programmatic header-row detection for messy food-DB Excel/CSV inputs.

Replaces the fragile `skiprows=N` / `iloc[2:]` pattern that drops rows
asymmetrically across multi-sheet workbooks. Searches the first N rows
for any of a known set of food-id column names and returns the matched
row index.
"""
from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd


# Known column-name needles ordered by specificity. First match wins.
DEFAULT_NEEDLES: tuple[str, ...] = (
    "Food Code", "Food code", "FoodCode",
    "Food Item ID", "Food ID", "FoodID", "Food id",
    "Food Number", "FoodNumber",
    "item_no", "Item No", "Item No.",
    "FDC_ID", "fdc_id",
    "Food Name", "Food name",
    "Code", "ID",
)


def find_header_row(
    excel_path: str,
    sheet: str | int = 0,
    needles: Iterable[str] = DEFAULT_NEEDLES,
    max_scan: int = 30,
) -> Optional[int]:
    """Scan the first `max_scan` rows of a sheet for any header needle.

    Returns the 0-indexed row to pass as `header=...`. None if not found.
    """
    try:
        raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=max_scan,
                            engine="openpyxl")
    except Exception:
        try:
            raw = pd.read_excel(excel_path, sheet_name=sheet, header=None, nrows=max_scan)
        except Exception:
            return None
    needle_set = {n.strip().lower() for n in needles}
    for i in range(min(max_scan, len(raw))):
        row_cells = [str(c).strip().lower() for c in raw.iloc[i].tolist() if pd.notna(c)]
        if any(cell in needle_set for cell in row_cells):
            return i
    return None


def find_header_row_csv(
    csv_path: str,
    needles: Iterable[str] = DEFAULT_NEEDLES,
    max_scan: int = 30,
    encoding: str = "utf-8",
    delimiter: str = ",",
) -> Optional[int]:
    """Same idea for CSV/TSV files."""
    try:
        raw = pd.read_csv(csv_path, header=None, nrows=max_scan,
                          encoding=encoding, sep=delimiter, dtype=str,
                          engine="python", on_bad_lines="skip")
    except Exception:
        return None
    needle_set = {n.strip().lower() for n in needles}
    for i in range(min(max_scan, len(raw))):
        row_cells = [str(c).strip().lower() for c in raw.iloc[i].tolist() if pd.notna(c)]
        if any(cell in needle_set for cell in row_cells):
            return i
    return None


def read_merged_header(
    excel_path: str,
    sheet: str | int,
    header_row: int,
    n_header_rows: int = 2,
) -> list[str]:
    """Read N header rows and concatenate them column-wise with ' / '.

    Used for BioFoodComp / AFCD where the upper row gives the nutrient family
    (e.g. 'Fatty acids') and the lower row gives the specific measurement
    ('C18:2 n-6 cis'). Returns column names = ['Fatty acids / C18:2 n-6 cis', ...].
    Forward-fills the upper row across merged cells.
    """
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None,
                        skiprows=header_row, nrows=n_header_rows, engine="openpyxl")
    if len(raw) < n_header_rows:
        return [str(c) for c in raw.iloc[0].tolist()]
    # Forward-fill the upper row (merged cells in xlsx appear as NaN in subsequent columns)
    upper = raw.iloc[0].ffill()
    cols: list[str] = []
    for j in range(raw.shape[1]):
        parts = []
        for i in range(n_header_rows):
            v = raw.iloc[i, j] if i > 0 else upper.iloc[j]
            if pd.notna(v) and str(v).strip() and str(v).lower() != "nan":
                parts.append(str(v).strip())
        # Dedup if upper == lower
        if len(parts) == 2 and parts[0] == parts[1]:
            parts = parts[:1]
        cols.append(" / ".join(parts) if parts else f"_col{j}")
    return cols
