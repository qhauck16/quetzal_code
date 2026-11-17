#!/usr/bin/env python3
"""Write a supplementary table TSV as an Excel workbook.

Sheet 1 holds the table (header row bold and frozen, columns auto-sized);
sheet 2 ("README") holds the title and the column descriptions passed in a
two-column TSV (column<TAB>description). Needs openpyxl.

    python3 write_xlsx.py table.tsv columns.tsv "Supplementary Table 1. ..." out.xlsx
"""
import csv
import sys

from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter


def as_number(v):
    if v in ("", "NA"):
        return v if v == "" else None
    try:
        f = float(v)
    except ValueError:
        return v
    return int(f) if f.is_integer() and "." not in v and "e" not in v.lower() else f


def main(table_tsv, columns_tsv, title, out_xlsx):
    wb = Workbook()
    ws = wb.active
    ws.title = "Table"
    with open(table_tsv, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    for i, row in enumerate(rows):
        ws.append(row if i == 0 else [as_number(v) for v in row])
    for c in ws[1]:
        c.font = Font(bold=True)
    ws.freeze_panes = "A2"
    for j, name in enumerate(rows[0], start=1):
        width = max(len(str(r[j - 1])) if j - 1 < len(r) else 0 for r in rows[:200])
        ws.column_dimensions[get_column_letter(j)].width = min(max(10, width + 2), 60)

    rd = wb.create_sheet("README")
    rd.append([title])
    rd["A1"].font = Font(bold=True)
    rd.append([])
    rd.append(["Column", "Description"])
    for c in rd[3]:
        c.font = Font(bold=True)
    with open(columns_tsv, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row:
                rd.append(row)
    rd.column_dimensions["A"].width = 30
    rd.column_dimensions["B"].width = 110
    wb.save(out_xlsx)
    print("wrote", out_xlsx)


if __name__ == "__main__":
    main(*sys.argv[1:5])
