#!/usr/bin/env python3
from __future__ import annotations

from datetime import datetime
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "analysis_runs"
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
DEST = RUNS / STAMP

MAP = {
    "figures": [ROOT / "outputs" / "figures", ROOT / "results" / "figures"],
    "tables": [ROOT / "outputs" / "tables", ROOT / "results" / "tables"],
    "logs": [ROOT / "outputs" / "logs", ROOT / "results" / "logs"],
    "objects": [ROOT / "outputs" / "objects"],
    "legends": [ROOT / "outputs" / "legends", ROOT / "results" / "legends"],
}

DEST.mkdir(parents=True, exist_ok=True)
manifest_lines = ["source\tcollected_path\tsize_bytes"]

for group, sources in MAP.items():
    out_dir = DEST / group
    out_dir.mkdir(parents=True, exist_ok=True)
    for src_root in sources:
        if not src_root.exists():
            continue
        for p in src_root.rglob("*"):
            if not p.is_file():
                continue
            rel = p.relative_to(src_root)
            suffix = src_root.name
            target = out_dir / f"{suffix}__{rel.as_posix().replace('/', '__')}"
            if target.exists():
                target = out_dir / f"{suffix}__dup__{rel.as_posix().replace('/', '__')}"
            shutil.copy2(p, target)
            manifest_lines.append(f"{p.relative_to(ROOT)}\t{target.relative_to(ROOT)}\t{target.stat().st_size}")

manifest = DEST / "manifest.tsv"
manifest.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

latest = RUNS / "LATEST"
latest.write_text(str(DEST.relative_to(ROOT)) + "\n", encoding="utf-8")

print(f"Collected outputs into: {DEST}")
print(f"Manifest: {manifest}")
