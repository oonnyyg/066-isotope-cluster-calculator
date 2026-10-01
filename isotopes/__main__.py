"""Command line entry point: ``python -m isotopes input.json --out dir``."""

import argparse
import csv
import json
import re
import sys

from .errors import IsotopeError
from .runner import run


def _safe_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._-")
    return cleaned or "molecule"


def _cluster_to_dict(cluster) -> dict:
    return {
        "formula": cluster.formula,
        "charge": cluster.charge,
        "element_counts": cluster.counts,
        "peak_count": len(cluster.peaks),
        "total_probability": cluster.total_probability,
        "mean_molecular_mass": cluster.mean_molecular_mass,
        "peaks": [
            {
                "mass_number": peak.mass_number,
                "probability": peak.probability,
                "mean_exact_mass": peak.mean_mass,
                "mz": peak.mz,
            }
            for peak in cluster.peaks
        ],
    }


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m isotopes",
        description="计算分子式的整数质量数同位素峰簇分布",
    )
    parser.add_argument("input", help="输入 JSON 文件路径")
    parser.add_argument("--out", required=True, help="输出目录")
    args = parser.parse_args(argv)

    try:
        with open(args.input, "r", encoding="utf-8-sig") as handle:
            data = json.load(handle)
        clusters = run(data)
    except (OSError, json.JSONDecodeError, IsotopeError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1

    import os

    os.makedirs(args.out, exist_ok=True)

    summary = []
    used: dict[str, int] = {}
    for cluster in clusters:
        base = _safe_name(f"{cluster.name}_z{cluster.charge}")
        used[base] = used.get(base, 0) + 1
        filename = base if used[base] == 1 else f"{base}_{used[base]}"
        csv_path = f"{filename}.csv"
        with open(os.path.join(args.out, csv_path), "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["mass_number", "probability", "mean_exact_mass_u", "mz"])
            for peak in cluster.peaks:
                writer.writerow([
                    peak.mass_number,
                    repr(peak.probability),
                    repr(peak.mean_mass),
                    repr(peak.mz),
                ])
        entry = _cluster_to_dict(cluster)
        entry["csv"] = csv_path
        summary.append(entry)

    report = {
        "charge": clusters[0].charge if clusters else data.get("charge"),
        "molecule_count": len(summary),
        "molecules": summary,
    }
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
