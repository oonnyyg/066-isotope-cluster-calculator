"""CLI: python -m isotopes input.json --out output_dir"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from .core import compute_cluster, summarize
from .errors import IsotopeError
from .formula import parse_formula
from .table import validate_table


def _load_request(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise IsotopeError(f"无法读取输入文件 {path}: {exc}") from None
    except json.JSONDecodeError as exc:
        raise IsotopeError(f"输入文件 {path} 不是合法 JSON: {exc}") from None
    if not isinstance(data, dict):
        raise IsotopeError("输入 JSON 顶层必须是对象")
    return data


def _get_charge(data: dict) -> int:
    charge = data.get("charge")
    if not isinstance(charge, int) or isinstance(charge, bool) or charge <= 0:
        raise IsotopeError(f"电荷数 charge 必须为正整数，得到 {charge!r}")
    return charge


def run(input_path: Path, out_dir: Path) -> dict:
    data = _load_request(input_path)
    if "elements" not in data:
        raise IsotopeError("输入 JSON 缺少 elements 元素同位素表")
    table = validate_table(data["elements"])
    formula = data.get("formula")
    counts = parse_formula(formula)
    charge = _get_charge(data)
    cluster = compute_cluster(counts, table)
    summary = summarize(cluster, charge)
    report = {
        "formula": formula,
        "charge": charge,
        "element_counts": dict(sorted(counts.items())),
        "total_probability": summary["total_probability"],
        "average_mass": summary["average_mass"],
        "average_mz": summary["average_mz"],
        "peak_count": len(summary["peaks"]),
        "peaks": summary["peaks"],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    with (out_dir / "peaks.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["mass_number", "probability", "mean_exact_mass", "mz"])
        for peak in summary["peaks"]:
            writer.writerow(
                [peak["mass_number"], repr(peak["probability"]),
                 repr(peak["mean_exact_mass"]), repr(peak["mz"])]
            )
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m isotopes",
        description="分子式同位素峰簇计算器",
    )
    parser.add_argument("input", type=Path, help="输入 JSON 文件路径")
    parser.add_argument("--out", type=Path, required=True, help="输出目录")
    args = parser.parse_args(argv)
    try:
        report = run(args.input, args.out)
    except IsotopeError as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1
    print(
        f"完成: {report['peak_count']} 个峰，"
        f"总概率 {report['total_probability']:.12g}，"
        f"平均分子质量 {report['average_mass']:.6f} u"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
