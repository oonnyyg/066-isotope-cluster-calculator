"""Element isotope table loading and validation."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .errors import IsotopeError

ABUNDANCE_TOL = 1e-9


@dataclass(frozen=True)
class Isotope:
    mass_number: int
    exact_mass: float
    abundance: float


def validate_table(raw: object) -> dict[str, list[Isotope]]:
    if not isinstance(raw, dict) or not raw:
        raise IsotopeError("元素同位素表必须是非空对象")
    table: dict[str, list[Isotope]] = {}
    for symbol, entries in raw.items():
        where = f"元素 {symbol!r}"
        if not isinstance(symbol, str) or not symbol.isalpha() or not symbol[0].isupper():
            raise IsotopeError(f"{where}: 元素符号非法")
        if not isinstance(entries, list) or not entries:
            raise IsotopeError(f"{where}: 同位素列表必须非空")
        isotopes: list[Isotope] = []
        seen: set[int] = set()
        total = 0.0
        for idx, entry in enumerate(entries):
            loc = f"{where} 第 {idx} 个同位素"
            if not isinstance(entry, dict):
                raise IsotopeError(f"{loc}: 条目必须是对象")
            try:
                mass_number = entry["mass_number"]
                exact_mass = entry["exact_mass"]
                abundance = entry["abundance"]
            except KeyError as exc:
                raise IsotopeError(f"{loc}: 缺少字段 {exc.args[0]}") from None
            if not isinstance(mass_number, int) or isinstance(mass_number, bool) or mass_number <= 0:
                raise IsotopeError(f"{loc}: 质量数必须为正整数")
            if mass_number in seen:
                raise IsotopeError(f"{loc}: 质量数 {mass_number} 在同元素内重复")
            seen.add(mass_number)
            if not isinstance(exact_mass, (int, float)) or isinstance(exact_mass, bool) \
                    or not math.isfinite(exact_mass) or exact_mass <= 0:
                raise IsotopeError(f"{loc}: 精确质量必须为有限正数")
            if not isinstance(abundance, (int, float)) or isinstance(abundance, bool) \
                    or not math.isfinite(abundance) or abundance < 0:
                raise IsotopeError(f"{loc}: 丰度必须为有限非负数")
            total += abundance
            isotopes.append(Isotope(mass_number, float(exact_mass), float(abundance)))
        if abs(total - 1.0) > ABUNDANCE_TOL:
            raise IsotopeError(
                f"{where}: 丰度和 {total!r} 与 1 的差超过 {ABUNDANCE_TOL}"
            )
        if total <= 0:
            raise IsotopeError(f"{where}: 丰度和必须为正")
        table[symbol] = [
            Isotope(i.mass_number, i.exact_mass, i.abundance / total)
            for i in isotopes
        ]
    return table
