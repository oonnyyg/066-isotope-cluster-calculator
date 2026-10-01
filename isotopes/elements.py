"""Loading and validation of the offline isotope abundance table."""

import math

from .errors import IsotopeError


class Isotope:
    __slots__ = ("mass_number", "exact_mass", "abundance")

    def __init__(self, mass_number: int, exact_mass: float, abundance: float):
        self.mass_number = mass_number
        self.exact_mass = exact_mass
        self.abundance = abundance


class ElementTable:
    """Maps element symbol -> list of :class:`Isotope`."""

    def __init__(self, raw: object):
        if not isinstance(raw, dict) or not raw:
            raise IsotopeError("元素同位素表必须是非空对象", where="elements")
        self.elements: dict[str, list[Isotope]] = {}
        for symbol, isotopes in raw.items():
            self.elements[symbol] = self._parse_element(symbol, isotopes)

    def symbols(self) -> "set[str]":
        return set(self.elements)

    @staticmethod
    def _parse_element(symbol: object, isotopes: object) -> list[Isotope]:
        if not isinstance(symbol, str) or not symbol:
            raise IsotopeError("元素符号必须是非空字符串", where=f"elements.{symbol!r}")
        if not isinstance(isotopes, list) or not isotopes:
            raise IsotopeError("元素至少需要一个同位素", where=f"元素 {symbol}")
        parsed: list[Isotope] = []
        seen: set[int] = set()
        abundance_sum = 0.0
        for index, item in enumerate(isotopes):
            where = f"元素 {symbol} 第 {index + 1} 个同位素"
            if not isinstance(item, dict):
                raise IsotopeError("同位素条目必须是对象", where=where)
            mass_number = item.get("mass_number")
            exact_mass = item.get("exact_mass")
            abundance = item.get("abundance")
            if not isinstance(mass_number, int) or isinstance(mass_number, bool) or mass_number <= 0:
                raise IsotopeError("质量数必须为正整数", where=where)
            if mass_number in seen:
                raise IsotopeError(f"质量数 {mass_number} 在该元素内重复", where=where)
            seen.add(mass_number)
            if not isinstance(exact_mass, (int, float)) or isinstance(exact_mass, bool):
                raise IsotopeError("精确质量必须是数字", where=where)
            exact_mass = float(exact_mass)
            if not math.isfinite(exact_mass) or exact_mass <= 0.0:
                raise IsotopeError("精确质量必须为有限正数", where=where)
            if not isinstance(abundance, (int, float)) or isinstance(abundance, bool):
                raise IsotopeError("丰度必须是数字", where=where)
            abundance = float(abundance)
            if not math.isfinite(abundance) or abundance < 0.0:
                raise IsotopeError("丰度必须为有限非负数", where=where)
            parsed.append(Isotope(mass_number, exact_mass, abundance))
            abundance_sum += abundance
        if not parsed:  # pragma: no cover - guarded above
            raise IsotopeError("元素至少需要一个同位素", where=f"元素 {symbol}")
        if abs(abundance_sum - 1.0) > 1e-9:
            raise IsotopeError(
                f"丰度和 {abundance_sum!r} 与 1 的差超过 1e-9",
                where=f"元素 {symbol}",
            )
        for isotope in parsed:
            isotope.abundance = isotope.abundance / abundance_sum
        return parsed
