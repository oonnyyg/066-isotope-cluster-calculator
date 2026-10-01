"""Molecular formula parser: symbols, integer subscripts, nested parens."""

from __future__ import annotations

from .errors import IsotopeError


def parse_formula(formula: str) -> dict[str, int]:
    if not isinstance(formula, str) or not formula:
        raise IsotopeError("分子式为空或不是字符串")
    counts, pos = _parse_group(formula, 0, top=True)
    if pos != len(formula):
        raise IsotopeError(
            f"分子式 {formula!r} 在位置 {pos} 处存在多余字符 {formula[pos]!r}"
        )
    if not counts:
        raise IsotopeError(f"分子式 {formula!r} 为空")
    return counts


def _parse_group(text: str, pos: int, top: bool) -> tuple[dict[str, int], int]:
    counts: dict[str, int] = {}
    n = len(text)
    while pos < n:
        ch = text[pos]
        if ch == "(":
            sub, pos = _parse_group(text, pos + 1, top=False)
            mult, pos = _parse_number(text, pos)
            for sym, num in sub.items():
                counts[sym] = counts.get(sym, 0) + num * mult
        elif ch == ")":
            if top:
                raise IsotopeError(f"分子式 {text!r} 位置 {pos} 处括号未配对")
            if not counts:
                raise IsotopeError(f"分子式 {text!r} 位置 {pos} 处括号内容为空")
            return counts, pos + 1
        elif ch.isupper():
            start = pos
            pos += 1
            while pos < n and text[pos].islower():
                pos += 1
            symbol = text[start:pos]
            mult, pos = _parse_number(text, pos)
            counts[symbol] = counts.get(symbol, 0) + mult
        else:
            raise IsotopeError(
                f"分子式 {text!r} 位置 {pos} 处存在非法字符 {ch!r}"
            )
    if not top:
        raise IsotopeError(f"分子式 {text!r} 存在未闭合的括号")
    return counts, pos


def _parse_number(text: str, pos: int) -> tuple[int, int]:
    start = pos
    n = len(text)
    while pos < n and text[pos].isdigit():
        pos += 1
    if start == pos:
        return 1, pos
    value = int(text[start:pos])
    if value <= 0:
        raise IsotopeError(f"分子式 {text!r} 位置 {start} 处下标不得为零")
    return value, pos
