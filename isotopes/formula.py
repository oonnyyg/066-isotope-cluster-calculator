"""Parser for molecular formulas such as ``Ca(OH)2`` or ``((CH3)2CH)2``."""

from .errors import IsotopeError


def parse_formula(text: str, known_elements: "set[str] | None" = None,
                  *, where: str = "") -> "dict[str, int]":
    """Parse *text* into merged ``{element: count}`` counts.

    Subscripts are positive integers; omitted subscripts mean 1.
    Parenthesised groups may nest. Raises :class:`IsotopeError` with a
    1-based character position for empty formulas/groups, zero or
    malformed subscripts, unbalanced brackets or unknown elements.
    """
    if not isinstance(text, str) or text == "":
        raise IsotopeError("分子式不能为空", where=where, position=1)

    pos = 0
    n = len(text)

    def error(message: str, at: int) -> IsotopeError:
        return IsotopeError(message, where=where, position=at + 1)

    def parse_group(stop: str = "") -> "dict[str, int]":
        nonlocal pos
        counts: dict[str, int] = {}
        group_start = pos
        while pos < n and text[pos] != stop:
            ch = text[pos]
            if ch == "(":
                open_pos = pos
                pos += 1
                inner = parse_group(")")
                if pos >= n:  # pragma: no cover - guarded by parse_group
                    raise error("括号未闭合", open_pos)
                pos += 1  # consume ')'
                mult = read_subscript()
                for symbol, amount in inner.items():
                    counts[symbol] = counts.get(symbol, 0) + amount * mult
            elif ch == ")":
                raise error("存在未匹配的右括号", pos)
            elif ch.isalpha() and ch.isupper():
                symbol, sym_pos = read_symbol()
                if known_elements is not None and symbol not in known_elements:
                    raise error(f"未知元素 '{symbol}'", sym_pos)
                mult = read_subscript()
                counts[symbol] = counts.get(symbol, 0) + mult
            elif ch.isdigit():
                raise error("数字下标前缺少元素或括号", pos)
            else:
                raise error(f"分子式中出现非法字符 {ch!r}", pos)
        if pos == group_start:
            raise error("括号内容不能为空", group_start)
        return counts

    def read_symbol() -> "tuple[str, int]":
        nonlocal pos
        start = pos
        pos += 1
        while pos < n and text[pos].isalpha() and text[pos].islower():
            pos += 1
        return text[start:pos], start

    def read_subscript() -> int:
        nonlocal pos
        start = pos
        digits = ""
        while pos < n and text[pos].isdigit():
            digits += text[pos]
            pos += 1
        if digits == "":
            return 1
        value = int(digits)
        if value == 0:
            raise error("下标不得为零", start)
        return value

    counts = parse_group()
    if pos != n:
        raise error("括号未闭合", pos)
    return counts
