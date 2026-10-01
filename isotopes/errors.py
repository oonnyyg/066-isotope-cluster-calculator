"""Errors carrying a location for invalid input."""


class IsotopeError(ValueError):
    """Raised for invalid formulas or element tables.

    ``where`` identifies the offending molecule or element, and
    ``position`` (1-based) points at a character inside the formula.
    """

    def __init__(self, message: str, *, where: str = "", position: int | None = None):
        self.message = message
        self.where = where
        self.position = position
        loc = ""
        if where:
            loc = f"[{where}]"
            if position is not None:
                loc = f"[{where} 位置 {position}]"
        elif position is not None:
            loc = f"[位置 {position}]"
        super().__init__(f"{loc} {message}".strip())
