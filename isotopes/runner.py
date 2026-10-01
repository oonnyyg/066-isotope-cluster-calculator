"""Turn validated input JSON into cluster results (no filesystem writes here)."""

from .cluster import compute_cluster
from .elements import ElementTable
from .errors import IsotopeError
from .formula import parse_formula


def load_spec(data: object) -> "tuple[ElementTable, list[dict]]":
    """Validate the parsed input document and build computation jobs."""
    if not isinstance(data, dict):
        raise IsotopeError("输入根节点必须是 JSON 对象")
    table = ElementTable(data.get("elements"))

    charge = data.get("charge")
    if not isinstance(charge, int) or isinstance(charge, bool) or charge <= 0:
        raise IsotopeError("电荷数必须是正整数", where="charge")

    formulas = data.get("formulas")
    single = isinstance(formulas, str)
    if single:
        formulas = [formulas]
    if not isinstance(formulas, list) or not formulas:
        raise IsotopeError("formulas 必须是非空字符串或非空数组", where="formulas")

    jobs = []
    names: set[str] = set()
    for index, entry in enumerate(formulas):
        if single:
            name = "molecule"
            formula = entry
            if not isinstance(formula, str):
                raise IsotopeError("formulas 必须是分子式字符串", where="formulas")
        else:
            if not isinstance(entry, dict):
                raise IsotopeError("每个分子式条目必须是对象", where=f"formulas[{index}]")
            formula = entry.get("formula")
            name = entry.get("name", f"molecule_{index + 1}")
            if not isinstance(formula, str):
                raise IsotopeError("formula 必须是字符串", where=f"formulas[{index}]")
            if not isinstance(name, str) or not name:
                raise IsotopeError("name 必须是非空字符串", where=f"formulas[{index}]")
        if name in names:
            raise IsotopeError(f"分子式名称 {name!r} 重复", where=f"formulas[{index}]")
        names.add(name)
        counts = parse_formula(formula, table.symbols(), where=f"分子式 {name}")
        jobs.append({"name": name, "formula": formula, "counts": counts})
    return table, jobs


def run(data: object) -> list:
    table, jobs = load_spec(data)
    charge = data["charge"]
    clusters = []
    for job in jobs:
        clusters.append(
            compute_cluster(job["name"], job["formula"], job["counts"], table, charge)
        )
    return clusters
