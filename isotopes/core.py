"""Isotope cluster convolution over integer total mass numbers."""

from __future__ import annotations

from .errors import IsotopeError
from .table import Isotope


def compute_cluster(
    counts: dict[str, int],
    table: dict[str, list[Isotope]],
) -> dict[int, tuple[float, float]]:
    """Return {total_mass_number: (probability, prob_weighted_exact_mass_sum)}."""
    dist: dict[int, list[float]] = {0: [1.0, 0.0]}
    for symbol, count in counts.items():
        if symbol not in table:
            raise IsotopeError(f"分子式中存在未知元素 {symbol!r}")
        for _ in range(count):
            dist = _convolve_atom(dist, table[symbol])
    return {k: (v[0], v[1]) for k, v in dist.items() if v[0] > 0.0}


def _convolve_atom(
    dist: dict[int, list[float]], isotopes: list[Isotope]
) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    for mass_num, (prob, mass_sum) in dist.items():
        for iso in isotopes:
            if iso.abundance == 0.0:
                continue
            key = mass_num + iso.mass_number
            p = prob * iso.abundance
            slot = out.setdefault(key, [0.0, 0.0])
            slot[0] += p
            slot[1] += mass_sum * iso.abundance + p * iso.exact_mass
    return out


def summarize(
    cluster: dict[int, tuple[float, float]], charge: int
) -> dict:
    total_prob = sum(p for p, _ in cluster.values())
    avg_mass = sum(m for _, m in cluster.values())
    peaks = [
        {
            "mass_number": k,
            "probability": p,
            "mean_exact_mass": m / p,
            "mz": (m / p) / charge,
        }
        for k, (p, m) in sorted(cluster.items())
    ]
    return {
        "total_probability": total_prob,
        "average_mass": avg_mass,
        "average_mz": avg_mass / charge,
        "peaks": peaks,
    }
