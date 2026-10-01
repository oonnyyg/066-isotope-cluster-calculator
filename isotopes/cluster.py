"""Independent-isotope convolution into an integer nominal-mass distribution."""

from dataclasses import dataclass

from .elements import ElementTable


@dataclass
class Peak:
    mass_number: int       # summed integer mass number
    probability: float     # probability of this nominal mass
    mean_mass: float       # probability-weighted mean exact molecular mass (u)
    charge: int = 1

    @property
    def mz(self) -> float:
        return self.mean_mass / self.charge


@dataclass
class Cluster:
    name: str
    formula: str
    charge: int
    counts: dict
    peaks: list

    @property
    def total_probability(self) -> float:
        return sum(peak.probability for peak in self.peaks)

    @property
    def mean_molecular_mass(self) -> float:
        return sum(peak.probability * peak.mean_mass for peak in self.peaks)


def _convolve(dist, mass_number, exact_mass, probability):
    """Merge one isotope choice into *dist* (mass_number -> [prob, mass_sum])."""
    merged = {}
    for nominal, (prob, mass_sum) in dist.items():
        key = nominal + mass_number
        new_prob = prob * probability
        new_mass_sum = (mass_sum + exact_mass * prob) * probability
        if key in merged:
            old_prob, old_mass_sum = merged[key]
            merged[key] = (old_prob + new_prob, old_mass_sum + new_mass_sum)
        else:
            merged[key] = (new_prob, new_mass_sum)
    return merged


def compute_cluster(name: str, formula: str, counts: dict, table: ElementTable,
                    charge: int) -> Cluster:
    """Compute the nominal-mass distribution assuming independent atoms.

    Combinations sharing a total mass number are merged; the stored exact
    mass is the probability-weighted mean over those combinations.
    """
    dist = {0: (1.0, 0.0)}
    for symbol, atom_count in counts.items():
        for _ in range(atom_count):
            choices = table.elements[symbol]
            nxt = {}
            for isotope in choices:
                step = _convolve(
                    dist,
                    isotope.mass_number,
                    isotope.exact_mass,
                    isotope.abundance,
                )
                for key, (prob, mass_sum) in step.items():
                    if key in nxt:
                        old_prob, old_mass_sum = nxt[key]
                        nxt[key] = (old_prob + prob, old_mass_sum + mass_sum)
                    else:
                        nxt[key] = (prob, mass_sum)
            dist = nxt

    peaks = []
    for nominal in sorted(dist):
        prob, mass_sum = dist[nominal]
        if prob > 0.0:
            peaks.append(Peak(nominal, prob, mass_sum / prob, charge))
    return Cluster(name, formula, charge, counts, peaks)
