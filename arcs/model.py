from __future__ import annotations

from collections import defaultdict
from itertools import product
import re
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict
import chempy
import numpy as np
import xarray as xr
from scipy import interpolate  # type: ignore


if TYPE_CHECKING:
    import numpy.typing as npt


type Table = dict[tuple[str, str], tuple[npt.NDArray[np.int16], npt.NDArray[np.float64]]]

MODEL_PATH = Path(__file__).parent.parent / "model"
BOLTZMANN_CONSTANT = np.float64(8.617333262e-5)  # Boltzmann constant in eV K⁻¹


class ReactionType(TypedDict):
    e: chempy.Equilibrium
    k: float
    g: float


DATASET = xr.open_dataset(MODEL_PATH / "dataset.nc", engine="h5netcdf")
TEMPERATURES: npt.NDArray[np.int64] = DATASET["temperature"].values
PRESSURES: npt.NDArray[np.int64] = DATASET["pressure"].values
SUBSTANCES: set[str] = set(DATASET.attrs["substances"])


_PARSER_PATTERN = re.compile(r"(\d+ )?((:?[A-Z][a-z]?\d*)+)")
_MOLECULE_PATTERN = re.compile(r"[A-Z][a-z]?(\d+)?")


def _parse_side(string: str) -> dict[str, int]:
    return {
        m[2]: int(m[1]) if m[1] else 1
        for m in _PARSER_PATTERN.finditer(string)
    }


def _parse_reaction(string: str) -> chempy.Equilibrium:
    reac, prod = string.split("=")
    return chempy.Equilibrium(
        _parse_side(reac), _parse_side(prod), checks=()
    )


@lru_cache
def _num_atoms(name: str) -> int:
    return sum(int(x) if x else 1 for x in _MOLECULE_PATTERN.findall(name))


def _cost_function(
    gibbs: float,
    temperature: float,
    compounds: dict[str, int],
) -> float:
    num_atoms = sum(_num_atoms(k) * v for k, v in compounds.items())
    return np.log(1 + (273 / temperature) * np.exp(gibbs / num_atoms))


@lru_cache
def get_reactions() -> list[chempy.Equilibrium]:
    return list(map(_parse_reaction, DATASET["reaction"].values))


@lru_cache
def get_gibbs(temperature: int, pressure: int) -> npt.NDArray[np.float64]:
    return DATASET["gibbs"].sel(temperature=temperature, pressure=pressure).values.ravel()


@lru_cache
def get_table(
    temperature: int, pressure: int, *, rank_small_reactions_higher: bool = True
) -> Table:
    gibbs = get_gibbs(temperature, pressure)
    table: dict[tuple[str, str], list[tuple[int, np.float64]]] = defaultdict(list)

    for i, (g, eq) in enumerate(zip(gibbs, get_reactions())):
        substances = set(eq.reac) | set(eq.prod)

        mult = 10 ** len(substances) if rank_small_reactions_higher else 1

        f_cost = np.float64(1 / _cost_function(g, temperature, eq.reac) * mult)
        b_cost = np.float64(1 / _cost_function(-g, temperature, eq.prod) * mult)
        for cost, srcs in [(f_cost, eq.reac), (b_cost, eq.prod)]:
            for src, dst in product(srcs, substances):
                if src == dst:
                    continue
                table[(src, dst)].append((i, cost))

    new_table: Table = {}
    for key, val in table.items():
        new_index = np.array([i for i, _ in val])
        new_gibbs = np.array([g for _, g in val])
        sorted_indices = np.argsort(new_gibbs)[::-1]
        new_table[key] = (
            new_index[sorted_indices].astype(np.int16),
            new_gibbs[sorted_indices],
        )

    return new_table


def get_equilibrium_constants(temperature: int, pressure: int) -> npt.NDArray[np.float64]:
    gibbs = get_gibbs(temperature, pressure)
    return np.exp(gibbs / (BOLTZMANN_CONSTANT * temperature))


def _interp_gibbs(
    temperature: int,
    pressure: int,
) -> npt.NDArray[np.float64]:
    """
    Calculates Gibbs contant using linear interpolation between 2 enclosing points

    Returns: Data required for reaction table
    """
    t_lower, t_upper = _find_enclosing(TEMPERATURES, temperature)
    p_lower, p_upper = _find_enclosing(PRESSURES, pressure)

    # Keeping pressure constant
    gibbs_p_lower = DATASET["gibbs"].sel(temperature=[t_lower, t_upper], pressure=p_lower).values.transpose()
    gibbs_p_upper = DATASET["gibbs"].sel(temperature=[t_lower, t_upper], pressure=p_upper).values.transpose()

    # Interpolate Gibbs values for temperature
    g_val_low_pressure = interpolate_gibbs_values(
        gibbs_p_lower, t_lower, t_upper, temperature
    )
    g_val_high_pressure = interpolate_gibbs_values(
        gibbs_p_upper, t_lower, t_upper, temperature
    )

    # Prepare for pressure interpolation
    gibbs_values_for_changing_pressure = np.array(
        [g_val_low_pressure, g_val_high_pressure]
    ).transpose()

    # Interpolate Gibbs values for pressure
    calculated_gibbs_values = interpolate_gibbs_values(
        gibbs_values_for_changing_pressure, p_lower, p_higher, pressure
    )

    return calculated_gibbs_values


def interpolate_gibbs_values(
    values: npt.NDArray[np.float64], point1: int, point2: int, target: int
) -> npt.NDArray[np.float64]:
    interpolation = [interpolate.interp1d([point1, point2], y) for y in values]
    return np.array([f(target) for f in interpolation])


def _find_enclosing(a: npt.NDArray[np.int64], x: int) -> tuple[np.int64, np.int64]:
    index = np.searchsorted(a, x)
    try:
        y = a[index]
        if index == 0 and x != y:
            raise IndexError()
        return (y, y) if x == y else (a[index - 1], y)
    except IndexError:
        raise IndexError(f"Search value must be between {a[0]} and {a[-1]}, is {x}")


@lru_cache
def get_reaction_compounds() -> dict[int, set[str]]:
    return {k: set(r.reac) | set(r.prod) for k, r in enumerate(get_reactions())}
