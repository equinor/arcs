"""
Ensure that the model/dataset.nc file is consistent
"""
from typing import Iterable
import pytest
import chempy
import numpy as np
import xarray as xr
from arcs.model import MODEL_PATH


@pytest.fixture
def ds() -> Iterable[xr.Dataset]:
    with xr.open_dataset(MODEL_PATH / "dataset.nc", engine="h5netcdf") as ds:
        yield ds


def test_has_variables(ds: xr.Dataset):
    assert set(ds.variables.keys()) == {"gibbs", "pressure", "temperature", "reaction"}


def test_has_units(ds: xr.Dataset):
    assert ds["gibbs"].attrs["unit"] == "eV"
    assert ds["temperature"].attrs["unit"] == "K"
    assert ds["pressure"].attrs["unit"] == "bar A"


def test_reactions_are_valid(ds: xr.Dataset):
    for value in ds["reaction"].values:
        # Chempy will check that it's balanced and such
        eq = chempy.Equilibrium.from_string(value)
        assert value == str(eq), value


def test_has_valid_substances(ds: xr.Dataset):
    assert len(ds.attrs["substances"]) > 0
    for value in ds.attrs["substances"]:
        s = chempy.Substance.from_formula(value)
        assert s.composition != {}, value
        assert s.charge == 0, f"{value} charge is not 0. (Consider updating this test)"


def test_gibbs_is_finite(ds: xr.Dataset):
    assert np.isfinite(ds["gibbs"]).all()


@pytest.mark.parametrize("key", ["temperature", "pressure"])
def test_is_sorted(ds: xr.Dataset, key: str):
    array = ds[key].values
    assert (array[:-1] <= array[1:]).all()
