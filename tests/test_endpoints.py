from __future__ import annotations

from typing import Iterable

import pytest
from fastapi.testclient import TestClient

from api.app import app


@pytest.fixture
def client() -> Iterable[TestClient]:
    with TestClient(app) as client:
        yield client


def test_run_simulation(client):
    payload = {
        "temperature": 300,
        "pressure": 10,
        "concs": {
            "CO2": 1,
            "H2O": 2e-05,
            "H2S": 3e-05,
            "SO2": 1e-05,
            "NO2": 5e-05,
        },
        "samples": 10,
    }

    response = client.post("/run_simulation", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert "results" in data
    assert "analysis" in data
