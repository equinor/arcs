# test_app.py
from fastapi.testclient import TestClient
from api.app import app

client = TestClient(app)


def test_run_simulation():
    payload = {
        "trange": [300],
        "prange": [10],
        "concs": {
            "CO2": 1,
            "H2O": 2e-05,
            "H2S": 3e-05,
            "SO2": 1e-05,
            "NO2": 5e-05,
        },
        "settings": {
            "nprocs": 1,
            "sample_length": 320,
            "max_rank": 10,
            "max_compounds": 5,
            "probability_threshold": 0.1,
            "path_depth": 5,
            "ceiling": 2000,
            "scale_highest": 0.2,
            "rank_small_reactions_higher": True,
        },
    }
    response = client.post("/run_simulation", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert "initfinaldiff" in data
    assert "data" in data
