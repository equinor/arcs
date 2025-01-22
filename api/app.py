from __future__ import annotations
import os
import pandas as pd
from fastapi import Depends, FastAPI
from fastapi.middleware.cors import CORSMiddleware
from api.authentication import authenticated_user_claims
from api.models import SimulationRequest
from arcs.analysis import AnalyseSampling
from arcs.traversal import traverse
from dotenv import load_dotenv
from api.simulation_runner import run_simulation
import httpx

load_dotenv()
app = FastAPI(dependencies=[Depends(authenticated_user_claims)])
app.swagger_ui_init_oauth = {
    "clientId": os.environ.get("CLIENT_ID"),
    "appName": "ARCS API",
    "usePkceWithAuthorizationCodeGrant": True,  # Enable PKCE
    "scope": os.environ.get("API_SCOPE"),
}

origins = [
    "http://localhost:5173",
    "https://frontend-acidwatch-dev.radix.equinor.com",
    "https://acidwatch.radix.equinor.com",
    "https://frontend-acidwatch-prod.radix.equinor.com",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)




@app.post("/run_simulation")
def run_simulation(form: SimulationRequest):


    result = run_simulation(form)
    return result

@app.post("/run_simulation_radix_job")
def run_simulation(form: SimulationRequest):
    url ="http://runsimulation:8000/api/v1"
    payload_path = "/runsimulation/args"
    params = {
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

    payload = {
        "payload": params
    }
    with httpx.Client() as client:
        response = client.post(
            url,
            json=payload,
            timeout=60.0)

        return response.json()


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="localhost", port=8000)
