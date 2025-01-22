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



if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="localhost", port=8000)
