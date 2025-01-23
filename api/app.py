from __future__ import annotations
import os

from fastapi import Depends, FastAPI, APIRouter
from fastapi.middleware.cors import CORSMiddleware
from api.authentication import authenticated_user_claims
from api.job_manager import JobManager
from api.models import SimulationRequest
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


job_manager = JobManager(baseurl="http://runsimulation:8000")
router = APIRouter()


@router.post("/run_simulation")
def run_simulation_endpoint(form: SimulationRequest):
    result = run_simulation(form)
    return result


@router.post("/start_radix_job")
async def start_job_endpoint(form: SimulationRequest):
    print(f"Starting job: {form}")
    return await job_manager.start_job(form)


@router.post("/cancel_radix_job")
async def cancel_job_endpoint(job_id: str):
    return await job_manager.cancel_job(job_id)


@router.get("/get_radix_job_status")
async def get_job_status_endpoint(job_id: str):
    return await job_manager.get_job_status(job_id)


app.include_router(router)

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="localhost", port=8000)
