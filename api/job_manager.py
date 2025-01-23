import json
import httpx
from fastapi import HTTPException
from api.models import SimulationRequest


class JobManager:
    def __init__(self, baseurl: str):
        self.baseurl = baseurl
        self.client = httpx.AsyncClient(base_url=baseurl)

    async def start_job(self, form: SimulationRequest):
        payload = {
            "payload": json.dumps(form.model_dump()),
            "resources": {
                "requests": {"memory": "4Gi", "cpu": "2000m"},
                "limits": {"memory": "8Gi", "cpu": "4000m"},
            },
        }

        try:
            r = await self.client.post("/api/v1/jobs", json=payload)
            r.raise_for_status()
        except httpx.HTTPError as e:
            raise HTTPException(
                status_code=500, detail=f"Could not create job: {str(e)}"
            )
        return r.json()

    async def cancel_job(self, job_id: str):
        try:
            r = await self.client.post(f"/api/v1/jobs/{job_id}/stop")
            r.raise_for_status()
        except httpx.HTTPError as e:
            raise HTTPException(
                status_code=500, detail=f"Could not cancel job: {str(e)}"
            )
        return r.json()

    async def get_job_status(self, job_id: str):
        try:
            r = await self.client.get(f"/api/v1/jobs/{job_id}")
            r.raise_for_status()
        except httpx.HTTPError as e:
            raise HTTPException(
                status_code=500, detail=f"Could not get job status: {str(e)}"
            )
        return r.json()
