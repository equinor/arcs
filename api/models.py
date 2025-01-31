from typing import Optional
from pydantic import BaseModel, Field


class ResourceRequests(BaseModel):
    memory: Optional[str] = "8Gi"
    cpu: Optional[str] = "4000m"


class ResourceLimits(BaseModel):
    memory: Optional[str] = "16Gi"
    cpu: Optional[str] = "8000m"


class Resources(BaseModel):
    requests: ResourceRequests
    limits: ResourceLimits


class SimulationRequest(BaseModel):
    project_id: Optional[str] = None
    scenario_id: Optional[str] = None
    temperature: int
    pressure: int
    concs: dict[str, float] = Field(default_factory=dict)
    samples: int

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "temperature": 300,
                    "pressure": 10,
                    "concs": {
                        "SO2": 10e-6,
                        "NO2": 50e-6,
                        "H2S": 30e-6,
                        "H2O": 20e-6,
                    },
                    "samples": 10,
                }
            ]
        }
    }
