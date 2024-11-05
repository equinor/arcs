from __future__ import annotations
from fastapi import FastAPI
from pydantic import BaseModel, Field
from arcs.traversal import Traversal


app = FastAPI()


class SimulationRequest(BaseModel):
    temperature: float
    pressure: float
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


@app.post("/run_simulation")
async def run_simulation(form: SimulationRequest):
    global GRAPH
    global REACTIONS

    traversal = Traversal()
    traversal.run(
        trange=[form.temperature],
        prange=[form.pressure],
        ic=form.concs,
        sample_length=10,
    )

    results = {"initfinaldiff": traversal.initfinaldiff, "data": traversal.data}
    return results
