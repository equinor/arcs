from __future__ import annotations

import os

import pandas as pd
from azure.monitor.opentelemetry import configure_azure_monitor
from dotenv import load_dotenv
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from opentelemetry import trace
from opentelemetry.instrumentation.fastapi import FastAPIInstrumentor
from opentelemetry.instrumentation.httpx import HTTPXClientInstrumentor
from opentelemetry.trace import get_tracer_provider
from pydantic import BaseModel, Field

from arcs.analysis import AnalyseSampling
from arcs.traversal import traverse

tracer = trace.get_tracer(__name__, tracer_provider=get_tracer_provider())

load_dotenv()
app = FastAPI()

if os.getenv("APPLICATIONINSIGHTS_CONNECTION_STRING"):
    configure_azure_monitor(
        connection_string=os.getenv("APPLICATIONINSIGHTS_CONNECTION_STRING")
    )

HTTPXClientInstrumentor().instrument()
FastAPIInstrumentor.instrument_app(app)


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


class SimulationRequest(BaseModel):
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


@app.post("/run_simulation")
def run_simulation(form: SimulationRequest):
    results = traverse(
        form.temperature,
        form.pressure,
        form.concs,
        samples=form.samples,
        nproc=4,  # We hardcode to avoid exhausting too early. Should be a better solution
    )

    analysis = AnalyseSampling(results.data)
    analysis.reaction_statistics()
    analysis.mean_sampling()
    analysis.reaction_paths()

    df_m_t = pd.DataFrame(analysis.mean_data).T
    df_m_t = df_m_t[df_m_t["value"] != 0]
    result_stats = pd.DataFrame(
        {
            "comps": list(df_m_t.T.keys()),
            "values": df_m_t["value"].values,
            "variance": df_m_t["variance"].values,
            "variance_minus": -df_m_t["variance"].values,
        }
    )

    return {"results": results, "analysis": analysis, "chart_data": result_stats}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="localhost", port=8002)
