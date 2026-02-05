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
from arcs.traversal import Traversal
from arcs.generate import GenerateInitialConcentrations
from arcs.generate import GraphGenerator
from arcs.analysis import AnalyseSampling

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
    graph = GraphGenerator().from_file(
        filename='../app/data/dft_data.json',
        temperature=form.temperature,
        pressure=form.pressure,
        max_reaction_length=5
    )

    gic = GenerateInitialConcentrations(graph=graph).update_ic(form.concs)
    
    t = Traversal(graph=graph)
    
    data = t.sample(initial_concentrations=gic,ncpus=4,nsamples=1000)    
    analysis = AnalyseSampling()
    stats = pd.Series(analysis.reaction_statistics(data)).sort_values(ascending=False)

    average_data = pd.DataFrame(analysis.average_sampling(data))
    average_data = average_data.loc[~(average_data==0).all(axis=1)]
    average_data.sort_values(by='diff',inplace=True)

    result_stats = pd.DataFrame(
        {
            "comps": list(average_data.index),
            "values": average_data["mean"].values,
            "variance": average_data["var"].values,
            "variance_minus": -average_data["var"].values,
        }
    )

    return {"results": results, "analysis": analysis, "chart_data": result_stats}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="localhost", port=8002)
