from __future__ import annotations
from fastapi import FastAPI
from pydantic import BaseModel, Field
import pandas as pd
from arcs.analysis import AnalyseSampling
from arcs.traversal import traverse
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

origins = [
    "http://localhost:5173",
    "https://frontend-acidwatch-dev.radix.equinor.com",
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
async def run_simulation(form: SimulationRequest):
    results = traverse(
        form.temperature,
        form.pressure,
        form.concs,
        sample_length=form.samples,
    )

    analysis = AnalyseSampling(results.data, markdown=True)
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

    uvicorn.run(app, host="127.0.0.1", port=8000)
