from __future__ import annotations
from contextlib import asynccontextmanager
from pathlib import Path
from typing import TypedDict
from fastapi import FastAPI
from pydantic import BaseModel, Field
import pickle
from arcs.traversal import Traversal
import chempy
import networkx as nx


class _ReactionDict(TypedDict):
    e: chempy.Equilibrium
    k: float
    g: float


DATA_DIR = Path(__file__).parent.parent / "arcs/dash_app/data"
REACTIONS: dict[float, dict[float, dict[int, _ReactionDict]]] = {}
GRAPHS: dict[float, dict[float, nx.MultiDiGraph]] = {}


@asynccontextmanager
async def lifespan(app: FastAPI):
    global REACTIONS
    global GRAPH

    with open(DATA_DIR / "SCAN_reactions.p", "rb") as stream:
        REACTIONS = pickle.load(stream)
    with open(DATA_DIR / "SCAN_graph.p", "rb") as stream:
        GRAPH = pickle.load(stream)

    yield


app = FastAPI(lifespan=lifespan)


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

    traversal = Traversal(
        graph=GRAPH,
        reactions=REACTIONS,
    )
    traversal.run(
        trange=[form.temperature],
        prange=[form.pressure],
        ic=form.concs,
        sample_length=10,
    )

    results = {"initfinaldiff": traversal.initfinaldiff, "data": traversal.data}
    return results
