import os
from fastapi import FastAPI
from pydantic import BaseModel
import pickle
from arcs.traversal import Traversal

app = FastAPI()


class SimulationRequest(BaseModel):
    trange: list = [300]
    prange: list = [10]
    concs: dict = {}
    settings: dict = {}


@app.post("/run_simulation")
async def run_simulation(request: SimulationRequest):
    trange = request.trange
    prange = request.prange
    concs = request.concs
    settings = request.settings

    graph = pickle.load(
        open(os.path.join(os.path.dirname(__file__), "../app/data/SCAN_graph.p"), "rb")
    )
    traversal = Traversal(
        graph=graph,
        reactions=os.path.join(
            os.path.dirname(__file__), "../app/data/SCAN_reactions.p"
        ),
    )
    traversal.run(trange=trange, prange=prange, save=False, ic=concs, **settings)

    results = {"initfinaldiff": traversal.initfinaldiff, "data": traversal.data}
    # TODO: return graph data?
    return results


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="127.0.0.1", port=8000)
