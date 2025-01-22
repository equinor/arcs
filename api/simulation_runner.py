import pandas as pd
from api.models import SimulationRequest
from arcs.analysis import AnalyseSampling
from arcs.traversal import traverse


def run_simulation(form: SimulationRequest):
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

    df_mean = pd.DataFrame(analysis.mean_data).T
    df_mean = df_mean[df_mean["value"] != 0]
    result_stats = pd.DataFrame(
        {
            "comps": list(df_mean.T.keys()),
            "values": df_mean["value"].values,
            "variance": df_mean["variance"].values,
            "variance_minus": -df_mean["variance"].values,
        }
    )

    return {"results": results, "analysis": analysis, "chart_data": result_stats}


if __name__ == "__main__":
    import json

    with open("/runsimulation/args/payload", "r", encoding="utf8") as f:
        payload = json.load(f)
        print(f"payload: {payload}")
        form = SimulationRequest(
            temperature=payload["temperature"],
            pressure=payload["pressure"],
            concs=payload["concs"],
            samples=payload["samples"],
        )
        output = run_simulation(form)
        print(json.dumps(output, default=str))
