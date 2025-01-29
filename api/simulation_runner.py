import os
import pandas as pd
from api.models import SimulationRequest
from arcs.analysis import AnalyseSampling
from arcs.traversal import traverse
import json
import httpx
from azure.identity import DefaultAzureCredential


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

    return {
        "results": results.to_dict(),
        "analysis": analysis.to_dict(),
        "chart_data": result_stats.to_dict(orient="records"),
    }

if __name__ == "__main__":

    scope = os.environ.get("API_SCOPE")
    credential = DefaultAzureCredential()
    token = credential.get_token(scope)

    with open("/runsimulation/args/payload", "r", encoding="utf8") as f:
        payload = json.load(f)
        project_id = payload["project_id"]
        scenario_id = payload["scenario_id"]
        form = SimulationRequest(
            temperature=payload["temperature"],
            pressure=payload["pressure"],
            concs=payload["concs"],
            samples=payload["samples"],
        )
        result = run_simulation(form)

        post_data = {
            "scenario_id": scenario_id,
            "raw_results": json.dumps(result),
       }

        httpx.post(
            f"https://backend-acidwatch-dev.radix.equinor.com/project/{project_id}/scenario/{scenario_id}/result",
            json=post_data,
            headers={"Authorization": "Bearer " + token.token},
        )





