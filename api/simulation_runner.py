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
    import argparse
    import json

    parser = argparse.ArgumentParser(description="Run the simulation.")
    parser.add_argument("--temperature", type=float, required=True, help="Temperature value.")
    parser.add_argument("--pressure", type=float, required=True, help="Pressure value.")
    parser.add_argument("--concs", type=json.loads, required=True, help="Concentrations in JSON format.")
    parser.add_argument("--samples", type=int, required=True, help="Sample length.")
    args = parser.parse_args()

    form = SimulationRequest(
        temperature=args.temperature,
        pressure=args.pressure,
        concs=args.concs,
        samples=args.samples,
    )
    output = run_simulation(form)
    print(json.dumps(output, default=str))
