#!/usr/env/bin python3
import warnings
import pickle
from tabulate import tabulate
from arcs.setup_functions import GraphGenerator
from arcs.traversal import Traversal
from monty.serialization import loadfn

if __name__ == "__main__":
    warnings.simplefilter("ignore")

    print("""   \   _ \  __|  __| 
  _ \    / (   \__ \ 
_/  _\_|_\\___|____/""")

    functional = "pbe"

    print("\n loading reactions...")

    reactions = pickle.load(
        open("../../reaction_data/{}_reactions.p".format(functional), "rb")
    )

    trange = [300]
    prange = [98]

    gs = GraphGenerator(reactions).multidigraph_from_t_and_p_range(trange, prange)

    ic = loadfn("../initial_concentrations.json")

    headers = ["compound", "conc /ppm"]

    print(tabulate(ic.items(), headers=headers))

    sample_length = 10000
    path_depth = 200

    print(
        "starting traversal... \n  -| path depth = {} \n -| samples = {}".format(
            path_depth, sample_length
        )
    )

    sets = Traversal(gs, reactions, ic, trange, prange)
    t_and_p_data = sets.graph_sampling_processes(
        sample_length=sample_length,
        path_depth=path_depth,
        nprocs=32,
        random_path_depth=True,
    )

    print("Writing to file...")

    pickle.dump(
        t_and_p_data,
        open(
            "concs_and_stats_{}_{}_{}.p".format(functional, sample_length, path_depth),
            "wb",
        ),
    )

    print(
        "Data written to concs_and_stats_{}_{}_{}.p".format(
            functional, sample_length, path_depth
        )
    )
