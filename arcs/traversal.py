import random
from chempy.equilibria import Equilibrium, EqSystem
from chempy import Substance
import copy
import networkx as nx
import itertools as it
from tqdm import tqdm
from numpy.random import choice
import platform
import psutil

import warnings

from datetime import datetime
import math
import numpy as np
import pandas as pd
from pathos.helpers import mp as pmp
import pickle


class Traversal:
    def __init__(self, graph, reactions):
        if isinstance(graph, str):
            self.graph = pickle.load(open(graph, "rb"))
        else:
            self.graph = graph
        if isinstance(reactions, str):
            self.reactions = pickle.load(open(reactions, "rb"))
        else:
            self.reactions = reactions
        # self.concs = copy.deepcopy(concs) # saves having to reload the class each time
        self.trange = list(self.graph)
        self.prange = list(self.graph[self.trange[0]])

        # default values:
        self.co2 = False
        self.max_compounds = 5
        self.probability_threshold = 0.05
        self.max_rank = 5
        self.sample_length = 1000
        self.path_depth = 20
        self.random_path_depth = False
        self.nprocs = 4
        self.ceiling = 2000
        self.scale_highest = 0.1
        self.rank_small_reactions_higher = True
        self.method = "Bellman-Ford"
        self.final_concs = {}
        self.initfinaldiff = {}

    def _get_weighted_random_compounds(
        self,
        T,
        P,
        init_concs=None,
        co2=False,  # should probably be "exclude_co2"
        max_compounds=5,
        probability_threshold=0.05,
        scale_highest=0.1,  # how much to scale the highest components
        ceiling=3000,
    ):  # ceiling percent larger than the median average
        [n for n in self.graph[T][P].nodes() if isinstance(n, str)]
        concs = copy.deepcopy(init_concs)  # don't modify the original
        if co2 is False:
            del concs["CO2"]  # CO2 will always be too large as it is the background
        # house keeping:
        num_not_zero = len([x for x in concs.values() if x > 0])
        if max_compounds > num_not_zero:
            max_compounds = num_not_zero

        # scale the probabilities accordingly based upon a ceiling percentage

        median_conc = np.median(
            [v for v in concs.values() if v > 0]
        )  # median > mean for this
        # new_concs = {}
        above_ceiling = {
            k: v for k, v in concs.items() if v > (median_conc * (1 + (ceiling / 100)))
        }
        # modify the ceiling by scaling it down to a suitable value
        # should still max out if concentrations become way to high
        for k, v in above_ceiling.items():
            concs[k] = v * scale_highest

        # get the probabilities based upon relative concentrations:
        p_1 = {k: v / sum(concs.values()) for k, v in concs.items()}
        # now filter based upon the probability threshold:
        p_2 = {k: v for k, v in p_1.items() if v > probability_threshold}
        p_3 = {k: v / sum(p_2.values()) for k, v in p_2.items()}
        # make a list of choices based upon the probabilities
        available = list(
            np.random.choice(list(p_3.keys()), 100, p=list(p_3.values()))
        )  # make this list length of the nodes
        # now make a list max_compounds long of random choices based on available
        choices = {}
        for c in range(max_compounds):
            if c == 0:
                c1 = np.random.choice(available)
                choices[c1] = p_3[c1]
            else:
                try:
                    for i in range(available.count(list(choices)[c - 1])):
                        available.remove(list(choices)[c - 1])
                    try:
                        c2 = np.random.choice(available)
                        choices[c2] = p_3[c2]
                    except Exception:
                        pass
                except Exception:
                    pass

        return choices

    def _length_multiplier(self, candidate):
        if self.rank_small_reactions_higher:
            return len(list(candidate))
        else:
            return 1

    def _get_weighted_reaction_rankings(
        self, T, P, choices, max_rank=20, method="Bellman-Ford"
    ):
        rankings = {}
        if len(choices) > 1:
            possibilities = list(
                nx.shortest_paths.all_shortest_paths(
                    self.graph[T][P], list(choices)[0], list(choices)[1], method=method
                )
            )

            for x in possibilities:
                candidates = list(self.graph[T][P][x[1]])
                if len(choices) > 2:
                    for c in list(choices)[2:]:
                        if c in candidates:
                            weight = self.graph[T][P].get_edge_data(x[0], x[1])[0][
                                "weight"
                            ] * 10 ** self._length_multiplier(self.graph[T][P][x[1]])
                            rankings[x[1]] = {
                                "candidates": candidates,
                                "weight": weight,
                            }
                else:
                    weight = self.graph[T][P].get_edge_data(x[0], x[1])[0][
                        "weight"
                    ] * 10 ** self._length_multiplier(self.graph[T][P][x[1]])
                    rankings[x[1]] = {"candidates": candidates, "weight": weight}
        if rankings:
            sorted_rankings = (
                pd.DataFrame(rankings).sort_values(by="weight", axis=1).to_dict()
            )
            topranks = [
                x for i, x in enumerate(sorted_rankings) if i <= max_rank
            ]  # need to sort first
            rankings = {x: rankings[x] for x in topranks}
            # return(pd.DataFrame(rankings).sort_values(by='weight',axis=1).to_dict()) # sort them
            # according to lowest weight first
            return rankings
        else:
            return None

    def _random_choice_unconnected(
        self, T, P, force_direct=False, co2=False
    ):  # currently randomly disjointed reactions that are weighted
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n, str)]
        if force_direct is True:
            pstring = [0, 1, 2]
            while len(pstring) > 2:
                source = self._get_weighted_random_compound(
                    T, P, co2=co2, force_selection=None
                )
                target = random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P], source, target, weight="weight")
                pstring = [n for n in p if isinstance(p, str)]
        else:
            source = self._get_weighted_random_compound(T, P)
            target = random.choice(nodes)
            p = nx.shortest_path(self.graph[T][P], source, target)
        return p

    def _random_choice_connected(
        self, T, P, force_direct=False, previous_index=None, co2=False
    ):  # this will be joined - I think we can also make a ranking of potential reactions based upon components in the stream as well
        if previous_index is None:
            raise ValueError("no previous compound selected")
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n, str)]
        if force_direct is True:
            pstring = [0, 1, 2]
            while len(pstring) > 2:
                present = [
                    c
                    for c in list(self.reactions[T][P][previous_index]["e"].reac)
                    + list(self.reactions[T][P][previous_index]["e"].prod)
                ]  # this should probably be weighted according to stoichiometry i.e. 2CO2 + H2O = [CO2, CO2, H2O]
                source = self._get_weighted_random_compound(
                    T, P, co2=co2, force_selection=present
                )
                target = random.choice(nodes)  # the next path will be random
                p = nx.shortest_path(self.graph[T][P], source, target, weight="weight")
                pstring = [n for n in p if isinstance(p, str)]
        else:
            source = self._get_weighted_random_compound(T, P)
            target = random.choice(nodes)
            p = nx.shortest_path(self.graph[T][P], source, target)
        return p

    def generate_eqsystem(self, index, T, P):
        charged_species = {
            "CO3H": -1,
            "NH4": +1,
            "NH2CO2": -1,
        }  # this needs to be added to the arguments
        rs = self.reactions[T][P][index]
        r = rs["e"].reac
        p = rs["e"].prod
        k = rs["k"]
        substances = {}
        for n in list(it.chain(*[list(r) + list(p)])):
            if n in list(charged_species.keys()):
                s = Substance.from_formula(n, **{"charge": charged_species[n]})
                substances[s.name] = s
            else:
                s = Substance.from_formula(
                    n, **{"charge": 0}
                )  # ,charge=0) # charge buggers up everything have removed for now....
                substances[s.name] = s
        eql = Equilibrium(reac=r, prod=p, param=k)
        try:
            return EqSystem(
                [eql], substances
            )  # might not just be able to try a return...
        except Exception:
            return None

    def equilibrium_concentrations(self, concs, eq):
        # something is going wrong here...
        fc = copy.deepcopy(concs)
        try:
            x, sol, sane = eq.root(fc)
            assert sol["success"] and sane
            for n, c in enumerate(x):
                fc[eq.substance_names()[n]] = c

            concs = fc
            eq = eq.string()
        except Exception:
            concs = fc
            eq = None
        return (concs, eq)

    def random_walk(
        self,
        T,
        P,
        probability_threshold=0.05,
        path_depth=50,
        # concs=None,
        max_compounds=5,
        max_rank=5,
        co2=False,
        scale_highest=1000,
        ceiling=3000,
        method="bellman-ford",
    ):
        final_concs = {0: copy.deepcopy(self.concs)}
        reactionstats = {0: None}

        for ip in range(1, path_depth + 1):
            fcs = copy.deepcopy(final_concs[ip - 1])
            try:
                choices = self._get_weighted_random_compounds(
                    T=T,
                    P=P,
                    init_concs=fcs,
                    max_compounds=max_compounds,
                    probability_threshold=probability_threshold,
                    co2=co2,
                    scale_highest=scale_highest,
                    ceiling=ceiling,
                )
            except Exception:
                path_depth = ip + 1
                break
            if len(choices) <= 1:  # not sure this is necessary....
                path_depth = ip + 1
                break
            rankings = self._get_weighted_reaction_rankings(
                T=T, P=P, choices=choices, max_rank=max_rank, method=method
            )
            if not rankings:
                break
            weights = {k: 1 / rankings[k]["weight"] for k in rankings}
            probabilities = {k: v / sum(weights.values()) for k, v in weights.items()}
            chosen_reaction = random.choice(
                [
                    choice(
                        list(probabilities.keys()),
                        len(probabilities),
                        p=list(probabilities.values()),
                    )
                ][0]
            )

            eqsyst = self.generate_eqsystem(chosen_reaction, T, P)
            # if reaction was previous reaction then break
            path_available = [r for r in reactionstats.values() if r is not None]
            if path_available:
                if (
                    eqsyst.string() == path_available[-1]
                    and eqsyst.string() == path_available[-1]
                ):
                    break  # extra break

            final_concs[ip], reactionstats[ip] = self.equilibrium_concentrations(
                fcs, eqsyst
            )

        return {
            "data": final_concs[list(final_concs)[-1]],
            "equation_statistics": [r for r in reactionstats.values() if r is not None],
            "path_length": len([r for r in reactionstats.values() if r is not None]),
        }

    def _queue_function(
        self,
        pbari,
        samples,
        T,
        P,
        probability_threshold,
        path_depth,
        max_compounds,
        max_rank,
        co2,
        scale_highest,
        ceiling,
        method,
        out_q,
    ):
        sample_data = {}
        with tqdm(
            total=len(samples),
            bar_format="progress: {desc:<10}|{bar:50}|",
            ascii=" >=",
            position=0,
            leave=False,
        ) as pbar:
            for sample in samples:
                sample_data[sample] = self.random_walk(
                    T=T,
                    P=P,
                    probability_threshold=probability_threshold,
                    path_depth=path_depth,
                    max_compounds=max_compounds,
                    max_rank=max_rank,
                    co2=co2,
                    scale_highest=scale_highest,
                    ceiling=ceiling,
                    method=method,
                )
                pbar.update(1)

        out_q.put(sample_data)

    def sampling_multiprocessing(self, T=None, P=None, **kw):
        init_concs = copy.deepcopy(self.concs)
        result_dict = {
            0: {"data": init_concs, "equation_statistics": [], "path_length": None}
        }
        # start the queue
        out_queue = pmp.Queue()
        samples = list(range(1, self.sample_length + 1, 1))
        data_chunks = [
            samples[chunksize * i : chunksize * (i + 1)]
            for i in range(self.nprocs)
            for chunksize in [int(math.ceil(len(samples) / float(self.nprocs)))]
        ]

        jobs = []
        for i, chunk in enumerate(data_chunks):
            process = pmp.Process(
                target=self._queue_function,
                args=(
                    i,
                    chunk,
                    T,
                    P,
                    self.probability_threshold,
                    self.path_depth,
                    self.max_compounds,
                    self.max_rank,
                    self.co2,
                    self.scale_highest,
                    self.ceiling,
                    self.method,
                    out_queue,
                ),
            )
            jobs.append(process)
            process.start()

        for proc in jobs:
            result_dict.update(out_queue.get())

        for proc in jobs:
            proc.terminate()

        for proc in jobs:
            proc.join()

        out_queue.close()

        return result_dict

    def sampling_serial(self, T=None, P=None, **kw):
        init_concs = copy.deepcopy(self.concs)
        result_dict = {
            0: {"data": init_concs, "equation_statistics": [], "path_length": None}
        }
        with tqdm(
            total=self.sample_length,
            bar_format="progress: {desc:<10}|{bar:50}|",
            ascii=" >=",
            position=0,
            leave=False,
        ) as pbar:
            for sample in range(self.sample_length):
                result_dict[sample + 1] = self.random_walk(
                    T=T,
                    P=P,
                    probability_threshold=self.probability_threshold,
                    path_depth=self.path_depth,
                    max_compounds=self.max_compounds,
                    max_rank=self.max_rank,
                    co2=self.co2,
                    scale_highest=self.scale_highest,
                    ceiling=self.ceiling,
                    method=self.method,
                )
                pbar.update(1)
        return result_dict

    def run(
        self,
        trange,
        prange,
        ic=None,
        save=False,
        savename=None,
        ignore_warnings=True,
        logging=False,
        **kw,
    ):
        if ignore_warnings:
            warnings.filterwarnings("ignore")

        """
        kwargs = sample_length,probability_threshold,max_compounds,max_rank,path_depth,nprocs,random_path_depth,co2=False
        """
        from loguru import logger
        from io import StringIO

        # setup logger
        stream = StringIO()
        logger.remove()
        logger.add(stream, format="{message}")

        num = 1
        total = len(trange) * len(prange)

        needed_args = self.__dict__
        for i in needed_args:
            if i in kw:
                self.__dict__[i] = kw[i]

        logger.info(
            """\n                                             
                                            
    // | |     //   ) )  //   ) )  //   ) ) 
   //__| |    //___/ /  //        ((        
  / ___  |   / ___ (   //           \\      
 //    | |  //   | |  //              ) )   
//     | | //    | | ((____/ / ((___ / /    
version:1.2
{}
        ->sample_length = {}
        ->probability_threshold = {}
        ->max_compounds = {}
        ->max_rank = {}
        ->path_depth = {}
        ->co2 = {}
        ->shortest path method = {}
        ->number of processes = {}
        ->concentration ceiling = {} %
        ->scale highest = {}
        ->rank smaller reactions higher = {}\n""".format(
                str(datetime.now()),
                self.sample_length,
                self.probability_threshold,
                self.max_compounds,
                self.max_rank,
                self.path_depth,
                self.co2,
                self.method,
                self.nprocs,
                self.ceiling,
                self.scale_highest,
                self.rank_small_reactions_higher,
            )
        )

        logger.info("initial concentrations (ppm):\n")
        self.concs = ic
        concstring = pd.Series({k: v for k, v in self.concs.items() if v > 0}) / 1e-6
        del concstring["CO2"]
        logger.info(concstring.to_string() + "\n")

        if logging:
            print(stream.get_value())

        path_lengths = []
        total_data = {}
        for T in trange:
            data_2 = {}
            final_concs_2 = {}
            initfinaldiff = {}
            for P in prange:
                start = datetime.now()
                logger.info(
                    "\n {}/{}: temperature = {}K, pressure = {}bar ".format(
                        num, total, T, P
                    ),
                    end="\n",
                )
                if self.nprocs > 1:
                    data_2[P] = self.sampling_multiprocessing(T, P, **kw)
                else:
                    data_2[P] = self.sampling_serial(T, P, **kw)

                finish = datetime.now() - start
                logger.info(
                    "-> completed in {} seconds".format(finish.total_seconds()),
                    end="\n",
                )
                reformatted = [
                    {x: v for x, v in data_2[P][i]["data"].items()} for i in data_2[P]
                ]
                mean = (
                    pd.Series(
                        {
                            k: v
                            for k, v in pd.DataFrame(reformatted).mean().items()
                            if v > 0.5e-6
                        }
                    ).drop("CO2")
                    / 1e-6
                )
                # mean = pd.Series({x:v for x,v in np.mean(pd.DataFrame(data_2[P][i]['data'] for i in data_2[P]).keys()) if v > 0.5e-6}).drop('CO2')/1e-6
                logger.info("\n final concentrations (>0.5ppm):\n")
                logger.info(mean.round(1).to_string())
                final_concs_2[P] = mean.to_dict()
                diff_concs = pd.Series(mean.to_dict()) - pd.Series(
                    {k: v / 1e-6 for k, v in self.concs.items()}
                )
                ift = pd.DataFrame(
                    [
                        {k: v / 1e-6 for k, v in self.concs.items() if v > 0},
                        mean.to_dict(),
                        diff_concs.to_dict(),
                    ],
                    index=["initial", "final", "change"],
                ).T
                initfinaldiff[P] = ift.dropna(how="all").fillna(0.0).to_dict()
                avgpathlength = np.median(
                    [
                        data_2[P][i]["path_length"]
                        for i in data_2[P]
                        if data_2[P][i]["path_length"] is not None
                    ]
                )

                logger.info("\n median path length: {}".format(avgpathlength))
                path_lengths.append(avgpathlength)
                num += 1
            total_data[T] = data_2
            self.final_concs[T] = final_concs_2
            self.initfinaldiff[T] = initfinaldiff

        if save:
            from monty.serialization import dumpfn

            if not savename:
                from datetime import date

                today = str(date.today())
                savename = "sampling_{}.json".format(today)
            dumpfn(total_data, savename, indent=4)

        self.metadata = {
            "arcs_version": "1.4.0",
            "avg_path_length": np.mean(path_lengths),
            "co2": self.co2,
            "max_compounds": self.max_compounds,
            "probability_threshold": self.probability_threshold,
            "shortest_path_method": self.method,
            "max_rank": self.max_rank,
            "sample_length": self.sample_length,
            "path_depth": self.path_depth,
            "random_path_depth": self.random_path_depth,
            "nprocs": self.nprocs,
            "ceiling": self.ceiling,
            "scale_highest": self.scale_highest,
            "rank_small_reactions_higher": self.rank_small_reactions_higher,
            "platform": platform.platform(),
            "python_version": platform.python_version(),
            "processor": platform.processor(),
            "available_cores": psutil.cpu_count(),
            "available_memory": str(
                int(psutil.virtual_memory()[0] / 1000 / 1000 / 1000)
            )
            + "Gb",
            "date": str(datetime.now()),
        }

        self.data = total_data
        if logging:
            print(stream.get_value())


# done
