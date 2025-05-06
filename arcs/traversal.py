from chempy.equilibria import Equilibrium,EqSystem
from chempy import Substance
import copy
import networkx as nx
import warnings
from datetime import datetime
import numpy as np
import pathos.multiprocessing as pmp 
import psutil
import time
import datetime
import os

class Traversal:

    def __init__(self,graph):
        self.graph = graph

        #default values:
        self.exclude_co2 = False
        self.max_compounds = 5
        self.discovery_threshold=5 # % percent
        self.maximum_reaction_numbe = 5
        self.nsamples=1000
        self.max_steps=10
        self.ncpus = 4    
        self.ceiling = 2000
        self.scale_largest=10
        self.rank_small_reactions_higher=True
        self.shortest_path_method='Djikstra'

    def length_multiplier(
            self,
            candidate_reaction:int,
            by_coefficients:bool = False
            ):
        """
        given a candidate reaction, if self.rank_small_reactions_higher == True, then return the length of the reaction as a multiplier
        i.e. H2 + 1/2 O2 = H2O has a length multiplier of 3  
        """
        if self.rank_small_reactions_higher:
            reaction_dict = self.graph.nodes[candidate_reaction]['reaction']
            if by_coefficients:
                num_reactants = np.sum(list(reaction_dict['reactants'].values()))
                num_products = np.sum(list(reaction_dict['products'].values()))
            else:
                num_reactants = len(reaction_dict['reactants'])
                num_products = len(reaction_dict['products'])                
            return(num_reactants + num_products)
              # should probably include coefficients
        else:
            return (1)

    @staticmethod
    def scale_large_concentrations(
        concentrations:dict,
        scale_largest:float,
        ceiling:float,
        )->dict:
        """
        function that takes a dict of concentrations and scales abnormally large concentrations (above ceiling DEFAULT = 3000%) and scales them by scale_highest (DEFAULT = 10% of original value)

        this is to be used with self.get_weighted_random_compounds
        """

        median_conc = np.median([v for v in concentrations.values() if v > 0])
        species_above_ceiling = {k:v for k,v in concentrations.items() if v > (median_conc * (1+(ceiling/100)))}
        #modify the ceiling by scaling it down to a suitable value 
        #should still max out if concentrations become way to high 
        for k,v in species_above_ceiling.items():
            concentrations[k] = v*1/scale_largest

        return(concentrations)

    def get_weighted_random_compounds(
            self,
            concentrations: dict,
            exclude_co2: bool = True,
            max_compounds: int = 5,
            discovery_threshold: float = 5,  # discovery threshold in %
            scale_largest: float = 10,  # how much to scale the highest components in %
            ceiling: float = 1000,  # ceiling percent larger than the median average in %
    ) -> dict:
        """
        given a dictionary of concentrations e.g. {'H2O':100,'NO2':50} a weighted ranking can be returned with probabilities given a discovery threshold DEFAULT = 5%. 

        exceedingly large concentrations (up to ceiling % DEFAULT = 1000% above the median concentration) that may occur are scaled using self.scale_large_concentrations (scaled with scale_largest DEFAULT = 10% of original value) such that reactions may continue even with very large concentrations of species up to a point. 

        returns a dictionary with length up to max_compounds depending on the discovery_threshold and scale_largest factors. 

        CO2 is by default excluded (exclude_co2 = True) as it is considered background, however this can be turned on if you want to test CO2 containing reactions.
        """

        concs = copy.deepcopy(concentrations)

        if exclude_co2:
            # CO2 will always be too large as it is the background
            del concs['CO2']

        # scale potential large concentrations
        concs = self.scale_large_concentrations(
            concentrations=concs, scale_largest=scale_largest, ceiling=ceiling)

        # get the probabilities based upon relative concentrations:
        p_1 = {k: v/sum(concs.values()) for k, v in concs.items()}
        # now filter based upon the probability threshold: (discovery)
        p_2 = {k: v for k, v in p_1.items() if v > discovery_threshold/100}
        # remake the probabilities
        p_3 = {k: v/sum(p_2.values()) for k, v in p_2.items()}
        # make a list of choices based upon the probabilities
        available = list(
            np.random.choice(
                a=list(p_3), size=len(concs)*10, p=list(p_3.values())
            )
        )
        #  now make a list max_compounds long of random choices based on available
        choices = {} 
        for i in range(max_compounds):
            try:
                compound = np.random.choice(available)
                choices[compound] = p_3[compound]

                available = list(
                    filter(
                        lambda a: a != list(choices)[i-1], available
                    )
                )
            except ValueError:
                pass

        return(choices)

    def get_weighted_reaction_rankings(
            self,
            weighted_random_compounds:dict,
            maximum_reaction_number:int = 20,
            shortest_path_method:str = 'Djikstra',
            ):
        """
        given a dictionary of weighted random compounds from self.get_weighted_random_compounds find the shortest path between: 
        1. compound 0 and compound 1 -> list
        2. filter based on available compounds 
        3. weight based on edge weight and (optional) length multiplier for unreasonable large reactions (if option selected)
        returns a dictionary of ranked reactions and their weighting. 
        """
        
        if len(weighted_random_compounds) == 1:
            return(None)

        rankings = {}

        shortest_paths = nx.shortest_paths.all_shortest_paths(
                G=self.graph,
                source=list(weighted_random_compounds)[0],
                target=list(weighted_random_compounds)[1],
                method=shortest_path_method)  # gets a list of shortest paths without weights first
                
        rankings = {}
        for path in shortest_paths:
            source,reaction,target = path
            reaction_compounds = list(self.graph[reaction])
            #find reactions with >3rd compound
            if len(weighted_random_compounds) > 2:
                for compound in list(weighted_random_compounds)[2:]:
                    if compound in reaction_compounds:
                        rankings[reaction] = self.graph.get_edge_data(
                            u=source,
                            v=reaction
                        )[0]['weight']*self.length_multiplier(reaction,by_coefficients=True) #need to play around with coefficients=True
            else:
                rankings[reaction] = self.graph.get_edge_data(
                    u=source,
                    v=reaction
                )[0]['weight']*self.length_multiplier(reaction,by_coefficients=True)
        #limit based on maximum_reaction_number:        
        rankings = {k:rankings[k] for k in list(rankings)[0:maximum_reaction_number]}
        return(rankings)
    
    @staticmethod
    def choose_reaction(ranked_reactions:dict)->int:
        """
        given a dictionary of ranked reactions from self.get_weighted_reaction_rankings
        chose a reaction based on weights and probabilities
        """
        weights = {k:1/v for k,v in ranked_reactions.items()} # here higher is better 
        probabilities = {k:v/sum(weights.values()) for k,v in weights.items()}
        chosen_reaction = np.random.choice(
                [
                    np.random.choice(
                        a=list(probabilities.keys()),
                        size=len(probabilities)*10,
                        p=list(probabilities.values())
                    )
                ][0]
            )
        return(chosen_reaction)
    
    def generate_chempy_eqsystem(
            self,
            index:int
    )->EqSystem:
        """
        given a reaction index form a chempy.equilibria.EqSystem 

        eventually this will be deprecated as it is a speed bottleneck

        needs to involve charged species

        """
        node_dict = self.graph.nodes[index]
        reactants = node_dict['reaction']['reactants']
        products = node_dict['reaction']['products']
        k = node_dict['equilibrium_constant']

#        substances = {}
#        for compound in list(
#            it.chain(
#                *[list(reactants)+list(products)]
#            )
#        ):
#            substances[compound] = Substance.from_formula(compound,**{'charge':0}) # use if substance_factory doesnt work

        equation = Equilibrium(reac=reactants,prod=products,param=k)
        try:
            return(EqSystem([equation],substance_factory=Substance.from_formula)) # might not just be able to try a return...
        except Exception:
            return(None)
        
    def chempy_equilibrium_concentrations(
            self,
            concentrations:dict,
            equilibrium_reaction:EqSystem
            )->dict:
        """
        generate equilibrium concentrations 
        if the reaction is "sane" and a "success" then it returns the equilibrium concentrations as a dict of concentrations 
        elsewise 
        return None
        """
        warnings.simplefilter('ignore')
        _concs = copy.deepcopy(concentrations)
        try:
            result = equilibrium_reaction.solve(init_concs=_concs)
            assert result.success and result.sane
            for compound,concentration in enumerate(result.conc):
                _concs[equilibrium_reaction.substance_names()[compound]] = concentration
            return(_concs)
        except Exception:
            return(None)
    
    def random_walk(
            self,
            initial_concentrations: dict,
            max_steps: int = 10,
            discovery_threshold: float = 5,
            max_compounds: int = 5,
            exclude_co2: bool = False,
            scale_largest: float = 10,  # in %
            ceiling: float = 1000,  # in %
            maximum_reaction_number: int = 5,
            shortest_path_method='djikstra'
    ) -> dict:
        """
        does a random sampling of the reaction network with max_steps  DEFAULT = 10.

        for each sample step:
        1. get weighted random compounds 
        2. get ranked reactions 
        3. choose a reaction
        3. generate a chempy eqsystem
        4. calculate the equilibrium concentrations 
        5. update the concentrations and reaction statistics 
        """

        concentrations = {0: initial_concentrations}
        reactionstats = {0: None}
        i = 0
        for step in range(1, max_steps+1):
            _concentrations = copy.deepcopy(concentrations[i])
            # 1 grab weighted_random_compounds
            weighted_random_compounds = self.get_weighted_random_compounds(
                concentrations=_concentrations,
                exclude_co2=exclude_co2,
                max_compounds=max_compounds,
                discovery_threshold=discovery_threshold,
                scale_largest=scale_largest,
                ceiling=ceiling
            )
            # 2 grab reaction rankings
            ranked_reactions = self.get_weighted_reaction_rankings(
                weighted_random_compounds=weighted_random_compounds,
                maximum_reaction_number=maximum_reaction_number,
                shortest_path_method=shortest_path_method,
            )
            # 3 if no reactions found then break the for loop
            if not ranked_reactions:
                break

            # 4 choose a reaction
            chosen_reaction_index = self.choose_reaction(
                ranked_reactions=ranked_reactions)
            # 5 generate a chgempy eqsystem
            eqsystem = self.generate_chempy_eqsystem(
                index=chosen_reaction_index)

            # 6 check that reaction wasn't previously done (can probably leave this out) (supposed to stop infinite loops, but that might be ok)
            # previous_reactions = [r for r in reactionstats.values() if r]
            # try:
            #    if eqsyst.string() == previous_reactions[-1]:
            #        break
            # except IndexError:
            #    pass

            # 6 get equilibrium_concentrations and update relevant dictionaries.
            final_concentrations = self.chempy_equilibrium_concentrations(
                concentrations=_concentrations,
                equilibrium_reaction=eqsystem
            )
            if final_concentrations:
                i += 1
                concentrations[i] = final_concentrations
                reactionstats[i] = self.graph.nodes[chosen_reaction_index]['reaction']['reaction_string']

        return (
            {'concentrations': concentrations,
                'equation_statistics': reactionstats
             }
        )
    
    def sampling_function(self,iterable):
        """
        sampling function to be multiprocessed - runs one random walk
        """
        initial_concentrations = self.initial_concentrations
        return(self.random_walk(initial_concentrations=initial_concentrations)) 
    

    def sample(
            self,
            initial_concentrations:dict,
            nsamples:int = 1000,
            ncpus:int = 4,
            **random_walk_kws:dict,
    )->dict:
        """
        samples the graph network nsamples DEFAULT = 1000
        multiprocessed with ncpus DEFAULT = 4  
        """
        self.initial_concentrations = initial_concentrations

        with pmp.Pool(processes=ncpus) as pool:
            data = pool.map(
                self.sampling_function, list(
                    range(
                        nsamples
                    )
                )
            )
        return(data)
