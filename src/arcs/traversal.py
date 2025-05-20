from chempy.equilibria import Equilibrium,EqSystem
from chempy import Substance
import copy
import networkx as nx
import warnings
import numpy as np
import tqdm_pathos
import itertools as it 
from arcs.generate import parse_molecule

class Traversal:

    def __init__(self,graph,max_reaction_length=5,**kws):
        self.graph = graph
        
        #default values:
        self.exclude_co2 = True
        self.max_compounds = 5
        self.discovery_threshold=5 # % percent
        self.maximum_reaction_number = 10
        self.max_steps=5
        self.ncpus = 4    
        self.ceiling = 2000
        self.scale_largest=10
        self.rank_small_reactions_higher = True
        self.rank_by_number_of_atoms = True
        self.shortest_path_method='Djikstra'
        self.__dict__.update(kws)

    def length_multiplier(
            self,
            candidate_reaction:int=None,
            **kws
            )->float:
        """
        given a candidate reaction, 
        
        if self.rank_small_reactions_higher == True, then return the sum of the coefficients as a multiplier

        i.e. H2 + 1/2 O2 = H2O has a length multiplier of 2.5 

        if self.rank_by_number_of_atoms = True, then rank by the number of atoms in the reactants (which = num_atoms_in_products)

        i.e. H2 + 1/2 O2 = H2O has a length multiplier of 3

        """
        self.__dict__.update(kws)

        if self.rank_small_reactions_higher:
            if self.rank_by_number_of_atoms:
                reaction_dict = self.graph.nodes[candidate_reaction]['reaction']
                reactants = reaction_dict['reactants']
                num_atoms = []
                for species,coefficient in reactants.items():
                    num_atoms.append(
                        np.sum(
                            list(parse_molecule(species).values())
                        ) * coefficient
                    )
                return(np.sum(num_atoms))
            else:
                reaction_dict = self.graph.nodes[candidate_reaction]['reaction']
                num_reactants = np.sum(list(reaction_dict['reactants'].values()))
                num_products = np.sum(list(reaction_dict['products'].values()))
                return((num_reactants + num_products))
        else:
            return (1)
        
    def check_reactant_atoms(
            self,
            reaction_index:int,
            weighted_random_compounds:list,
            **kws
            )->bool:
        """
        takes a reaction index, chosen random_compounds, and checks for atom balance
        returns bool
        """
        self.__dict__.update(kws)

        reac_atoms = list(
            dict.fromkeys(
                list(
                    it.chain(
                        *[
                            list(parse_molecule(x)) for x in self.graph.nodes()[reaction_index]['reaction']['reactants']
                        ]
                    )
                )
            )
        )
        random_compounds_atoms = list(
            dict.fromkeys(
                list(
                    it.chain(
                        *[
                            list(parse_molecule(x)) for x in weighted_random_compounds
                        ]
                    )
                )
            )
        )
        if sorted(reac_atoms) == sorted(random_compounds_atoms):
             return(True)
        else:
             return(False)

    def scale_large_concentrations(self, concentrations: dict, **kws) -> dict:
        """
        function that takes a dict of concentrations and scales abnormally large concentrations (above ceiling DEFAULT = 3000%) and scales them by scale_highest (DEFAULT = 10% of original value)

        this is to be used with self.get_weighted_random_compounds
        """
        self.__dict__.update(kws)

        median_conc = np.median([v for v in concentrations.values() if v > 0])
        species_above_ceiling = {k:v for k,v in concentrations.items() if v > (median_conc * (1+(self.ceiling/100)))}
        #modify the ceiling by scaling it down to a suitable value 
        #should still max out if concentrations become way to high 
        for k,v in species_above_ceiling.items():
            concentrations[k] = v*1/self.scale_largest

        return(concentrations)

    def get_weighted_random_compounds(
            self,
            concentrations: dict,
            **kws
    ) -> list:
        """
        given a dictionary of concentrations e.g. {'H2O':100,'NO2':50} a weighted ranking can be returned with probabilities given a discovery threshold DEFAULT = 5%. 

        exceedingly large concentrations (up to ceiling % DEFAULT = 1000% above the median concentration) that may occur are scaled using self.scale_large_concentrations (scaled with scale_largest DEFAULT = 10% of original value) such that reactions may continue even with very large concentrations of species up to a point. 

        returns a list with length up to max_compounds depending on the discovery_threshold and scale_largest factors. 

        CO2 is by default excluded (exclude_co2 = True) as it is considered background, however this can be turned on if you want to test CO2 containing reactions.
        """
        self.__dict__.update(kws)

        concs = copy.deepcopy(concentrations)

        if self.exclude_co2:
            # CO2 will always be too large as it is the background
            del concs['CO2']

        # scale potential large concentrations
        concs = self.scale_large_concentrations(
            concentrations=concs, scale_largest=self.scale_largest, ceiling=self.ceiling
        )

        # get the probabilities based upon relative concentrations:
        p_1 = {k: v/sum(concs.values()) for k, v in concs.items()}
        # now filter based upon the probability threshold: (discovery)
        p_2 = {k: v for k, v in p_1.items() if v >= self.discovery_threshold/100}
        if not p_2:
            return([])
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
        for i in range(self.max_compounds):
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
        return(list(choices)[0:np.random.randint(2,self.max_compounds)])
    
    def get_weighted_reaction_rankings(
            self,
            weighted_random_compounds: list,
            **kws
    ):
        """
        returns a dictionary of {<reaction_index>:<weighting>} given a list of weighted_random_compounds from self.get_weighted_random_compounds (needs at least 2 to give a result)

        algorithm follows:
        given weighted_random_compounds = ['NO2','H2O','O2']
        1. generates combinations of useable compounds 
            i.e. [['NO2','H2O'],['NO2','O2'],['H2O','O2']]
        2. for each combination, generate a list of shortest paths using networkx.shortest_paths.all_shortest_paths
        3. check that the reactant 

        """
        self.__dict__.update(kws)

        #return None if there isn't enough to make a reaction
        if len(weighted_random_compounds) <=1:
            return(None)
        # 1. generate possible combinations
        combinations = list(it.combinations(weighted_random_compounds, 2))
        # 2. generate shortest path possibilities from the combinations
        possibilities = []
        for compounds in combinations:
            possibilities.extend([
                x[1] for x in list(
                    nx.shortest_paths.all_shortest_paths(
                        G=self.graph, source=compounds[0], target=compounds[1],method=self.shortest_path_method)
                )
            ])

        # 3. check that all reactant atoms (=product atoms) are accounted for. 
        # this is so that self.get_chempy_equilibrium_concentrations gives a result.
        possibilities = [i for i in possibilities if self.check_reactant_atoms(
                reaction_index=i,
                weighted_random_compounds=weighted_random_compounds,
                **kws
                )]
        # 4. rank the possibilities based on edge weight and a length_multiplier
        # the length_multiplier is based on number of reaction coefficients
        # idea: perhaps by number of atoms ? 
        rankings = {}
        for i,reaction in enumerate(possibilities):
            for compound in weighted_random_compounds:
                try:
                    weight = self.graph.get_edge_data(
                        u=compound,
                        v=reaction
                        )[0]['weight']*self.length_multiplier(reaction,**kws)
                except TypeError:
                    pass
            rankings[reaction] = weight # this should be ammended to return None as well

        rankings = dict(sorted(rankings.items(),key=lambda item: item[1])[0:self.maximum_reaction_number])
        
        return(rankings)
    
    @staticmethod
    def choose_reaction(ranked_reactions:dict)->int:
        """
        given a dictionary of ranked reactions from self.get_weighted_reaction_rankings
        chose a reaction based on weights and probabilities
        """
        weights = {k:1/v**2 for k,v in ranked_reactions.items()} # here higher is better 
        #added a square multiplier to force more the larger coefficients
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

        substances = {}
        for compound in list(
            it.chain(
                *[list(reactants)+list(products)]
            )
        ):
            substances[compound] = Substance.from_formula(compound,**{'charge':0})  

        equation = Equilibrium(reac=reactants,prod=products,param=k)
        try:
            return(EqSystem([equation],substances=substances))
        except Exception:
            return(None)
    
    @staticmethod 
    def chempy_equilibrium_concentrations(
            concentrations:dict,
            equilibrium_reaction:EqSystem,
            chempy_sane=False,
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
            #assert result.success and result.sane
            if chempy_sane:
                assert result.success and result.sane
            else:
                assert result.success
            for compound,concentration in enumerate(result.conc):
                _concs[equilibrium_reaction.substance_names()[compound]] = concentration
            return(_concs)
        except Exception:
            return(None)
    
    def random_walk(
            self,
            initial_concentrations: dict,
            chempy_sane=True, # typically for very large 
            **kws
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
        self.__dict__.update(kws)

        concentrations = {0: initial_concentrations}
        reactionstats = {0: None}
        i = 0
        for step in range(1, self.max_steps+1):
            _concentrations = copy.deepcopy(concentrations[i])
            # 1 grab weighted_random_compounds
            weighted_random_compounds = self.get_weighted_random_compounds(
                concentrations=_concentrations,
                exclude_co2=self.exclude_co2,
                max_compounds=self.max_compounds,
                discovery_threshold=self.discovery_threshold,
                scale_largest=self.scale_largest,
                ceiling=self.ceiling,
                **kws
            )
            # 2 grab reaction rankings
            ranked_reactions = self.get_weighted_reaction_rankings(
                weighted_random_compounds=weighted_random_compounds,
                **kws
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

            # 6 get equilibrium_concentrations and update relevant dictionaries.
            final_concentrations = self.chempy_equilibrium_concentrations(
                concentrations=_concentrations,
                equilibrium_reaction=eqsystem,
                chempy_sane=chempy_sane
            )
            if final_concentrations:
                i += 1
                concentrations[i] = final_concentrations
                reactionstats[i] = {
                    'reaction': self.graph.nodes[chosen_reaction_index]['reaction'],
                    'equilibrium_constant': self.graph.nodes[chosen_reaction_index]['equilibrium_constant']
                }

        return (
            {'concentrations': concentrations,
            'reaction_statistics': reactionstats
             }
        )
    
    def sampling_function(self,iterable,**kws):
        """
        sampling function to be multiprocessed - runs one random walk
        """
        self.__dict__.update(kws)
        initial_concentrations = self.initial_concentrations
        return(
            self.random_walk(initial_concentrations=initial_concentrations,**kws)
            ) 
    

    def sample(
            self,
            initial_concentrations:dict,
            nsamples:int,
            ncpus:int,
            tqdm_kws:dict = {},
            **kws
    )->dict:
        """
        samples the graph network nsamples DEFAULT = 1000
        multiprocessed with ncpus DEFAULT = 4  
        """
        self.__dict__.update(kws)

        self.initial_concentrations = initial_concentrations
        
        data = tqdm_pathos.map(
                self.sampling_function,
                list(range(nsamples)),
                n_cpus=ncpus,
                tqdm_kwargs = tqdm_kws,
                **kws
            )

        return(data)
