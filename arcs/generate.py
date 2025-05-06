import os
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from ase.thermochemistry import IdealGasThermo
from scipy.constants import Boltzmann, e
import numpy as np 
import networkx as nx
from ase.atoms import Atoms
import ase
import re 
from collections import defaultdict
from monty.serialization import loadfn 

def get_compound_directory(base,compound,size):
    return(os.path.join(base,compound,size))

class GetEnergyandVibrationsVASP:
    """Class to get the Total Energy and Vibrations from a directory containing a calculations
    """
    def __init__(self,relax_directory,vibrations_directory):
        self.relax = relax_directory
        self.vibrations = vibrations_directory
    
    @staticmethod
    def get_initial_magnetic_moments(aseatoms:ase.Atoms)->list:
        magmoms = []
        for atomic_number in aseatoms.get_atomic_numbers():
            magmoms.append([0 if atomic_number %2 == 0 else 1][0])
        return(magmoms)
        
    def get_atoms(self)->ase.Atoms:
        """
        generates an ASE.Atoms object from a given POSCAR with magnetic moments
        """

        aseatoms = read('{}/POSCAR'.format(self.relax))

        aseatoms.set_initial_magnetic_moments(
            self.get_initial_magnetic_moments(aseatoms=aseatoms)
        )

        return(aseatoms)
        
    def get_energy(self)->float:
        """
        grabs the total energy from a VASP OUTCAR file (assumes one formula unit per cell)
        """
        outcar = open('{}/OUTCAR'.format(self.relax),'r').readlines()

        for line in outcar:
            if 'y=' in line:
                energy = float(line.split()[-1])
        if len(
            list(
                dict.fromkeys(self.get_atoms().get_atomic_numbers())
            )
        ) == 1:
            if not any(
                x in self.get_atoms().symbols.get_chemical_formula() for x in ['H2', 'O2', 'N2','S8'] # elemental species...
            ):
                energy = energy / self.atoms().get_global_number_of_atoms()
    
        return(energy)
    
    def get_spin(self)->int:
        """
        determines the spin of a system 
        """
        if self.get_atoms().get_chemical_formula() in ['O2','CO3']:
            return(1)
        else:    
            outcar = open('{}/OUTCAR'.format(self.relax),'r').readlines()
            for line in outcar:
                if 'NELECT' in line:
                    nelect = float(line.split()[2])
                    return([0 if nelect %2 == 0 else 0.5][0])
            
    def get_pointgroup(self)->str:
        """
        uses pymatgen's pymatgen.symmetry.analyzer.PointGroupAnalyzer class to determine the point group of an ase.Atoms object
        returns a string
        """
        aseatoms = self.get_atoms()
        pg = PointGroupAnalyzer(
            AseAtomsAdaptor.get_molecule(aseatoms)
        ).get_pointgroup()

        return(pg.sch_symbol)
    
    def islinear(self)->str:
        """
        determines whether the molecule is linear - can determine this from the point group (if * the molecule is linear)
        """
        num_at = self.get_atoms()
        if num_at.get_global_number_of_atoms() == 1:
            return('monatomic')
        else:
            pg = self.get_pointgroup()
            if '*' in pg:
                return('linear')
            else:
                return('nonlinear')

    def get_rotation_num(self)->int:
        """
        determines the rotational number from the point group
        returns an integer
        """
        pg = [x for x in self.get_pointgroup()]
        if pg[0] == 'C':
            if pg[-1] == 'i':
                rot = 1
            elif pg[-1] == 's':
                rot = 1
            elif pg[-1] == 'h':
                rot = int(pg[1])
            elif pg[-1] == 'v':
                if pg[1] == '*':
                    rot = 1
                else:
                    rot = int(pg[1])
            elif len(pg) == 2:
                rot = int(pg[-1])
                
        elif pg[0] == 'D':
            if pg[-1] == 'h':
                if pg[1] == '*':
                    rot = 2
                else:
                    rot = 2*int(pg[1])
            elif pg[-1] == 'd':
                rot = 2*int(pg[1])
            elif len(pg) == 2:
                rot = 2*int(pg[1])
            
        elif pg[0] == 'T':
            rot = 12
        
        elif pg[0] == 'O':
            rot = 24
        
        elif pg[0] == 'I':
            rot = 60        
                    
        return(rot)          
        
    def get_vibrations(self)->list:
        """
        gets the vibrational frequencies from a DFPT calculation as outputted in a VASP OUTCAR
        returns as list of both imaginary and normal modes of vibration in eV
        """
        outcar = open('{}/OUTCAR'.format(self.vibrations),'r').readlines()
        frequencies = []
        #i_frequencies = []
        for line in outcar:
            if 'THz' in line: 
                if 'f/i' not in line: 
                    frequencies.append(float(line.split()[-2])/1000)

                else: # imaginary frequencies are returned separately as a real number
                    frequencies.append(-float(line.split()[-2])/1000) # need to go back over this
        return(frequencies)
    
    def as_dict(self):
        return({'atoms':self.get_atoms().todict(),
                'pointgroup':self.get_pointgroup(),
                'spin':self.get_spin(),
                'rotation_num':self.get_rotation_num(),
                'islinear':self.islinear(),
                'energy':self.get_energy(),
                'vibrations':self.get_vibrations()})

class ReactionGibbsandEquilibrium: 
    """
    class that takes the relevant first principles data generated by GetEnergyandVibrationsVASP and converts it into Gibbs Free Energies and equilibrium constants for a given reaction 
    """
    
    def __init__(
        self,
            reaction_input: dict
    ):

        self.reaction_input = reaction_input
    
    @staticmethod
    def Gibbs(
            dft_dict:dict,
            temperature: float, # in K
            pressure: float # in bar
    ) -> float:
        """
        uses ASE's IdealGasThermo to calculate the relevant Gibbs Free Energy of a given compound at a specific temperature (K) and pressure (bar converted to Pa automatically)
        """
        igt = IdealGasThermo(
            vib_energies=dft_dict['vibrations'],
            geometry=dft_dict['islinear'],
            potentialenergy=dft_dict['energy'],
            atoms=Atoms.fromdict(dft_dict['atoms']),
            symmetrynumber=dft_dict['rotation_num'],
            spin=dft_dict['spin'],
            natoms=Atoms.fromdict(dft_dict['atoms']).get_global_number_of_atoms(),
            ignore_imag_modes=True
        )
        gibbs_free_energy = igt.get_gibbs_energy(temperature,pressure*100000,verbose=False)

        return(gibbs_free_energy)
    
    def reaction_gibbs(
            self,
            reaction: dict,#chempy.Equilibrium,
            pressure: float,  # in in bar
            temperature: float,  #  in K
    ) -> float:
        """
        returns the Gibbs Free Energy of Reaction for a given reaction (dict).
        example reaction:
        {'reaction_string': '1 H2O + 1 H2CO = 1 H2 + 1 CH2O2',
        'reactants': {'H2O': 1, 'H2CO': 1},
        'products': {'H2': 1, 'CH2O2': 1}},
        """
        products = reaction['products']
        reactants = reaction['reactants']
        reaction_compounds = list(products)+list(reactants)

        gibbs = {
            compound:self.Gibbs(
                dft_dict=self.reaction_input[compound],
                temperature=temperature,
                pressure=pressure,
                ) for compound in reaction_compounds
                }
        prod_sum = np.sum(
            [
                gibbs[compound]*products[compound] for compound in gibbs if compound in products
                ]
            )
        reac_sum = np.sum(
            [
                gibbs[compound]*reactants[compound] for compound in gibbs if compound in reactants
                ]
            )

        return(float(prod_sum - reac_sum))
    
    @staticmethod
    def equilibrium_constant(
        gibbs_free_energy: float,
        temperature=float,  # in K
    )->float:
        """
        returns an equilibrium constant, K, from the Gibbs Free Energy of Reaction (in eV)
        from DeltaG = -k_B T ln(K)
        """
        K = np.exp(
                -(gibbs_free_energy*e)/(Boltzmann*temperature)
                )
        return(K)
    
    def get_reaction_gibbs_and_equilibrium(self,
                                           reaction,
                                           temperature=float,# in K
                                           pressure=float,#in bar
                                           )->dict:
        """
        given a reaction dictionary object (from reactit), the function generates a dictionary with both gibbs free energy and equilibrium constant

        example reaction dictionary: 

        {'reaction_string': '1 H2O + 1 H2CO = 1 H2 + 1 CH2O2',
        'reactants': {'H2O': 1, 'H2CO': 1},
        'products': {'H2': 1, 'CH2O2': 1}},

        """
        if isinstance(reaction,str):
            from reactit import ReactionGenerator
            try:
                reaction_dict = ReactionGenerator.get_reactants_products(reaction)
                reaction = {'reaction_string':reaction,'reactants':reaction_dict[0],'products':reaction_dict[1]}
            except Exception as e:
                print(e)

        gibbs_free_energy = self.reaction_gibbs(reaction=reaction,temperature=temperature,pressure=pressure)
        equilibrium_constant = self.equilibrium_constant(gibbs_free_energy=gibbs_free_energy,temperature=temperature)
        return (
            {'g': gibbs_free_energy,
             'k': equilibrium_constant}
        )

class GraphGenerator:    
    """
    generates an nx.multidigraph object with weightings from the dft data gibbs free energy and equilibrium constants.
    """
    
    def __init__(self):
        """
        GraphGenerator
        """
        #self.applied_reactions = applied_reactions

    @staticmethod
    def parse_molecule(formula:str)->dict: 
        """
        parses a molecule string i.e. 'H2O' into a dictionary broken down into elemental counts 
        i.e. {'H':2,'O':1} 
        taken with permission from https://github.com/badw/reactit.git
        """
        # Regular expression to match elements and their counts 
        pattern = r'([A-Z][a-z]?)(\d*)' 
        matches = re.findall(pattern, formula) 
         
        atom_count = defaultdict(int) 
     
        for element, count in matches: 
            if count == '': 
                count = 1  # Default count is 1 if not specified 
            else: 
                count = int(count)  # Convert count to integer 
             
            atom_count[element] += count 
     
        return dict(atom_count)      

    def cost_function(
            self,
        gibbs_free_energy: float,  # in eV
        temperature: float,  # in K
        reactants: dict
    ) -> float:
        """
        takes a gibbs free energy of reaction and normalises it using a cost function taken from:
        https://www.nature.com/articles/s41467-021-23339-x.pdf which takes the form:
        cost = ln(1+(273/temperature)*exp(G/num_reactant_atoms))
        
        (this is normalised per reactant atom)

        i.e. CO2 + H2 = H2O + CO 

        num_reactant_atoms = 5

        """

        compounds = []
        for reactant,coefficient in reactants.items():
            for i in range(coefficient):
                compounds.append(reactant)

        num_atoms = np.sum(
            [
                np.sum(
                    [
                        y for x, y in self.parse_molecule(compound).items() 
                    ]
                )
                for compound in compounds
            ]
        )

        return(
            np.log(1+(273/temperature)*np.exp(gibbs_free_energy/num_atoms/1))
            )

    def generate_multidigraph(
            self,
            temperature:float, #in K
            applied_reactions:list,
            )-> nx.multidigraph:
        """
        This function generates reaction graph in networkx weighted using the self.costfunction 
        returns an nx.multidigraph object
        """

        graph = nx.MultiDiGraph(directed=True)

        for i,reaction in enumerate(applied_reactions):
            forward_cost = self.cost_function(
                gibbs_free_energy=reaction['g'],
                temperature=temperature,
                reactants=reaction['r']['reactants']
            )
            backward_cost = self.cost_function(
                gibbs_free_energy=-reaction['g'],
                temperature=temperature,
                reactants=reaction['r']['products']
                )

            graph.add_weighted_edges_from(
                [compound,i,forward_cost] for compound in reaction['r']['reactants']
                ) #reactants -> reaction
            graph.add_weighted_edges_from(
                [i,compound,backward_cost] for compound in reaction['r']['reactants']
                ) #reaction -> reactants
            graph.add_weighted_edges_from(
                [i,compound,forward_cost] for compound in reaction['r']['products']
                ) #reaction -> products
            graph.add_weighted_edges_from(
                [compound,i,backward_cost] for compound in reaction['r']['products']
                ) #products -> reaction

        return(graph)
    

    def from_file(
            self,
            filename:str,
            temperature:float, # in K
            pressure:float # in bar
            )->nx.MultiDiGraph:
        """
        generates a networkx.multidigraph from a .json file of reactions and dft_dict generated with reactit and GetEnergyandVibrationsVASP
        needs a sanity check on the file though

        file takes the format of a dictionary with: 
        {
    '<compound>': {
        'atoms': ase.Atoms.asdict(),
        'pointgroup': <pointgroup>,
        'spin': <spin>,
        'rotation_num': <rotation_num>,
        'islinear': <islinear>,
        'energy': <energy>,
        'vibrations': <vibrations
        }
    'reactions':list({
        'reaction_string':<reaction_string>,
                      'reactants':dict(reactants),
                      'products':dict(products)}
                      )
        }

        reactions can be generated with https://github.com/badw/reactit.git 
        """
        dft_dict = loadfn(filename)
        rge = ReactionGibbsandEquilibrium(dft_dict)
        applied_reactions = [] 
        for reaction in dft_dict['reactions'].values():
            g_k_dict = rge.get_reaction_gibbs_and_equilibrium(reaction=reaction,temperature=temperature,pressure=pressure)
            g_k_dict['r'] = reaction
            applied_reactions.append(g_k_dict)  

        graph = self.generate_multidigraph(applied_reactions=applied_reactions,temperature=temperature)
        return(graph)
    
    def from_dict(
            self,
            dft_dict:str,
            temperature:float, # in K
            pressure:float # in bar
            )->nx.MultiDiGraph:
        """
        generates a networkx.multidigraph from a dict representation of the .json file of reactions and dft_dict generated with reactit and GetEnergyandVibrationsVASP
        needs a sanity check on the file though

        dictionary takes the format of:
                {
    '<compound>': {
        'atoms': ase.Atoms.asdict(),
        'pointgroup': <pointgroup>,
        'spin': <spin>,
        'rotation_num': <rotation_num>,
        'islinear': <islinear>,
        'energy': <energy>,
        'vibrations': <vibrations
        }
    'reactions':list({
        'reaction_string':<reaction_string>,
                      'reactants':dict(reactants),
                      'products':dict(products)}
                      )
        }

        reactions can be generated with https://github.com/badw/reactit.git 
        """
        rge = ReactionGibbsandEquilibrium(dft_dict)
        applied_reactions = [] 
        for reaction in dft_dict['reactions'].values():
            g_k_dict = rge.get_reaction_gibbs_and_equilibrium(reaction=reaction,temperature=temperature,pressure=pressure)
            g_k_dict['r'] = reaction
            applied_reactions.append(g_k_dict)  

        graph = self.generate_multidigraph(applied_reactions=applied_reactions,temperature=temperature)
        return(graph)

class GenerateInitialConcentrations:
    """
    a class for generating initial concentrations in a standardised way
    typically in ppm or 1e-6

    this should be used to update concentrations and log them as well in future
    """

    def __init__(self,graph:nx.MultiDiGraph):
        self.graph = graph

    def all_random(self,include_co2=True)->dict:
        """
        generate fully random concentrations
        returns a dictionary of compounds with a random concentration 
        """
        compounds = [node for node in self.graph.nodes() if isinstance(node,str)]
        ic = {c:np.random.random()/1e5 for c in compounds}
        if not include_co2:
            ic['CO2'] = 1
        self.ic = ic
        return(self.ic)
    
    def all_zero(self,include_co2=True)->dict:
        """
        makes a dictionary of blank concentrations
        returns a dictionary of compounds with 0 concentration 
        """
        compounds = [node for node in self.graph.nodes() if isinstance(node,str)]
        ic = {c:0 for c in compounds}
        if not include_co2:
            ic['CO2'] = 1
        self.ic = ic
        return(self.ic)

    def specific_random(self,compounds=list,include_co2=True)->dict:
        """
        makes a dictionary of compounds with only certain compounds given random concentrations 
        i.e. 
        specific_random(compounds=['CO2','H2O','O2'])
        """
        full_list = [n for n in self.graph.nodes() if isinstance(n,str)]
        ic = {}
        for c in full_list:
            if c in compounds:
                ic[c] = np.random.random()/1e6
            else:
                ic[c] = 0 
        if not include_co2:
            ic['CO2'] = 1
        self.ic = ic
        return(self.ic)
    
    def update_ic(self,update_dict,include_co2=True)->dict:
        """
        update dict = {'CO2':1e-6,'H2O':300e-5} etc.
        used to update dictionary items
        need to have this log previous concentrations as well. 
        """
        try:
            hasattr(self.ic,'ic')
        except Exception:
            # if not initial concentrations created previously 
            self.all_zero(include_co2=include_co2)
        for k,v in update_dict.items():
            self.ic[k] = v
        return(self.ic)
