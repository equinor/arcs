import functools 
import os
import gzip
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from ase.thermochemistry import IdealGasThermo
from scipy.constants import Boltzmann, e
from monty.serialization import loadfn
import numpy as np 
from chempy.reactionsystem import Substance
from tqdm import tqdm
import networkx as nx
from pathos.helpers import mp as pmp
import math
import pickle
from ase.atoms import Atoms
from chempy import Equilibrium
import chempy 

def get_compound_directory(base,compound,size):
    return(os.path.join(base,compound,size))

class GetEnergyandVibrationsVASP:
    '''Class to get the Total Energy and Vibrations from a directory containing a calculations'''
    def __init__(self,relax_directory,vibrations_directory):
        self.relax = relax_directory
        self.vibrations = vibrations_directory
        
    def atoms(self):
        def _get_initial_magnetic_moments(ats):
            magmoms = []
            for atnum in ats.get_atomic_numbers():
                magmoms.append([0 if atnum %2 == 0 else 1][0])
            return(magmoms)
        structure = read('{}/POSCAR.gz'.format(self.relax))
        structure.set_initial_magnetic_moments(_get_initial_magnetic_moments(structure))
        return(structure)
        
    def energy(self):
        outcar = gzip.open('{}/OUTCAR.gz'.format(self.relax),'tr').readlines()
        for line in outcar:
            if 'y=' in line:
                energy = float(line.split()[-1])
        if len(list(dict.fromkeys(self.atoms().get_atomic_numbers()))) == 1:
            if not any(x in self.atoms().symbols.get_chemical_formula() for x in ['H2','O2','N2']):
                energy = energy / self.atoms().get_global_number_of_atoms()
    
        return(energy)
    
    def spin(self):
        if self.atoms().get_chemical_formula() in ['O2','CO3']:
            return(1)
        else:    
            outcar = gzip.open('{}/OUTCAR.gz'.format(self.relax),'tr').readlines()
            for line in outcar:
                if 'NELECT' in line:
                    nelect = float(line.split()[2])
                    return([0 if nelect %2 == 0 else 0.5][0])
            
    def pointgroup(self):
        atoms = self.atoms()
        pg = PointGroupAnalyzer(AseAtomsAdaptor.get_molecule(atoms)).get_pointgroup()
        return(pg.sch_symbol)
    
    def islinear(self):
        num_at = self.atoms()
        if num_at.get_global_number_of_atoms() == 1:
            return('monatomic')
        else:
            pg = self.pointgroup()
            if '*' in pg:
                return('linear')
            else:
                return('nonlinear')

    def rotation_num(self):
        pg = [x for x in self.pointgroup()]
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
        
    def get_vibrations(self):
        outcar = gzip.open('{}/OUTCAR.gz'.format(self.vibrations),'tr').readlines()
        frequencies = []
        for line in outcar:
            if 'THz' in line: #thermochemistry now has an option to ignore imaginary modes...
                if 'f/i' not in line: # we ignore imaginary modes
                    ev = float(line.split()[-2]) / 1000
                    frequencies.append(ev)
                else:
                    ev = -float(line.split()[-2]) /1000
                    frequencies.append(ev)
        return(frequencies)      
    
    def as_dict(self):
        return({'atoms':self.atoms().todict(),
                'pointgroup':self.pointgroup(),
                'spin':self.spin(),
                'rotation_num':self.rotation_num(),
                'islinear':self.islinear(),
                'energy':self.energy(),
                'vibrations':self.get_vibrations()})
    
    def get_vibrations(self):
        daltonout = open('{}/output.out'.format(self.vibrations),'r').readlines()
        datalines = []
        for i,line in enumerate(daltonout):
            if 'frequency' in line and 'mode' in line:
                datalines.append(i+4)
            elif 'Normal Coordinates' in line:
                datalines.append(i)
        recip_cm = [float(x.split()[2])/8100 for x in daltonout[datalines[0]:datalines[1]] if len(x.split()) > 1] # 8100 conversion cm-1 -> eV
        return(recip_cm) 

    def as_dict(self):
        return({'atoms':self.atoms(),
                'spin':self.spin(),
                'rotation_num':self.rotation_num(),
                'islinear':self.islinear(),
                'energy':self.energy(),
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
            reaction: chempy.Equilibrium,
            pressure: float,  # in in bar
            temperature: float,  #  in K
    ) -> float:
        """
        returns the Gibbs Free Energy of Reaction for a given reaction (chempy.Equilibrium object).
        currently does no checks on charge neutrality or mass balance. 
        """
        prod = reaction.prod
        reac = reaction.reac
        reaction_compounds = list(prod)+list(reac)
        gibbs = {
            compound:self.Gibbs(
                dft_dict=self.reaction_input[compound],
                temperature=temperature,
                pressure=pressure,
                ) for compound in reaction_compounds
                }
        prod_sum = np.sum(
            [gibbs[compound]*prod[compound] for compound in gibbs if compound in prod]
            )
        reac_sum = np.sum(
            [gibbs[compound]*reac[compound] for compound in gibbs if compound in reac]
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
                                           reaction:chempy.Equilibrium,
                                           temperature=float,# in K
                                           pressure=float,#in bar
                                           )->dict:
        """
        given a chempy.Equilibrium reaction object, the function generates a dictionary with both gibbs free energy and equilibrium constant
        """
        gibbs_free_energy = self.reaction_gibbs(reaction=reaction,temperature=temperature,pressure=pressure)
        equilibrium_constant = self.equilibrium_constant(gibbs_free_energy=gibbs_free_energy,temperature=temperature)
        return (
            {'g': gibbs_free_energy,
             'k': equilibrium_constant}
        )

class ApplyDataToReaction:
    ''' this class applies the Gibbs data to a specific reaction'''
    
    def __init__(self,trange,prange,data,nprocs):
        self.trange = trange
        self.prange = prange
        reactions = data['reactions']
        try:
            self.reactions = {i:Equilibrium.from_string(r) for i,r in enumerate(reactions)}
        except Exception:
            self.reactions = {i:r for i,r in enumerate(reactions)}
        self.compound_data = {k:data[k] for k in data.keys() if not k == 'reactions'}
        self.nprocs = nprocs
        self.barformat = '{desc:<20}{percentage:3.0f}%|{bar:10}{r_bar}'
        
#    def _generate_data_serial(self,t,p): #serial
#        reactions = {i:{'e':r,
#            'k':ReactionGibbsandEquilibrium(t,p,self.compound_data).#equilibrium_constant(r),
#            'g':ReactionGibbsandEquilibrium(t,p,self.compound_data).#reaction_energy(r)} 
#                     for i,r in self.reactions.items()}
#        return(reactions)
    
    def _generate_data_serial(self,t,p):
        rge = ReactionGibbsandEquilibrium(t,p,self.compound_data)
        reactions = {i:rge.as_dict(r) for i,r in self.reactions.items()}
        return(reactions)

    def generate_data(self,t,p): #multiprocessed

        manager = pmp.Manager()
        queue = manager.Queue()
        
        def mp_function(reaction_keys,out_q):

            data = {}
            for r in reaction_keys:
                rge = ReactionGibbsandEquilibrium(t,p,self.compound_data)
                data[r] = rge.as_dict(self.reactions[r])
            out_q.put(data)

        resultdict = {}
        r_keys = list(self.reactions.keys())
        chunksize = int(math.ceil(len(self.reactions)/float(self.nprocs)))
        processes = []

        for i in range(self.nprocs):
            pr = pmp.Process(target=mp_function,
                            args=(r_keys[chunksize*i:chunksize*(i+1)],queue))
            processes.append(pr)
            pr.start()

        for i in range(self.nprocs):
            resultdict.update(queue.get(timeout=1800))

        for pr in processes:
            pr.join()

        return(resultdict)


    def apply(self,serial=False):
        data = {}
        for t in self.trange:
            pdat = {}
            for p in self.prange:
                if serial:
                    pdat[p] = self._generate_data_serial(t,p)
                else:
                    pdat[p] = self.generate_data(t,p)
            data[t] = pdat
        self.data = data
        return(self.data) 
    
    def save(self,filename='applied_reactions.p'):
        pickle.dump(self.data,open(filename,'wb'))
        print('data saved to: {}'.format(filename))

class GraphGenerator:    
    
    def __init__(self,applied_reactions,ncores=4):
        try:
            self.applied_reactions = pickle.load(open(applied_reactions,'rb'))
        except Exception:
            self.applied_reactions = applied_reactions 
        self.trange = list(self.applied_reactions)
        self.prange = list(self.applied_reactions[self.trange[0]])
        self.ncores = ncores

    def _cost_function(self,gibbs,T,reactants):
        '''takes the cost function that is used in https://www.nature.com/articles/s41467-021-23339-x.pdf
        this is normalised per reactant atom'''

        comps = []
        for r,n in reactants.items():
            for i in range(n):
                comps.append(r)

        num_atoms = np.sum([np.sum([y 
                             for x,y in Substance.from_formula(c).composition.items()]) 
                     for c in comps]) 

        return(np.log(1+(273/T)*np.exp(gibbs/num_atoms/1)))

    def multidigraph_cost(self,T,P):
        ''' this will weight the graph in terms of a cost function which makes it better for a Djikstra algorithm to work'''
        t = nx.MultiDiGraph(directed=True)
        for i,reac in tqdm(self.applied_reactions[T][P].items()):
            f_cost = self._cost_function(reac['g'],T,reac['e'].reac) #forward cost
            b_cost = self._cost_function(-reac['g'],T,reac['e'].prod) #backward cost
            r = list(reac['e'].reac)
            p = list(reac['e'].prod)
            t.add_weighted_edges_from([c,i,f_cost] for c in r) #reactants -> reaction
            t.add_weighted_edges_from([i,c,b_cost] for c in r) #reaction -> reactants
            t.add_weighted_edges_from([i,c,f_cost] for c in p) #reaction -> products
            t.add_weighted_edges_from([c,i,b_cost] for c in p) #products -> reaction
        return(t)
    
    def multidigraph_cost_mp(self,T,P):

        def mp_func(rr,t):
            i,reac = rr
            f_cost = self._cost_function(reac['g'],T,reac['e'].reac) #forward cost
            b_cost = self._cost_function(-reac['g'],T,reac['e'].prod) #backward cost
            r = list(reac['e'].reac)
            p = list(reac['e'].prod)
            t.add_weighted_edges_from([c,i,f_cost] for c in r) #reactants -> reaction
            t.add_weighted_edges_from([i,c,b_cost] for c in r) #reaction -> reactants
            t.add_weighted_edges_from([i,c,f_cost] for c in p) #reaction -> products
            t.add_weighted_edges_from([c,i,b_cost] for c in p) #products -> reaction
        
        import tqdm_pathos as ptqdm

        t = nx.MultiDiGraph(directed=True)

        rr = list(self.applied_reactions[T][P].items())

        ptqdm.map(mp_func, rr,t, tqdm_kwargs={'disable':False},**{'n_cpus':self.ncores}) # this currently doesn't work...

        return(t)


    def multidigraph(self,T,P):
        t = nx.MultiDiGraph(directed=True)
        for i,reac in self.applied_reactions[T][P].items():
            r = list(reac['e'].reac)
            p = list(reac['e'].prod)
            k = reac['k'] # maybe check equilibrium.as_reactions ( gives forward and backward reactions!)
            if k <= 1: #favours reactants
                t.add_weighted_edges_from([c,i,1/k] for c in r)
                t.add_weighted_edges_from([i,c,k] for c in r)
                t.add_weighted_edges_from([i,c,1/k] for c in p)
                t.add_weighted_edges_from([c,i,k] for c in p)
            elif k >= 1: #favours products
                t.add_weighted_edges_from([c,i,1/k] for c in r)
                t.add_weighted_edges_from([i,c,k] for c in r)
                t.add_weighted_edges_from([i,c,1/k] for c in p)
                t.add_weighted_edges_from([c,i,k] for c in p)
        return(t)

    def generatemultidigraph(self,cost_function=True,mp = False):
        graphs = {}
        for T in self.trange:
            pdict = {}
            for P in self.prange:
                if cost_function:
                    if mp:
                        pdict[P] = self.multidigraph_cost_mp(T,P)
                    else:
                        pdict[P] = self.multidigraph_cost(T, P)
                else:
                    pdict[P] = self.multidigraph(T, P)
            graphs[T] = pdict
        self.graph = graphs
    
    def save(self,filename='graph.p'):
        pickle.dump(self.graph,open(filename,'wb'))
        print('graph saved to: {}'.format(filename))

class GenerateInitialConcentrations:
    
    '''all concentrations are given in x10^-6 (ppm)'''

    def __init__(self,graph=None,compounds=None):
        if graph:
            self.graph = graph
            self.T = list(graph)[0] # dummy temp
            self.P = list(graph[self.T])[0] # dummy pressure
        elif compounds:
            self.compounds = compounds
        else:
            print('need graph or list of compounds')

    def all_random(self,include_co2=True):
        try:
            hasattr(self.compounds,'compounds')
            compounds = self.compounds
        except Exception:
            compounds = [node for node in self.graph[self.T][self.P].nodes() if isinstance(node,str)]
        ic = {c:np.random.random()/1e5 for c in compounds}
        if not include_co2:
            ic['CO2'] = 1
        self.ic = ic
    
    def all_zero(self,include_co2=True):
        try:
            hasattr(self.compounds,'compounds')
            compounds = self.compounds
        except Exception:
            compounds = [node for node in self.graph[self.T][self.P].nodes() if isinstance(node,str)]
        ic = {c:0 for c in compounds}
        if not include_co2:
            ic['CO2'] = 1
        self.ic = ic

    def specific_random(self,compounds=None):
        try:
            hasattr(self.compounds,'compounds')
        except Exception:
            full_list = [n for n in self.graph[self.T][self.P].nodes() if isinstance(n,str)]
        ic = {}
        for c in full_list:
            if c  in self.compounds:
                ic[c] = np.random.random()/1e6
            else:
                ic[c] = 0 
        ic['CO2'] = 1
        self.ic = ic
    
    def update_ic(self,update_dict,include_co2=True):
        ''' update dict = {'CO2':1e-6,'H2O':300e-5} etc.'''
        try:
            hasattr(self.ic,'ic')
        except Exception:
            self.all_zero(include_co2=include_co2)
        for k,v in update_dict.items():
            self.ic[k] = v

    def from_file(self,file_name):
        nodes = [n for n in self.graph[self.T][self.P].nodes() if isinstance(n,str)]
        file_concentrations = loadfn(file_name)
        loaded_compounds = list(file_concentrations.keys())
        ic = {}
        for c in nodes:
            if c not in loaded_compounds:
                ic[c] = 0
            else:
                ic[c] = file_concentrations[c]
        ic['CO2'] = 1
        self.ic = ic
        return(ic)
