import os
import gzip
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from ase.thermochemistry import IdealGasThermo
from scipy.constants import Boltzmann, e
import numpy as np
from chempy.reactionsystem import Substance
from tqdm import tqdm
import networkx as nx
from pathos.helpers import mp as pmp
import math
import pickle


def get_compound_directory(base, compound, size):
    return os.path.join(base, compound, size)


class GetEnergyandVibrationsVASP:
    """Class to get the Total Energy and Vibrations from a directory containing a calculations"""

    def __init__(self, relax_directory, vibrations_directory):
        self.relax = relax_directory
        self.vibrations = vibrations_directory

    def atoms(self):
        def _get_initial_magnetic_moments(ats):
            magmoms = []
            for atnum in ats.get_atomic_numbers():
                magmoms.append([0 if atnum % 2 == 0 else 1][0])
            return magmoms

        structure = read("{}/POSCAR.gz".format(self.relax))
        structure.set_initial_magnetic_moments(_get_initial_magnetic_moments(structure))
        return structure

    def energy(self):
        outcar = gzip.open("{}/OUTCAR.gz".format(self.relax), "tr").readlines()
        for line in outcar:
            if "y=" in line:
                energy = float(line.split()[-1])
        if len(list(dict.fromkeys(self.atoms().get_atomic_numbers()))) == 1:
            if not any(
                x in self.atoms().symbols.get_chemical_formula()
                for x in ["H2", "O2", "N2"]
            ):
                energy = energy / self.atoms().get_global_number_of_atoms()

        return energy

    def spin(self):
        if self.atoms().get_chemical_formula() in ["O2", "CO3"]:
            return 1
        else:
            outcar = gzip.open("{}/OUTCAR.gz".format(self.relax), "tr").readlines()
            for line in outcar:
                if "NELECT" in line:
                    nelect = float(line.split()[2])
                    return [0 if nelect % 2 == 0 else 0.5][0]

    def pointgroup(self):
        atoms = self.atoms()
        pg = PointGroupAnalyzer(AseAtomsAdaptor.get_molecule(atoms)).get_pointgroup()
        return pg.sch_symbol

    def islinear(self):
        num_at = self.atoms()
        if num_at.get_global_number_of_atoms() == 1:
            return "monatomic"
        else:
            pg = self.pointgroup()
            if "*" in pg:
                return "linear"
            else:
                return "nonlinear"

    def rotation_num(self):
        pg = [x for x in self.pointgroup()]
        if pg[0] == "C":
            if pg[-1] == "i":
                rot = 1
            elif pg[-1] == "s":
                rot = 1
            elif pg[-1] == "h":
                rot = int(pg[1])
            elif pg[-1] == "v":
                if pg[1] == "*":
                    rot = 1
                else:
                    rot = int(pg[1])
            elif len(pg) == 2:
                rot = int(pg[-1])

        elif pg[0] == "D":
            if pg[-1] == "h":
                if pg[1] == "*":
                    rot = 2
                else:
                    rot = 2 * int(pg[1])
            elif pg[-1] == "d":
                rot = 2 * int(pg[1])
            elif len(pg) == 2:
                rot = 2 * int(pg[1])

        elif pg[0] == "T":
            rot = 12

        elif pg[0] == "O":
            rot = 24

        elif pg[0] == "I":
            rot = 60

        return rot

    def get_vibrations(self):
        outcar = gzip.open("{}/OUTCAR.gz".format(self.vibrations), "tr").readlines()
        frequencies = []
        for line in outcar:
            if "THz" in line:
                if "f/i" not in line:  # we ignore imaginary modes
                    ev = float(line.split()[-2]) / 1000
                    frequencies.append(ev)
        return frequencies

    def as_dict(self):
        return {
            "atoms": self.atoms(),
            "pointgroup": self.pointgroup(),
            "spin": self.spin(),
            "rotation_num": self.rotation_num(),
            "islinear": self.islinear(),
            "energy": self.energy(),
            "vibrations": self.get_vibrations(),
        }


class GetEnergyandVibrationsDalton:
    def __init__(self, relax_directory, vibrations_directory):
        self.relax = relax_directory
        self.vibrations = vibrations_directory

    def atoms(self):
        structure = read("{}/output.xyz".format(self.relax))
        return structure

    def energy(self):
        daltonout = open("{}/output.out".format(self.relax), "r").readlines()
        for line in daltonout:
            if "Final" in line and "energy" in line:
                energy = float(line.split()[-1])
        return energy * 27.211396  # a.u -> eV conversion

    def spin(self):
        nelect = np.sum(
            self.atoms().get_atomic_numbers()
        )  # not sure if this is correct
        return [0 if nelect % 2 == 0 else 1][0]

    def pointgroup(self):
        atoms = self.atoms()
        pg = PointGroupAnalyzer(AseAtomsAdaptor.get_molecule(atoms)).get_pointgroup()
        return pg.sch_symbol

    def islinear(self):
        num_at = self.atoms()
        if num_at.get_global_number_of_atoms() == 1:
            return "monatomic"
        else:
            pg = self.pointgroup()
            if "*" in pg:
                return "linear"
            else:
                return "nonlinear"

    def rotation_num(self):
        pg = [x for x in self.pointgroup()]
        if pg[0] == "C":
            if pg[-1] == "i":
                rot = 1
            elif pg[-1] == "s":
                rot = 1
            elif pg[-1] == "h":
                rot = int(pg[1])
            elif pg[-1] == "v":
                if pg[1] == "*":
                    rot = 1
                else:
                    rot = int(pg[1])
            elif len(pg) == 2:
                rot = int(pg[-1])

        elif pg[0] == "D":
            if pg[-1] == "h":
                if pg[1] == "*":
                    rot = 2
                else:
                    rot = 2 * int(pg[1])
            elif pg[-1] == "d":
                rot = 2 * int(pg[1])
            elif len(pg) == 2:
                rot = 2 * int(pg[1])

        elif pg[0] == "T":
            rot = 12

        elif pg[0] == "O":
            rot = 24

        elif pg[0] == "I":
            rot = 60

        return rot

    def get_vibrations(self):
        daltonout = open("{}/output.out".format(self.vibrations), "r").readlines()
        datalines = []
        for i, line in enumerate(daltonout):
            if "frequency" in line and "mode" in line:
                datalines.append(i + 4)
            elif "Normal Coordinates" in line:
                datalines.append(i)
        recip_cm = [
            float(x.split()[2]) / 8100
            for x in daltonout[datalines[0] : datalines[1]]
            if len(x.split()) > 1
        ]  # 8100 conversion cm-1 -> eV
        return recip_cm

    def as_dict(self):
        return {
            "atoms": self.atoms(),
            "spin": self.spin(),
            "rotation_num": self.rotation_num(),
            "islinear": self.islinear(),
            "energy": self.energy(),
            "vibrations": self.get_vibrations(),
        }


class ReactionGibbsandEquilibrium:
    def __init__(self, reaction, temperature, pressure, reaction_input):
        self.reaction = reaction
        self.temperature = temperature
        self.pressure = pressure * 100000  # pressure in bar -> Pa
        self.reaction_input = reaction_input

    def Gibbs(self, c):
        data = self.reaction_input[c]
        igt = IdealGasThermo(
            vib_energies=data.get_vibrations(),
            geometry=data.islinear(),
            potentialenergy=data.energy(),
            atoms=data.atoms(),
            symmetrynumber=data.rotation_num(),
            spin=data.spin(),
            natoms=data.atoms().get_global_number_of_atoms(),
        )
        G = igt.get_gibbs_energy(self.temperature, self.pressure, verbose=False)
        H = igt.get_enthalpy(self.temperature, verbose=False)
        S = igt.get_entropy(self.temperature, self.pressure, verbose=False)
        Z = igt.get_entropy(self.temperature, self.pressure, verbose=False)
        return {"G": G, "H": H, "S": S, "Z": Z}

    def reaction_energy(self):
        prod = self.reaction.prod
        reac = self.reaction.reac
        reaction_compounds = list(prod) + list(reac)
        # need to add a charge neutrality condition and mass balance violation
        gibbs = {c: self.Gibbs(c) for c in reaction_compounds}
        prod_sum = np.sum([self.Gibbs(c)["G"] * prod[c] for c in gibbs if c in prod])
        reac_sum = np.sum([self.Gibbs(c)["G"] * reac[c] for c in gibbs if c in reac])
        return float(prod_sum - reac_sum)

    def equilibrium_constant(self):
        K = np.exp(-(self.reaction_energy() * e) / (Boltzmann * self.temperature))
        return K

    def as_dict(self):
        return {
            "G_react": self.reaction_energy(),
            "K_react": self.equilibrium_constant(),
        }


class ApplyDataToReaction:
    # future functionality:
    # 1. interpolation - should speed things up
    """this class applies the gibbs data to a specific reaction"""

    def __init__(self, trange, prange, reactions, compound_data, nprocs):
        self.trange = trange
        self.prange = prange
        self.reactions = {
            i: r for i, r in enumerate(pickle.load(open(reactions, "rb")))
        }
        self.compound_data = compound_data
        self.nprocs = nprocs
        self.barformat = "{desc:<20}{percentage:3.0f}%|{bar:10}{r_bar}"

    def _generate_data_serial(self, t, p):  # serial
        reactions = {
            i: {
                "e": r,
                "k": ReactionGibbsandEquilibrium(
                    r, t, p, self.compound_data
                ).equilibrium_constant(),
                "g": ReactionGibbsandEquilibrium(
                    r, t, p, self.compound_data
                ).reaction_energy(),
            }
            for i, r in tqdm(self.reactions.items())
        }
        return reactions

    def generate_data(self, t, p):  # multiprocessed
        manager = pmp.Manager()
        queue = manager.Queue()

        def mp_function(reaction_keys, out_q):
            data = {}
            for r in reaction_keys:
                rge = ReactionGibbsandEquilibrium(
                    self.reactions[r], t, p, self.compound_data
                )
                data[r] = {
                    "e": self.reactions[r],
                    "k": rge.equilibrium_constant(),
                    "g": rge.reaction_energy(),
                }
            out_q.put(data)

        resultdict = {}
        r_keys = list(self.reactions.keys())
        chunksize = int(math.ceil(len(self.reactions) / float(self.nprocs)))
        processes = []

        for i in range(self.nprocs):
            pr = pmp.Process(
                target=mp_function,
                args=(r_keys[chunksize * i : chunksize * (i + 1)], queue),
            )
            processes.append(pr)
            pr.start()

        for i in range(self.nprocs):
            resultdict.update(queue.get(timeout=1800))

        for pr in processes:
            pr.join()

        return resultdict

    def apply(self, serial=False):
        data = {}
        with tqdm(
            total=len(self.trange) * len(self.prange), bar_format=self.barformat
        ) as pbar:
            for t in self.trange:
                pdat = {}
                for p in self.prange:
                    if serial is True:
                        pdat[p] = self._generate_data_serial(t, p)
                    else:
                        pdat[p] = self.generate_data(t, p)
                    pbar.update(1)
                data[t] = pdat
        self.data = data
        return self.data

    def save(self, filename="applied_reactions.p"):
        pickle.dump(self.data, open(filename, "wb"))
        print("data saved to: {}".format(filename))


class GraphGenerator:
    def __init__(self, applied_reactions):
        self.applied_reactions = pickle.load(open(applied_reactions, "rb"))
        self.trange = list(self.applied_reactions)
        self.prange = list(self.applied_reactions[self.trange[0]])
        self.barformat = "{desc:<20}{percentage:3.0f}%|{bar:10}{r_bar}"

    def _cost_function(self, gibbs, T, reactants):
        """takes the cost function that is used in https://www.nature.com/articles/s41467-021-23339-x.pdf"""

        comps = []
        for r, n in reactants.items():
            for i in range(n):
                comps.append(r)

        num_atoms = np.sum(
            [
                np.sum([y for x, y in Substance.from_formula(c).composition.items()])
                for c in comps
            ]
        )

        return np.log(1 + (273 / T) * np.exp(gibbs / num_atoms / 1))

    def multidigraph_cost(self, T, P):
        """this will weight the graph in terms of a cost function which makes it better for a Djikstra algorithm to work"""
        t = nx.MultiDiGraph(directed=True)
        for i, reac in self.applied_reactions[T][P].items():
            f_cost = self._cost_function(reac["g"], T, reac["e"].reac)  # forward cost
            b_cost = self._cost_function(-reac["g"], T, reac["e"].prod)  # backward cost
            r = list(reac["e"].reac)
            p = list(reac["e"].prod)
            t.add_weighted_edges_from(
                [c, i, f_cost] for c in r
            )  # reactants -> reaction
            t.add_weighted_edges_from(
                [i, c, b_cost] for c in r
            )  # reaction -> reactants
            t.add_weighted_edges_from([i, c, f_cost] for c in p)  # reaction -> products
            t.add_weighted_edges_from([c, i, b_cost] for c in p)  # products -> reaction
        return t

    def multidigraph(self, T, P):
        t = nx.MultiDiGraph(directed=True)
        for i, reac in self.applied_reactions[T][P].items():
            r = list(reac["e"].reac)
            p = list(reac["e"].prod)
            k = reac[
                "k"
            ]  # maybe check equilibrium.as_reactions ( gives forward and backward reactions!)
            if k <= 1:  # favours reactants
                t.add_weighted_edges_from([c, i, 1 / k] for c in r)
                t.add_weighted_edges_from([i, c, k] for c in r)
                t.add_weighted_edges_from([i, c, 1 / k] for c in p)
                t.add_weighted_edges_from([c, i, k] for c in p)
            elif k >= 1:  # favours products
                t.add_weighted_edges_from([c, i, 1 / k] for c in r)
                t.add_weighted_edges_from([i, c, k] for c in r)
                t.add_weighted_edges_from([i, c, 1 / k] for c in p)
                t.add_weighted_edges_from([c, i, k] for c in p)
        return t

    def generatemultidigraph(self, cost_function=True):
        graphs = {}
        with tqdm(
            total=len(self.trange) * len(self.prange), bar_format=self.barformat
        ) as pbar:
            pbar.set_description("generating graph")
            for T in self.trange:
                pdict = {}
                for P in self.prange:
                    if cost_function is True:
                        pdict[P] = self.multidigraph_cost(T, P)
                    elif cost_function is False:
                        pdict[P] = self.multidigraph(T, P)
                    pbar.update(1)
                graphs[T] = pdict
        self.graph = graphs

    def save(self, filename="graph.p"):
        pickle.dump(self.graph, open(filename, "wb"))
        print("graph saved to: {}".format(filename))
