import numpy as np 
import json 
from monty.serialization import loadfn
from chempy import Equilibrium 
from arcs.setup_functions import ReactionGibbsandEquilibrium


def test_json_data_loader():
    json.load(open('test_dft_data.json'))

def test_Gibbs():
    dft_data = loadfn('test_dft_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=dft_data)
    gibbs_free_energy = rge.Gibbs(dft_dict=dft_data['CO2'],temperature=100,pressure=1)
    assert gibbs_free_energy == -31.974842281859573

def test_reaction_gibbs():
    dft_data = loadfn('test_dft_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=dft_data)
    reaction = 'CO + H2O = H2 + CO2'
    eq = Equilibrium.from_string(reaction)
    reaction_gibbs = rge.reaction_gibbs(
        reaction=eq,
        pressure=1,
        temperature=100
    )
    assert reaction_gibbs == -0.6184489112909475

def test_equilibrium_constant():
    dft_data = loadfn('test_dft_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=dft_data)
    equilibrium_constant = rge.equilibrium_constant(gibbs_free_energy=-0.6184489112909475,
                                                    temperature=100)
    assert equilibrium_constant == 1.4738501117040624e+31

def test_get_reaction_gibbs_and_equilibrium():
    dft_data = loadfn('test_dft_data.json')
    reaction = 'CO + H2O = H2 + CO2'
    rge = ReactionGibbsandEquilibrium(reaction_input=dft_data)
    _dict = rge.get_reaction_gibbs_and_equilibrium(reaction=Equilibrium.from_string(reaction),temperature=100,pressure=1)

    assert _dict == {'g': -0.6184489112909475, 'k': 1.4738501117040624e+31}





    

