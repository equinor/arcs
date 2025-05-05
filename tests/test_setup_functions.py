import numpy as np 
import json 
from monty.serialization import loadfn
from chempy import Equilibrium 
from arcs.setup_functions import ReactionGibbsandEquilibrium


def test_json_data_loader():
    json.load(open('hse06_dft_data.json'))

def test_Gibbs():
    dft_data = loadfn('hse06_dft_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=dft_data)
    gibbs_free_energy = rge.Gibbs(compound='CO2',temperature=100,pressure=1)
    assert gibbs_free_energy == -31.974842281859573

def test_reaction_gibbs():
    dft_data = loadfn('hse06_dft_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=dft_data)
    reaction = 'CO + H2O = H2 + CO2'
    eq = Equilibrium.from_string(reaction)
    reaction_gibbs = rge.reaction_gibbs(
        reaction=eq,
        pressure=1,
        temperature=100
    )
    assert reaction_gibbs == -0.6184489112909475



    

