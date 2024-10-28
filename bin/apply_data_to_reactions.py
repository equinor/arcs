from arcs.setup_functions import GetEnergyandVibrationsVASP, ApplyDataToReaction
import warnings
import os

nprocs = 32
functional = "SCAN"
trange = range(10, 100, 4)
prange = range(1, 100, 4)

if __name__ == "__main__":
    dft_directory = "../../../{}/20/new/".format(functional)
    compounds = [
        "CO2",
        "H2O",
        "H2S",
        "CO",
        "NO2",
        "H2CO",
        "CH3CHO",
        "O2",
        "H2",
        "CH4",
        "NH3",
        "H2SO4",
        "S8",
        "H2SO3",
        "HNO3",
        "HNO2",
        "NO",
        "NOHSO4",
        "SO2",
        "CH3OH",
        "CH3CH2OH",
        "N2",
        "CH2O2",
        "CH3COOH",
    ]
    dft_scraped_data = {
        c: GetEnergyandVibrationsVASP(
            os.path.join(dft_directory, c, "relax"),
            os.path.join(dft_directory, c, "vibrations"),
        )
        for c in compounds
    }
    warnings.simplefilter("ignore")

    adr = ApplyDataToReaction(
        trange=trange,
        prange=prange,
        reactions="reactions.p",
        compound_data=dft_scraped_data,
        nprocs=nprocs,
    )
    adr.apply()
    adr.save(filename=open("{}_reactions.p".format(functional), "wb"))
