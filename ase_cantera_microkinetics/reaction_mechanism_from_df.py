# -------------------------------------------------------------------------------------
# IMPORTS
# -------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import cantera as ct

from ase_cantera_microkinetics import units
from ase_cantera_microkinetics.cantera_utils import get_std_gibbs_dict
from ase_cantera_microkinetics.reaction_mechanism import (
    NameAnalyzer,
    get_species_from_coeffs_NASA_dict,
    change_reference_energies,
    get_species_from_g0_dict_fixed_T,
    get_surf_reactions_from_g0_act_dict_fixed_T,
)

# -------------------------------------------------------------------------------------
# GET GAS SPECIES FROM DF FIXED T
# -------------------------------------------------------------------------------------

def get_gas_species_from_df_fixed_T(
    filename: str,
    energy_ref_funs,
    temperature_fixed: float,
    structure: str = "gas",
    sheet_g0_gas: str = "g0_gas",
    name_analyzer = NameAnalyzer(),
    units_energy = units.eV / units.molecule,
    change_reference = True,
):
    """
    Get gas species from excel file at fixed temperature.
    """
    # Read Gibbs free energies from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_g0_gas, index_col=0)
    df_dict = df.to_dict()
    g0_dict = df_dict[structure]
    # Get gas species.
    gas_species = get_species_from_g0_dict_fixed_T(
        g0_dict=g0_dict,
        temperature_fixed=temperature_fixed,
        name_analyzer=name_analyzer,
        units_energy=units_energy,
    )
    # Change reference energies.
    if change_reference is True:
        gas_species = change_reference_energies(
            species=gas_species,
            energy_ref_funs=energy_ref_funs,
            name_analyzer=name_analyzer,
        )
    # Return gas species.
    return gas_species

# -------------------------------------------------------------------------------------
# GET ADS SPECIES FROM DF FIXED T
# -------------------------------------------------------------------------------------

def get_ads_species_from_df_fixed_T(
    filename: str,
    structure: str,
    temperature_fixed: float,
    sheet_g0_ads: str = "g0_ads",
    name_analyzer: object = NameAnalyzer(),
    units_energy: float = units.eV/units.molecule,
):
    """
    Get ads species from excel file at fixed temperature.
    """
    # Read Gibbs free energies from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_g0_ads, index_col=0)
    df_dict = df.to_dict()
    g0_dict = df_dict[structure]
    # Get reaction intermediates species.
    ads_species = get_species_from_g0_dict_fixed_T(
        g0_dict=g0_dict,
        temperature_fixed=temperature_fixed,
        name_analyzer=name_analyzer,
        units_energy=units_energy,
    )
    # Return ads species.
    return ads_species

# -------------------------------------------------------------------------------------
# GET TS SPECIES FROM DF FIXED T
# -------------------------------------------------------------------------------------

def get_ts_species_from_df_fixed_T(
    filename: str,
    structure: str,
    temperature_fixed: float,
    sheet_g0_ts: str = "g0_ts",
    units_energy: float = units.eV/units.molecule,
):
    """
    Get transition state species from excel file at fixed temperature.
    """
    # Read Gibbs free energies from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_g0_ts, index_col=0)
    df_dict = df.to_dict()
    g0_dict = df_dict[structure]
    # Get transition state species.
    ts_species = get_species_from_g0_dict_fixed_T(
        g0_dict=g0_dict,
        temperature_fixed=temperature_fixed,
        name_analyzer=None,
        units_energy=units_energy,
    )
    # Return transition state species.
    return ts_species

# -------------------------------------------------------------------------------------
# GET STICKING REACTIONS
# -------------------------------------------------------------------------------------

def get_sticking_reactions(
    e_act_dict: dict,
):
    """
    Get sticking reactions from activation energies dictionary.
    """
    sticking_reactions = []
    for name in e_act_dict:
        if isinstance(e_act_dict[name], str) and "sticking" in e_act_dict[name]:
            sticking_reactions.append(name)
            if "e_act=" in e_act_dict[name]:
                e_act_dict[name] = float(e_act_dict[name].split("e_act=")[1])
            else:
                e_act_dict[name] = 0.
    return sticking_reactions

# -------------------------------------------------------------------------------------
# GET SURF REACTIONS FROM DF FIXED T
# -------------------------------------------------------------------------------------

def get_surf_reactions_from_df_fixed_T(
    filename: str,
    gas: ct.Solution,
    site_density: float,
    structure: str,
    temperature_fixed: float,
    sheet_g0_ads: str = "g0_ads",
    sheet_g0_ts: str = "g0_ts",
    name_analyzer: object = NameAnalyzer(),
    units_energy: float = units.eV / units.molecule,
):
    """
    Get surface reactions from excel file at fixed temperature.
    """
    # Read Gibbs free energies of transition states from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_g0_ts, index_col=0)
    df_dict = df.to_dict()
    g0_ts_dict = df_dict[structure]
    # Read Gibbs free energies of reaction intermediates from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_g0_ads, index_col=0)
    df_dict = df.to_dict()
    g0_ads_gas_dict = df_dict[structure]
    g0_ads_gas_dict.update(get_std_gibbs_dict(
        phase=gas,
        units_energy=units_energy,
    ))
    # Get activation Gibbs free energies dictionary.
    g0_act_dict = {}
    for name in g0_ts_dict:
        g0_act_dict[name] = g0_ts_dict[name]
        if isinstance(g0_ts_dict[name], (int, float, np.integer, np.floating)):
            for reactant in name_analyzer.get_reactants(name=name):
                g0_act_dict[name] -= g0_ads_gas_dict[reactant]
    # Get sticking reactions.
    sticking_reactions = get_sticking_reactions(g0_act_dict)
    # Get surface reactions.
    surf_reactions = get_surf_reactions_from_g0_act_dict_fixed_T(
        gas=gas,
        g0_act_dict=g0_act_dict,
        temperature_fixed=temperature_fixed,
        name_analyzer=name_analyzer,
        site_density=site_density,
        sticking_reactions=sticking_reactions,
        units_energy=units_energy,
    )
    return surf_reactions

# -------------------------------------------------------------------------------------
# GET GAS SPECIES FROM DF COEFFS NASA
# -------------------------------------------------------------------------------------

def get_gas_species_from_df_coeffs_NASA(
    filename: str,
    energy_ref_funs: dict,
    names: list = "all",
    sheet_coeffs: str = "coeffs_gas",
    name_analyzer: object = NameAnalyzer(),
    units_energy: float = units.eV / units.molecule,
    change_reference: bool = True,
):
    """
    Get gas species from excel file with NASA polynomial coefficients.
    """
    # Read NASA coefficients from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_coeffs, index_col=0)
    df_dict = df.to_dict(orient="index")
    if names == "all":
        names = list(df_dict.keys())
    coeffs_NASA_dict = {name: list(df_dict[name].values()) for name in names}
    # Get gas species.
    gas_species = get_species_from_coeffs_NASA_dict(
        coeffs_NASA_dict=coeffs_NASA_dict,
        e_form_dict=None,
        units_energy=units_energy,
        name_analyzer=name_analyzer,
    )
    # Change reference energies.
    if change_reference is True:
        gas_species = change_reference_energies(
            species=gas_species,
            energy_ref_funs=energy_ref_funs,
            name_analyzer=name_analyzer,
        )
    # Return gas species.
    return gas_species

# -------------------------------------------------------------------------------------
# GET ADS SPECIES FROM DF COEFFS NASA
# -------------------------------------------------------------------------------------

def get_ads_species_from_df_coeffs_NASA(
    filename: str,
    structure: str,
    sheet_coeffs: str = "coeffs_ads",
    sheet_e_ads: str = "e_ads",
    names: list = "all",
    name_analyzer: object = NameAnalyzer(),
    units_energy: float = units.eV / units.molecule,
):
    """
    Get reaction intermediates species from excel file with NASA coefficients.
    """
    # Read NASA coefficients from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_coeffs, index_col=0)
    df_dict = df.to_dict(orient="index")
    coeffs_NASA_dict = {
        name: list(df_dict[name].values()) 
        for name in df_dict if names == "all" or name in names
    }
    # Get adsorption energies from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_e_ads, index_col=0)
    df_dict = df.to_dict()
    e_ads_dict = df_dict[structure]
    # Get reaction intermediates species.
    ads_species = get_species_from_coeffs_NASA_dict(
        coeffs_NASA_dict=coeffs_NASA_dict,
        e_form_dict=e_ads_dict,
        units_energy=units_energy,
        name_analyzer=name_analyzer,
    )
    return ads_species

# -------------------------------------------------------------------------------------
# GET TS SPECIES FROM DF COEFFS NASA
# -------------------------------------------------------------------------------------

def get_ts_species_from_df_coeffs_NASA(
    filename: str,
    structure: str,
    sheet_coeffs: str = "coeffs_ts",
    sheet_e_ts: str = "e_ts",
    names: list = "all",
    units_energy: float = units.eV / units.molecule,
):
    """
    Get transition state species from excel file with NASA coefficients.
    """
    # Read NASA coefficients from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_coeffs, index_col=0)
    df_dict = df.to_dict(orient="index")
    coeffs_NASA_dict = {
        name: list(df_dict[name].values()) 
        for name in df_dict if names == "all" or name in names
    }
    # Get transition state energies from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_e_ts, index_col=0)
    df_dict = df.to_dict()
    e_ts_dict = df_dict[structure]
    # Get transition state species.
    ts_species = get_species_from_coeffs_NASA_dict(
        coeffs_NASA_dict=coeffs_NASA_dict,
        e_form_dict=e_ts_dict,
        units_energy=units_energy,
        name_analyzer=None,
    )
    return ts_species

# -------------------------------------------------------------------------------------
# GET SURF REACTIONS FROM DF COEFFS NASA
# -------------------------------------------------------------------------------------

def get_surf_reactions_from_df_coeffs_NASA(
    filename: str,
    gas: ct.Solution,
    site_density: float,
    structure: str,
    sheet_e_ads: str = "e_ads",
    sheet_e_ts: str = "e_ts",
    name_analyzer: object = NameAnalyzer(),
    units_energy: float = units.eV/units.molecule,
):
    """
    Get surface reactions from excel file with NASA coefficients.
    """
    # Get trantisiton state energies from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_e_ts, index_col=0)
    df_dict = df.to_dict()
    e_ts_dict = df_dict[structure]
    # Get adsorption energies from excel file.
    df = pd.read_excel(filename, sheet_name=sheet_e_ads, index_col=0)
    df_dict = df.to_dict()
    e_form_dict = df_dict[structure]
    for spec in gas.species():
        e_form_dict[spec.name] = spec.thermo.coeffs[6] * ct.gas_constant / units_energy
    # Get activation energies dictionary.
    e_act_dict = {}
    for name in e_ts_dict:
        e_act_dict[name] = e_ts_dict[name]
        if isinstance(e_ts_dict[name], (int, float, np.integer, np.floating)):
            for reactant in name_analyzer.get_reactants(name=name):
                e_act_dict[name] -= e_form_dict[reactant]
    # Get sticking reactions.
    sticking_reactions = get_sticking_reactions(e_act_dict)
    # Get surface reactions.
    surf_reactions = get_surf_reactions_from_g0_act_dict_fixed_T(
        gas=gas,
        g0_act_dict=e_act_dict,
        temperature_fixed=0.,
        name_analyzer=name_analyzer,
        site_density=site_density,
        sticking_reactions=sticking_reactions,
        units_energy=units_energy,
    )
    return surf_reactions

# -------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------