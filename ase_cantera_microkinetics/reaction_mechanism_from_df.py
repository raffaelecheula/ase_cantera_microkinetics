# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

import pandas as pd
import cantera as ct
from ase_cantera_microkinetics import units
from ase_cantera_microkinetics.cantera_utils import get_std_gibbs_dict
from ase_cantera_microkinetics.reaction_mechanism import (
    NameAnalyzer,
    get_species_from_coeffs_NASA_dict,
    coeffs_NASA_dict_from_dataframe,
    change_reference_energies,
    get_species_from_g0_dict_fixed_T,
    get_surf_reactions_from_g0_act_dict_fixed_T,
)

# -----------------------------------------------------------------------------
# GET GAS SPECIES FROM DF FIXED T
# -----------------------------------------------------------------------------

def get_gas_species_from_df_fixed_T(
    filename,
    energy_ref_funs,
    temperature_fixed,
    structure = 'gas',
    sheet_g0_gas = 'g0_gas',
    name_analyzer = NameAnalyzer(),
    units_energy = units.eV/units.molecule,
):

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_g0_gas)

    df_dict = df.to_dict()
    g0_dict = df_dict[structure]

    gas_species = get_species_from_g0_dict_fixed_T(
        g0_dict = g0_dict,
        temperature_fixed = temperature_fixed,
        name_analyzer = name_analyzer,
        composition_dict = 'auto',
        size_dict = 'auto',
        units_energy = units_energy,
    )

    gas_species = change_reference_energies(
        species = gas_species,
        energy_ref_funs = energy_ref_funs,
        name_analyzer = name_analyzer,
    )

    return gas_species

# -----------------------------------------------------------------------------
# GET ADS SPECIES FROM DF FIXED T
# -----------------------------------------------------------------------------

def get_ads_species_from_df_fixed_T(
    filename,
    structure,
    temperature_fixed,
    sheet_g0_ads = 'g0_ads',
    name_analyzer = NameAnalyzer(),
    units_energy = units.eV/units.molecule,
):

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_g0_ads)

    df_dict = df.to_dict()
    g0_dict = df_dict[structure]

    ads_species = get_species_from_g0_dict_fixed_T(
        g0_dict = g0_dict,
        temperature_fixed = temperature_fixed,
        name_analyzer = name_analyzer,
        composition_dict = 'auto',
        size_dict = 'auto',
        units_energy = units_energy,
    )

    return ads_species

# -----------------------------------------------------------------------------
# GET TS SPECIES FROM DF FIXED T
# -----------------------------------------------------------------------------

def get_ts_species_from_df_fixed_T(
    filename,
    structure,
    temperature_fixed,
    sheet_g0_ts = 'g0_ts',
    units_energy = units.eV/units.molecule,
):

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_g0_ts)

    df_dict = df.to_dict()
    g0_dict = df_dict[structure]

    ts_species = get_species_from_g0_dict_fixed_T(
        g0_dict = g0_dict,
        temperature_fixed = temperature_fixed,
        composition_dict = None,
        size_dict = None,
        units_energy = units_energy,
    )

    return ts_species

# -----------------------------------------------------------------------------
# GET SURF REACTIONS FROM DF FIXED T
# -----------------------------------------------------------------------------

def get_surf_reactions_from_df_fixed_T(
    filename,
    gas,
    site_density,
    structure,
    temperature_fixed,
    sheet_g0_ads = 'g0_ads',
    sheet_g0_ts = 'g0_ts',
    name_analyzer = NameAnalyzer(),
    units_energy = units.eV/units.molecule,
):

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_g0_ts)

    df_dict = df.to_dict()
    g0_ts_dict = df_dict[structure]

    sticking_reactions = [
        name for name in df_dict['sticking'] if df_dict['sticking'][name]
    ]

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_g0_ads)

    df_dict = df.to_dict()
    g0_ads_gas_dict = df_dict[structure]

    g0_ads_gas_dict.update(get_std_gibbs_dict(
        phase = gas,
        units_energy = units_energy,
    ))

    g0_act_dict = {}
    for name in g0_ts_dict:
        g0_act_dict[name] = g0_ts_dict[name]
        for reactant in name_analyzer.get_reactants(name = name):
            g0_act_dict[name] -= g0_ads_gas_dict[reactant]

    surf_reactions = get_surf_reactions_from_g0_act_dict_fixed_T(
        gas = gas,
        g0_act_dict = g0_act_dict,
        temperature_fixed = temperature_fixed,
        name_analyzer = name_analyzer,
        pre_exp_dict = 'auto',
        site_density = site_density,
        sticking_reactions = sticking_reactions,
        units_energy = units_energy,
    )

    return surf_reactions

# -----------------------------------------------------------------------------
# GET GAS SPECIES FROM DF COEFFS NASA
# -----------------------------------------------------------------------------

def get_gas_species_from_df_coeffs_NASA(
    filename,
    energy_ref_funs,
    names = 'all',
    sheet_coeffs = 'coeffs_gas',
    species_key = 'species',
    coeff_key = 'NASA',
    name_analyzer = NameAnalyzer(),
    units_energy = units.eV/units.molecule,
):

    df_coeffs = pd.read_excel(filename, sheet_name = sheet_coeffs)

    coeffs_NASA_dict = coeffs_NASA_dict_from_dataframe(
        dataframe = df_coeffs,
        names = names,
        species_key = species_key,
        coeff_key = coeff_key,
    )

    gas_species = get_species_from_coeffs_NASA_dict(
        coeffs_NASA_dict = coeffs_NASA_dict,
        e_form_dict = None,
        units_energy = units_energy,
        name_analyzer = name_analyzer,
        composition_dict = 'auto',
        size_dict = 'auto',
    )

    gas_species = change_reference_energies(
        species = gas_species,
        energy_ref_funs = energy_ref_funs,
        name_analyzer = name_analyzer,
    )

    return gas_species

# -----------------------------------------------------------------------------
# GET ADS SPECIES FROM DF COEFFS NASA
# -----------------------------------------------------------------------------

def get_ads_species_from_df_coeffs_NASA(
    filename,
    structure,
    sheet_coeffs = 'coeffs_ads',
    sheet_e_ads = 'e_ads',
    names = 'all',
    species_key = 'species',
    coeff_key = 'NASA',
    name_analyzer = NameAnalyzer(),
    units_energy = units.eV/units.molecule,
):

    df_coeffs = pd.read_excel(filename, sheet_name = sheet_coeffs)

    coeffs_NASA_dict = coeffs_NASA_dict_from_dataframe(
        dataframe = df_coeffs,
        names = names,
        species_key = species_key,
        coeff_key = coeff_key,
    )

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_e_ads)

    df_dict = df.to_dict()
    e_ads_dict = df_dict[structure]
    
    ads_species = get_species_from_coeffs_NASA_dict(
        coeffs_NASA_dict = coeffs_NASA_dict,
        e_form_dict = e_ads_dict,
        units_energy = units_energy,
        name_analyzer = name_analyzer,
        composition_dict = 'auto',
        size_dict = 'auto',
    )

    return ads_species

# -----------------------------------------------------------------------------
# GET TS SPECIES FROM DF COEFFS NASA
# -----------------------------------------------------------------------------

def get_ts_species_from_df_coeffs_NASA(
    filename,
    structure,
    sheet_coeffs = 'coeffs_ts',
    sheet_e_ts = 'e_ts',
    names = 'all',
    species_key = 'species',
    coeff_key = 'NASA',
    units_energy = units.eV/units.molecule,
):

    df_coeffs = pd.read_excel(filename, sheet_name = sheet_coeffs)

    coeffs_NASA_dict = coeffs_NASA_dict_from_dataframe(
        dataframe = df_coeffs,
        names = names,
        species_key = species_key,
        coeff_key = coeff_key,
    )

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_e_ts)

    df_dict = df.to_dict()
    e_ts_dict = df_dict[structure]

    ts_species = get_species_from_coeffs_NASA_dict(
        coeffs_NASA_dict = coeffs_NASA_dict,
        e_form_dict = e_ts_dict,
        units_energy = units_energy,
        composition_dict = None,
        size_dict = None,
    )

    return ts_species

# -----------------------------------------------------------------------------
# GET SURF REACTIONS FROM DF COEFFS NASA
# -----------------------------------------------------------------------------

def get_surf_reactions_from_df_coeffs_NASA(
    filename,
    gas,
    site_density,
    structure,
    sheet_e_ads = 'e_ads',
    sheet_e_ts = 'e_ts',
    name_analyzer = NameAnalyzer(),
    units_energy = units.eV/units.molecule,
):

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_e_ts)

    df_dict = df.to_dict()
    e_ts_dict = df_dict[structure]

    sticking_reactions = [
        name for name in df_dict['sticking'] if df_dict['sticking'][name]
    ]

    df = pd.read_excel(filename, index_col = 0, sheet_name = sheet_e_ads)

    df_dict = df.to_dict()
    e_form_dict = df_dict[structure]

    for spec in gas.species():
        e_form_dict[spec.name] = (
            spec.thermo.coeffs[6]*ct.gas_constant/units_energy
        )

    e_act_dict = {}
    for name in e_ts_dict:
        e_act_dict[name] = e_ts_dict[name]
        for reactant in name_analyzer.get_reactants(name = name):
            e_act_dict[name] -= e_form_dict[reactant]

    surf_reactions = get_surf_reactions_from_g0_act_dict_fixed_T(
        gas = gas,
        g0_act_dict = e_act_dict,
        temperature_fixed = 0.,
        name_analyzer = name_analyzer,
        pre_exp_dict = 'auto',
        site_density = site_density,
        sticking_reactions = sticking_reactions,
        units_energy = units_energy,
    )

    return surf_reactions

# -----------------------------------------------------------------------------
# END
# -----------------------------------------------------------------------------