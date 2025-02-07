# -------------------------------------------------------------------------------------
# IMPORTS
# -------------------------------------------------------------------------------------

import yaml
import numpy as np
import cantera as ct

from ase_cantera_microkinetics import units
from ase_cantera_microkinetics.cantera_utils import reactions_from_cat_ts
from ase_cantera_microkinetics.reaction_mechanism import (
    get_spec_NasaPoly2,
    get_surf_reactions_from_g0_act_dict_fixed_T,
)

# -------------------------------------------------------------------------------------
# GET SPECIES FROM YAML NASA7
# -------------------------------------------------------------------------------------

def get_species_from_yaml_NASA7(
    species_list,
    e_form_dict={},
    units_energy=units.eV/units.molecule,
):
    # Get species from YAML dictionary.
    species_cantera = []
    for species_dict in species_list:
        T_list = species_dict["thermo"]["temperature-ranges"]
        species = get_spec_NasaPoly2(
            name=species_dict["name"],
            composition=species_dict["composition"],
            size=species_dict.get("size", None),
            e_form=e_form_dict.get(species_dict["name"], None),
            coeffs_NASA=np.concatenate(species_dict["thermo"]["data"]),
            units_energy=units_energy,
            T_low=T_list[0],
            T_mid=T_list[1] if len(T_list) == 3 else np.mean(T_list),
            T_high=T_list[-1],
        )
        species_cantera.append(species)
    return species_cantera

# -------------------------------------------------------------------------------------
# GET MECHANISM FROM YAML
# -------------------------------------------------------------------------------------

def get_mechanism_from_yaml(
    yaml_file='mechanism.yaml',
    e_form_dict=None,
    site_density=1e-8,
    temperature=298.15,
    units_energy=units.eV/units.molecule,
):
    # Read the reaction mechanism.
    with open(yaml_file, 'r') as fileobj:
        mechanism = yaml.safe_load(fileobj)
    # Get gas species.
    gas_species = get_species_from_yaml_NASA7(
        species_list=mechanism["species-gas"],
        e_form_dict=e_form_dict,
        units_energy=units_energy,
    )
    # Get adsorbates species.
    ads_species = get_species_from_yaml_NASA7(
        species_list=mechanism["species-adsorbates"],
        e_form_dict=e_form_dict,
        units_energy=units_energy,
    )
    # Get transition state species.
    ts_species = get_species_from_yaml_NASA7(
        species_list=mechanism["species-reactions"],
        e_form_dict=e_form_dict,
        units_energy=units_energy,
    )
    # Get sticking reactions.
    sticking_reactions = [
        species["name"] for species in mechanism["species-reactions"]
        if species["thermo"].get("sticking", False)
    ]
    # Initialize gas phase.
    gas = ct.Solution(
        name='gas',
        species=gas_species,
        thermo='ideal-gas',
        transport='Mix',
    )
    # Initialize surface reactions.
    surf_reactions = get_surf_reactions_from_g0_act_dict_fixed_T(
        gas=gas,
        g0_act_dict={spec.name: 0. for spec in ts_species},
        temperature_fixed=temperature,
        site_density=site_density,
        sticking_reactions=sticking_reactions,
    )
    # Initialize catalyst phase (adsorbates).
    cat = ct.Interface(
        name='cat',
        species=ads_species,
        reactions=surf_reactions,
        thermo='ideal-surface',
        kinetics='surface',
        adjacent=[gas],
    )
    # Initialize catalyst phase (transition states).
    cat_ts = ct.Interface(
        name='cat_ts',
        species=ts_species,
        thermo='ideal-surface',
        kinetics='surface',
        adjacent=[gas],
    )
    cat.site_density = site_density
    cat_ts.site_density = site_density
    # Get reactions from transition states.
    reactions_from_cat_ts(gas=gas, cat=cat, cat_ts=cat_ts)
    return gas, cat, cat_ts

# -------------------------------------------------------------------------------------
# GET E FORM FROM ASE ATOMS
# -------------------------------------------------------------------------------------

def get_e_form_from_ase_atoms(
    yaml_file,
    atoms_ads_list,
    atoms_ts_list,
    free_site,
    e_form_key="E_form",
    e_index=None,
    species_key="species_cantera",
    e_form_missing=1.,
):
    # Read the reaction mechanism.
    with open(yaml_file, 'r') as fileobj:
        mechanism = yaml.safe_load(fileobj)
    # Get species lists.
    adsorbates_list = [
        species["name"] for species in mechanism["species-adsorbates"]
    ]
    reactions_list = [
        species["name"] for species in mechanism["species-reactions"]
        if not species["thermo"].get("sticking", False)
    ]
    # Get formation energies.
    e_form_dict = {free_site: 0.}
    for species_name in [name for name in adsorbates_list if name != free_site]:
        atoms_list = [
            atoms for atoms in atoms_ads_list
            if atoms.info[species_key] == species_name
        ]
        e_form_list = [
            atoms.info[e_form_key] 
            if e_index is None else atoms.info[e_form_key][e_index]
            for atoms in atoms_list
        ]
        e_form = np.min(e_form_list) if len(e_form_list) > 0 else e_form_missing
        e_form_dict[species_name] = e_form
    for species_name in reactions_list:
        atoms_list = [
            atoms for atoms in atoms_ts_list
            if atoms.info[species_key] == species_name
        ]
        e_form_list = [
            atoms.info[e_form_key] 
            if e_index is None else atoms.info[e_form_key][e_index]
            for atoms in atoms_list
        ]
        e_form = np.min(e_form_list) if len(e_form_list) > 0 else e_form_missing
        e_form_dict[species_name] = e_form
    return e_form_dict

# -------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------