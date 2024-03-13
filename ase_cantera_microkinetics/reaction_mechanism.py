# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

import numpy as np
from ase.formula import Formula
from ase.data import atomic_masses, atomic_numbers
import cantera as ct
from . import units
from .cantera_utils import modify_enthalpy

# -----------------------------------------------------------------------------
# NAME ANALYZER
# -----------------------------------------------------------------------------

class NameAnalyzer():

    def __init__(
        self,
        site_separators = "(,)",
        react_separator = "<=>",
        text_separator = "_",
    ):
        
        self.site_separators = site_separators
        self.react_separator = react_separator
        self.text_separator = text_separator

    def get_names_reactants_products(self, name):

        # Split equation into names of reactants and products.
        if self.react_separator in name:
            names = name.split(f" {self.react_separator} ")
        else:
            names = [name]
    
        for name in names:
            if name[0] == " " or name[-1] == " " or "  " in name:
                raise RuntimeError(f'Error in name: {name}.')
    
        return names

    def get_names_without_mult_integers(self, name):

        # Convert multiplying integers into the corresponding pieces.
        names_new = []
        for name in self.get_names_reactants_products(name):
            pieces = name.split(" + ")
            pieces_new = []
            for piece in pieces:
                if " " in piece:
                    mm, pp = piece.split(" ")
                    pieces_new += [pp]*int(mm)
                else:
                    pieces_new += [piece]
            names_new.append(" + ".join(pieces_new))

        return names_new

    def get_n_pieces_gas_ads(self, name, index = 0):

        name_new = self.get_names_without_mult_integers(name = name)[index]
        pieces = name_new.split(" + ")

        n_pieces_gas = 0
        n_pieces_ads = 0
        for piece in pieces:
            if self.site_separators[0] in piece:
                n_pieces_ads += 1
            else:
                n_pieces_gas += 1

        return n_pieces_gas, n_pieces_ads

    def get_reactants(self, name):

        name_new = self.get_names_without_mult_integers(name = name)[0]
        reactants = name_new.split(" + ")

        return reactants

    def get_products(self, name):
    
        name_new = self.get_names_without_mult_integers(name = name)[1]
        reactants = name_new.split(" + ")

        return reactants

    def get_composition(self, name, include_site = True, index = 0):
        
        name_new = self.get_names_without_mult_integers(name = name)[index]
        composition = {}
        name_new = name_new.replace(" ", "")
        spec_list = name_new.split("+")
        for spec in spec_list:
            if self.site_separators[0] in spec:
                spec, site = spec.split(self.site_separators[0])
                site = site.split(self.site_separators[2])[0]
                site_list = site.split(self.site_separators[1])
                site_str = "".join(site_list)
            else:
                site_str = ""
            spec = spec.split(self.text_separator)[0]
            if include_site is True:
                spec = site_str+spec
            comp_spec = Formula(spec, strict=True).count()
            for elem in comp_spec:
                if elem in composition:
                    composition[elem] += comp_spec[elem]
                else:
                    composition[elem] = comp_spec[elem]
        
        return composition
        
    def get_size(self, name, index = 0):
        name_new = self.get_names_without_mult_integers(name = name)[index]
        size = (
            name_new.count(self.site_separators[0])
            +name_new.count(self.site_separators[1])
        )
        return size

    def get_composition_and_size(self, name, index = 0):
        composition = self.get_composition(name, index=index)
        size = self.get_size(name, index=index)
        return composition, size

    def get_molar_mass(self, name):
        comp_dict = self.get_composition(name)
        return sum(
            [atomic_masses[atomic_numbers[ii]]*comp_dict[ii] for ii in comp_dict]
        ) # [kg/kmol]

# -----------------------------------------------------------------------------
# GET SPEC CONSTANTCP
# -----------------------------------------------------------------------------

def get_spec_ConstantCp(
    name,
    h0,
    s0,
    cP,
    temperature,
    composition = None,
    size = None,
    units_energy = units.eV/units.molecule,
    T_low = 200.00,
    T_high = 2000.00,
    P_ref = ct.one_atm,
):
    
    spec = ct.Species(
        name = name,
        composition = composition,
        size = size,
    )
    coeffs = [
        temperature,
        h0*units_energy,
        s0*units_energy,
        cP*units_energy, 
    ]
    spec.thermo = ct.ConstantCp(
        T_low = T_low,
        T_high = T_high,
        P_ref = P_ref,
        coeffs = coeffs,
    )

    return spec

# -----------------------------------------------------------------------------
# GET SPEC NASAPOLY2
# -----------------------------------------------------------------------------

def get_spec_NasaPoly2(
    name,
    coeffs_NASA,
    e_form = None,
    composition = None,
    size = None,
    units_energy = units.eV/units.molecule,
    T_low = 200.00,
    T_mid = 1000.00,
    T_high = 2000.00,
    P_ref = ct.one_atm,
):
    
    spec = ct.Species(
        name = name,
        composition = composition,
        size = size,
    )
    if len(coeffs_NASA) == 7:
        coeffs = [T_mid]+list(coeffs_NASA)*2
    elif len(coeffs_NASA) == 14:
        coeffs = [T_mid]+list(coeffs_NASA)
    else:
        coeffs = list(coeffs_NASA)
    
    if e_form is not None:
        coeffs[6] += e_form*units_energy/ct.gas_constant
        coeffs[13] += e_form*units_energy/ct.gas_constant
    
    spec.thermo = ct.NasaPoly2(
        T_low = T_low,
        T_high = T_high,
        P_ref = P_ref,
        coeffs = coeffs,
    )

    return spec

# -----------------------------------------------------------------------------
# GET SPECIES FROM G0 DICT FIXED T
# -----------------------------------------------------------------------------

def get_species_from_g0_dict_fixed_T(
    g0_dict,
    temperature_fixed,
    units_energy = units.eV/units.molecule,
    name_analyzer = NameAnalyzer(),
    composition_dict = None,
    size_dict = None,
):

    species = []
    for name in g0_dict:

        g0 = g0_dict[name]
        if not isinstance(g0, (int, float, np.integer, np.floating)):
            g0 = np.nan
        elif np.isnan(g0):
            continue
        
        composition, size = None, None
        if name_analyzer is not None:
            composition, size = name_analyzer.get_composition_and_size(name = name)
        
        if composition_dict is not None:
            composition = composition_dict[name]
        
        if size_dict is not None:
            size = size_dict[name]

        spec = get_spec_ConstantCp(
            name = name,
            composition = composition,
            size = size,
            h0 = g0,
            s0 = 0.00,
            cP = 0.00,
            temperature = temperature_fixed,
            units_energy = units_energy,
        )
        species.append(spec)

    return species

# -----------------------------------------------------------------------------
# GET SPECIES FROM COEFFS NASA DICT
# -----------------------------------------------------------------------------

def get_species_from_coeffs_NASA_dict(
    coeffs_NASA_dict,
    e_form_dict = None,
    units_energy = units.eV/units.molecule,
    name_analyzer = NameAnalyzer(),
    composition_dict = None,
    size_dict = None,
    T_low = 200.00,
    T_mid = 1000.00,
    T_high = 2000.00,
    P_ref = ct.one_atm,
):

    species = []
    for name in coeffs_NASA_dict:

        coeffs_NASA = coeffs_NASA_dict[name]
        if not isinstance(coeffs_NASA[0], (np.integer, np.floating, int, float)):
            coeffs_NASA = [np.nan]*len(coeffs_NASA)
        elif np.isnan(coeffs_NASA[0]):
            continue

        composition, size = None, None
        if name_analyzer is not None:
            composition, size = name_analyzer.get_composition_and_size(name = name)
        
        if composition_dict is not None:
            composition = composition_dict[name]
        
        if size_dict is not None:
            size = size_dict[name]

        if e_form_dict is None:
            e_form = None
        else:
            e_form = e_form_dict[name]
            if not isinstance(e_form, (int, float, np.integer, np.floating)):
                e_form = np.nan
            elif np.isnan(e_form):
                continue

        spec = get_spec_NasaPoly2(
            name = name,
            composition = composition,
            size = size,
            e_form = e_form,
            coeffs_NASA = coeffs_NASA,
            units_energy = units_energy,
            T_low = T_low,
            T_mid = T_mid,
            T_high = T_high,
            P_ref = P_ref,
        )
        species.append(spec)

    return species

# -----------------------------------------------------------------------------
# GET SURF REACT FROM E ACT
# -----------------------------------------------------------------------------

def get_surf_react_from_e_act(
    gas,
    name,
    e_act,
    pre_exp,
    units_energy = units.eV/units.molecule,
    allow_negative_e_act = True,
    allow_sticking_e_act = True,
):

    if e_act < 0. and allow_negative_e_act is False:
        e_act = 0.

    if pre_exp == 'sticking':
        if allow_sticking_e_act is False:
            e_act = 0.
        rate = ct.StickingArrheniusRate(
            A = +1.0,
            b = +0.0,
            Ea = e_act*units_energy,
        )
    else:
        rate = ct.InterfaceArrheniusRate(
            A = pre_exp,
            b = +1.0,
            Ea = e_act*units_energy,
        )
    react = ct.Reaction(
        equation = name,
        rate = rate,
        kinetics = gas,
    )
    react.update_user_data({'name': name})

    return react

# -----------------------------------------------------------------------------
# PRE EXP FROM N REACTANTS GAS ADS
# -----------------------------------------------------------------------------

def pre_exp_from_n_reactants_gas_ads(
    n_reactants_gas,
    n_reactants_ads,
    site_density,
    s0_act,
    temperature,
    pressure_ref = ct.one_atm,
):

    pre_exp = units.kB/units.hP*np.exp(s0_act/units.Rgas)
    pre_exp *= (units.Rgas*temperature/pressure_ref)**n_reactants_gas
    pre_exp *= site_density**(1-n_reactants_ads)

    return pre_exp

# -----------------------------------------------------------------------------
# GET SURF REACTIONS FROM G0 ACT DICT FIXED T
# -----------------------------------------------------------------------------

def get_surf_reactions_from_g0_act_dict_fixed_T(
    gas,
    g0_act_dict,
    temperature_fixed,
    name_analyzer = NameAnalyzer(),
    pre_exp_dict = None,
    site_density = None,
    sticking_reactions = [],
    units_energy = units.eV/units.molecule,
    P_ref = ct.one_atm,
):

    surf_reactions = []
    for name in g0_act_dict:

        g0_act = g0_act_dict[name]
        if not isinstance(g0_act, (int, float, np.integer, np.floating)):
            g0_act = np.nan
        elif np.isnan(g0_act):
            continue

        if name_analyzer is not None:
            n_pieces_gas, n_pieces_ads = name_analyzer.get_n_pieces_gas_ads(
                name = name,
            )
            pre_exp = pre_exp_from_n_reactants_gas_ads(
                n_reactants_gas = n_pieces_gas,
                n_reactants_ads = n_pieces_ads,
                site_density = site_density,
                s0_act = 0.0,
                temperature = temperature_fixed,
                pressure_ref = P_ref,
            )
        elif pre_exp_dict is not None:
            pre_exp_dict[name]
        else:
            raise RuntimeError("name_analyzer or pre_exp_dict should be specified")

        if name in sticking_reactions:
            pre_exp = 'sticking'

        react = get_surf_react_from_e_act(
            gas = gas,
            name = name,
            e_act = g0_act,
            pre_exp = pre_exp,
            units_energy = units.eV/units.molecule,
        )
        surf_reactions.append(react)

    return surf_reactions

# -----------------------------------------------------------------------------
# CHANGE REFERENCE ENERGIES
# -----------------------------------------------------------------------------
    
def change_reference_energies(
    species,
    energy_ref_funs,
    name_analyzer = NameAnalyzer(),
    composition_dict = None,
    units_energy = units.eV/units.molecule,
):

    energy = {}
    for spec in species:
        if isinstance(spec.thermo, ct.ConstantCp):
            energy[spec.name] = spec.thermo.coeffs[1]/units_energy
        elif isinstance(spec.thermo, ct.NasaPoly2):
            energy[spec.name] = (
                spec.thermo.coeffs[6]/units_energy*ct.gas_constant
            )
        elif isinstance(spec.thermo, ct.ShomatePoly2):
            energy[spec.name] = (
                spec.thermo.coeffs[6]/(units.kiloJoule/units.mole)*ct.gas_constant
            )
        else:
            raise NotImplementedError("thermo class not implemented.")

    for name in energy_ref_funs:
        energy[name] = energy_ref_funs[name](energy = energy)
    
    for spec in species:
        if composition_dict is None:
            composition = name_analyzer.get_composition(name = spec.name)
        e_form = energy[spec.name]
        for elem in composition:
            e_form -= energy[elem]*composition[elem]
        spec = modify_enthalpy(
            spec = spec,
            e_form = e_form,
            units_energy = units_energy,
        )

    return species

# -----------------------------------------------------------------------------
# END
# -----------------------------------------------------------------------------