# -------------------------------------------------------------------------------------
# IMPORTS
# -------------------------------------------------------------------------------------

import timeit
import cantera as ct
import numpy as np
from ase_cantera_microkinetics import units
from ase_cantera_microkinetics.cantera_utils import (
    get_Y_dict,
    advance_sim_to_steady_state,
)
from ase_cantera_microkinetics.reaction_mechanism_from_df import (
    get_gas_species_from_df_coeffs_NASA,
    get_ads_species_from_df_fixed_T, 
    get_ts_species_from_df_fixed_T, 
    get_surf_reactions_from_df_fixed_T,
)

# -------------------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------------------

def main():

    # Set the name of the structure to investigate.
    structure = 'Ga'
    
    # Set temperature and pressure of the simulation.
    temperature_celsius = 300. # [C]
    temperature = units.Celsius_to_Kelvin(temperature_celsius) # [K]
    pressure = 2 * units.mega*units.Pa # [Pa]
    
    # Set molar fractions of gas species.
    gas_molfracs_inlet = {
        'CO2': 0.24,
        'H2': 0.72,
        'H2O': 0.00,
        'CH3OH': 0.00,
        'CH2O': 0.00,
        'CO': 0.00,
        'N2': 0.04,
    }
    
    # Set the initial coverages of the free catalytic sites.
    cat_coverages_inlet = {'(Zr)': 8/18, '(Ga)': 1/18, '(O)': 9/18}

    # Set number of cstr for the discretization of the pfr.
    n_cstr = 100

    # Set the catalyst parameters.
    cat_mass = 0.10 * units.gram # [kg]
    cat_mol_weight = 123.22 * units.gram/units.mole # [kg/kmol]
    cat_density = 6.10 * units.gram / units.centimeter**3 # [kg/m^3]
    cat_site_density = 2 * 1.41e-9 * units.mole/units.centimeter**2 # [kmol/m^2]
    area_spec = 42 * units.meter**2/units.gram # [m^2/kg]

    # Calculate catalyst dispersion.
    cat_dispersion = area_spec*cat_site_density*cat_mol_weight

    # Set the volumetric flow rate from gas hour space velocity.
    gas_hour_space_vel = 24000 * 1/units.hour # [1/s]
    vol_flow_rate = gas_hour_space_vel*(cat_mass/cat_density) # [m^3/s]
    
    # Set the cross section and the total volume of the reactor.
    reactor_length = 1.00 * units.centimeter # [m]
    diameter = 1.00 * units.millimeter # [m]

    # Calculate reactor volume and gas velocity.
    cross_section = np.pi*(diameter**2)/4. # [m^2]
    reactor_volume = cross_section*reactor_length # [m^3]
    gas_velocity = vol_flow_rate/cross_section # [m/s]

    # Calculate the catalyst active area per unit of reactor volume.
    cat_moles = cat_mass/cat_mol_weight # [kmol]
    cat_area = cat_moles*cat_dispersion/cat_site_density # [m^2]
    alpha_cat = cat_area/reactor_volume # [m^2/m^3]

    # Reaction mechanism parameters.
    excel_file = 'reaction_mechanism_ZrO2.xlsx'
    energy_ref_funs = {
        'H' : lambda energy: energy['H2']/2.,
        'O' : lambda energy: energy['H2O']-2*energy['H'],
        'C' : lambda energy: energy['CO2']-2*energy['O'],
        'N' : lambda energy: energy['N2']/2.,
    }

    # Create a solution representing the gas phase.
    gas_species = get_gas_species_from_df_coeffs_NASA(
        filename = excel_file,
        energy_ref_funs = energy_ref_funs,
        names = gas_molfracs_inlet.keys(),
    )
    gas = ct.Solution(
        name = 'gas',
        species = gas_species,
        thermo = 'ideal-gas',
        transport = 'Mix',
    )
    gas.TPX = temperature, pressure, gas_molfracs_inlet

    # Create interface representing the catalyst surface.
    ads_species = get_ads_species_from_df_fixed_T(
        filename = excel_file,
        structure = structure,
        temperature_fixed = temperature,
    )
    surf_reactions = get_surf_reactions_from_df_fixed_T(
        filename = excel_file,
        gas = gas,
        site_density = cat_site_density,
        structure = structure,
        temperature_fixed = temperature,
    )
    cat = ct.Interface(
        name = 'cat',
        species = ads_species,
        reactions = surf_reactions,
        thermo = 'ideal-surface',
        kinetics = 'surface',
        adjacent = [gas],
    )
    cat.TP = temperature, pressure
    cat.site_density = cat_site_density
    cat.coverages = cat_coverages_inlet
    
    # Create interface representing the transition states.
    ts_species = get_ts_species_from_df_fixed_T(
        filename = excel_file,
        structure = structure,
        temperature_fixed = temperature,
    )
    cat_ts = ct.Interface(
        name = 'cat_ts',
        species = ts_species,
        thermo = 'ideal-surface',
        kinetics = 'surface',
        adjacent = [gas],
    )
    cat_ts.TP = temperature, pressure
    cat_ts.site_density = cat_site_density

    # Create an ideal reactor.
    cstr_length = reactor_length/n_cstr
    cstr_volume = cross_section*cstr_length
    cstr_cat_area = alpha_cat*cstr_volume # [m^2]
    mass_flow_rate = vol_flow_rate*gas.density # [kg/s]
    cstr = ct.IdealGasReactor(gas, energy = 'off', name = 'cstr')
    cstr.volume = cstr_volume
    surf = ct.ReactorSurface(cat, cstr, A = cstr_cat_area)
    
    # Add mass flow and pressure controllers.
    upstream = ct.Reservoir(gas, name = 'upstream')
    master = ct.MassFlowController(
        upstream = upstream,
        downstream = cstr,
        mdot = mass_flow_rate,
    )
    downstream = ct.Reservoir(gas, name = 'downstream')
    pcontrol = ct.PressureController(
        upstream = cstr,
        downstream = downstream,
        master = master,
        K = 1e-6,
    )

    # Define the parameters of the simulation.
    sim = ct.ReactorNet([cstr])
    sim.rtol = 1.0e-12
    sim.atol = 1.0e-18
    sim.max_steps = 1e9
    sim.max_err_test_fails = 1e9
    sim.max_time_step = 1.0

    # Integrate the plug-flow reactor along the z (reactor length) axis.
    gas_array = ct.SolutionArray(gas, extra = ['z_reactor'])
    cat_array = ct.SolutionArray(cat, extra = ['z_reactor'])
    for ii in range(n_cstr+1):
        z_reactor = ii*cstr_length
        gas_array.append(state = gas.state, z_reactor = z_reactor)
        cat_array.append(state = cat.state, z_reactor = z_reactor)
        # Intergate the microkinetic model.
        if ii < n_cstr:
            advance_sim_to_steady_state(sim = sim, n_try_max = 1000)
        # Set the gas composition of the inlet of the next cstr equal to the 
        # outlet of the previous cstr.
        gas.TDY = cstr.thermo.TDY
        upstream.syncState()

    # Calculate the conversion of CO2 and the productivity of CH3OH.
    massfracs_in = get_Y_dict(gas_array[0])
    massfracs_out = get_Y_dict(gas_array[-1])
    conversion_CO2 = (massfracs_in['CO2']-massfracs_out['CO2'])/massfracs_in['CO2']
    productiv_CH3OH = massfracs_out['CH3OH']*mass_flow_rate/cat_mass
    productiv_CH3OH /= units.milligram/units.gram/units.hour
    print('')
    print(f'Conversion of CO2     = {conversion_CO2*100:+12.6f} [%]')
    print(f'Productivity of CH3OH = {productiv_CH3OH:+12.6f} [mg/g_cat/hour]')

# -------------------------------------------------------------------------------------
# IF NAME MAIN
# -------------------------------------------------------------------------------------

if __name__ == '__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print(f'\nExecution time = {stop-start:6.3} s\n')

# -------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------
