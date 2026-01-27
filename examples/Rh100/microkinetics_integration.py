# -------------------------------------------------------------------------------------
# IMPORTS
# -------------------------------------------------------------------------------------

import yaml
import timeit
import cantera as ct
import numpy as np
from tqdm import tqdm
from ase.db import connect
from ase_cantera_microkinetics import units
from ase_cantera_microkinetics.cantera_utilities import (
    get_Y_dict,
    reactions_from_cat_ts,
    advance_sim_to_steady_state,
    degree_rate_control,
    molar_balance_of_element,
)
from ase_cantera_microkinetics.ase_utilities import (
    get_mechanism_from_yaml,
    get_atoms_list_from_db,
)

# -------------------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------------------

def main():

    # Parameters.
    ensemble = True
    calculate_DRC = True
    e_form_zero_yaml = "e_form_zero.yaml"
    deltae_form_yaml = "deltae_form.yaml"
    mechanism_yaml = "mechanism.yaml"
    results_yaml = "results.yaml" if ensemble else None
    verbose = False
    n_calculations = 2000
    index_list = list(range(n_calculations)) if ensemble else [0]
    
    # Read input yaml files.
    with open(e_form_zero_yaml, "r") as fileobj:
        e_form_dict_zero = yaml.safe_load(fileobj)
    with open(deltae_form_yaml, "r") as fileobj:
        deltae_dict = yaml.safe_load(fileobj)
    # Run calculations.
    results_list = []
    for index in index_list if verbose else tqdm(index_list, ncols=100):
        # Get formation energies.
        if ensemble is True:
            e_form_dict = {
                name: e_form_dict_zero[name] + deltae_dict[name][index]
                for name in e_form_dict_zero
            }
        else:
            e_form_dict = e_form_dict_zero.copy()
        # Integrate the microkinetic model.
        results = integrate_microkinetic_model(
            e_form_dict=e_form_dict,
            mechanism_yaml=mechanism_yaml,
            calculate_DRC=calculate_DRC,
            verbose=verbose,
        )
        results_list.append(results)
    # Save results to yaml.
    if ensemble is True:
        results_all = {
            key: [results[key] for results in results_list]
            for key in results_list[0].keys()
        }
        save_results_to_yaml(results=results_all, results_yaml=results_yaml)

# -------------------------------------------------------------------------------------
# INTEGRATE MICROKINETIC MODEL
# -------------------------------------------------------------------------------------

def integrate_microkinetic_model(
    e_form_dict: dict,
    mechanism_yaml: str,
    calculate_DRC: bool = False,
    elements_molar_balance: bool = False,
    verbose: bool = True,
) -> dict:
    """
    Integrate the microkinetic model.
    """
    # Set temperature and pressure of the simulation.
    temperature_celsius = 500. # [C]
    temperature = units.Celsius_to_Kelvin(temperature_celsius) # [K]
    pressure = 1 * units.atm # [Pa]
    # Set molar fractions of gas species.
    gas_molfracs_inlet = {
        "CO2": 0.28,
        "H2": 0.28,
        "H2O": 0.02,
        "CO": 0.02,
        "N2": 0.40,
    }
    # Set the initial coverages of the free catalytic sites.
    cat_coverages_inlet = {"(Rh)": 1.0}
    # Select species for the calculation of TOF and conversion.
    gas_reac = "CO2"
    gas_prod = "CO"
    mol_weight_prod = 28 # [kg/kmol]
    # Set number of cstr for the discretization of the PFR.
    n_cstr = 1
    # Set the catalyst and reactor parameters.
    alpha_cat = 20000.0 # [m^2/m^3]
    cat_site_density = 1e-8 # [kmol/m^2]
    vol_flow_rate = 1e-6 # [m^3/s]
    # Set the cross section and the total volume of the reactor.
    reactor_length = 1.00 * units.centimeter # [m]
    diameter = 1.00 * units.millimeter # [m]
    cross_section = np.pi * (diameter ** 2) / 4 # [m^2]
    reactor_volume = cross_section * reactor_length # [m^3]
    cat_sites_tot = cat_site_density * alpha_cat * reactor_volume # [kmol]
    # Get reaction mechanism.
    gas, cat, cat_ts = get_mechanism_from_yaml(
        yaml_file=mechanism_yaml,
        e_form_dict=e_form_dict,
        site_density=cat_site_density,
        temperature=temperature,
        units_energy=units.eV / units.molecule,
    )
    gas.TPX = temperature, pressure, gas_molfracs_inlet
    cat.TP = temperature, pressure
    cat.coverages = cat_coverages_inlet
    cat_ts.TP = temperature, pressure
    # Calculate gas space velocity.
    molecular_weight = gas.mean_molecular_weight # [kg/kmol]
    molar_flow_rate = vol_flow_rate * gas.density / molecular_weight # [kmol/s]
    cat_moles = cat_site_density * alpha_cat * reactor_volume # [kmol]
    gas_space_velocity = molar_flow_rate / cat_moles # [1/s]
    # Update the reactions from the transition states.
    reactions_from_cat_ts(gas=gas, cat=cat, cat_ts=cat_ts)
    # Create an ideal reactor.
    cstr_length = reactor_length / n_cstr # [m]
    cstr_volume = cross_section * cstr_length # [m^3]
    cstr_cat_area = alpha_cat * cstr_volume # [m^2]
    mass_flow_rate = vol_flow_rate * gas.density # [kg/s]
    cstr = ct.IdealGasReactor(gas, energy="off", name="cstr")
    cstr.volume = cstr_volume
    surf = ct.ReactorSurface(cat, cstr, A=cstr_cat_area)
    # Add mass flow and pressure controllers.
    upstream = ct.Reservoir(gas, name="upstream")
    master = ct.MassFlowController(
        upstream=upstream,
        downstream=cstr,
        mdot=mass_flow_rate,
    )
    downstream = ct.Reservoir(gas, name="downstream")
    pcontrol = ct.PressureController(
        upstream=cstr,
        downstream=downstream,
        master=master,
        K=1e-6,
    )
    # Define the parameters of the simulation.
    sim = ct.ReactorNet([cstr])
    sim.rtol = 1.0e-15
    sim.atol = 1.0e-18
    sim.max_steps = 1e9
    sim.max_err_test_fails = 1e9
    sim.max_time_step = 1.0
    # Print the reaction mechanism energies.
    if verbose is True:
        print("\n- Reaction mechanism energies.")
        for name, e_form in e_form_dict.items():
            print(f"{name}: {e_form:+6.3f} [eV]")
    # Integrate the plug-flow reactor along the z (reactor length) axis.
    gas_array = ct.SolutionArray(gas, extra=["z_reactor"])
    cat_array = ct.SolutionArray(cat, extra=["z_reactor"])
    # Print the header of the table.
    if verbose is True:
        print("\n- Reactor integration.")
        string = "distance[m]".rjust(14)
        for spec in gas.species_names:
            string += ("x_"+spec+"[-]").rjust(12)
        print(string)
    # Integrate the reactor.
    for ii in range(n_cstr+1):
        z_reactor = ii * cstr_length
        gas_array.append(state=gas.state, z_reactor=z_reactor)
        cat_array.append(state=cat.state, z_reactor=z_reactor)
        # Print the state of the reactor.
        string = f"  {z_reactor:12f}"
        for spec in gas.species_names:
            string += f"  {gas[spec].X[0]:10f}"
        if verbose is True:
            print(string)
        # Intergate the microkinetic model.
        if ii < n_cstr:
            advance_sim_to_steady_state(sim=sim, n_try_max=1000)
        # Set the gas composition of the inlet of the next cstr equal to the 
        # outlet of the previous cstr.
        gas.TDY = cstr.thermo.TDY
        upstream.syncState()
    # Calculate the conversion and the TOF.
    massfracs_in = get_Y_dict(gas_array[0])
    massfracs_out = get_Y_dict(gas_array[-1])
    conversion_reac = (
        (massfracs_in[gas_reac] - massfracs_out[gas_reac]) / massfracs_in[gas_reac]
    )
    mass_flow_prod = (
        (massfracs_out[gas_prod] - massfracs_in[gas_prod]) * mass_flow_rate # [kg/s]
    )
    tof_prod = mass_flow_prod / mol_weight_prod / cat_sites_tot # [1/s]
    if verbose is True:
        print("\n- Catalytic activity.")
        print(f"Conversion of {gas_reac} = {conversion_reac * 100:+12.6f} [%]")
        print(f"Turnover frequency of {gas_prod} = {tof_prod:+12.6e} [1/s]")
    # Get the coverages.
    coverages = {spec: float(cat[spec].coverages[0]) for spec in cat.species_names}
    if verbose is True:
        print("\n- Coverages.")
        for spec, coverage in coverages.items():
            print(f"{spec:30s} {coverage:7.4e}")
    # Calculate elements molar balance.
    if elements_molar_balance is True:
        delta_moles_fracts = molar_balance_of_element(
            gas_in=gas_array[0],
            gas_out=gas_array[-1],
            mass=mass_flow_rate,
            return_fraction=True,
        )
        if verbose is True:
            print("\n- Molar balances of elements.")
            for elem, delta_molfract in delta_moles_fracts.items():
                print(f"Balance of {elem} = {delta_molfract * 100:+7.2f} [%]")
    # Calculate the degree of rate control.
    if calculate_DRC is True:
        DRC_dict = degree_rate_control(
            gas=gas,
            cat=cat,
            sim=sim,
            mdot=mass_flow_rate,
            upstream=upstream,
            gas_spec_target=gas_prod,
            multip_value=1.05,
            return_dict=True,
        )
        # Print DRC results.
        if verbose is True:
            print("\n- Degree of rate control.")
            for ii, react in enumerate(DRC_dict):
                print(f"{react:50s} {DRC_dict[react]:+7.4f}")
            sum_DRC = np.sum([DRC_dict[react] for react in DRC_dict])
            print(f"Sum DRC = {sum_DRC:+7.4f}")
    else:
        DRC_dict = {}
    # Return results.
    results = {
        "TOF": float(tof_prod),
        "conversion": float(conversion_reac),
        **{f"theta {key}": float(value) for key, value in coverages.items()},
        **{f"DRC {key}": float(value) for key, value in DRC_dict.items()},
    }
    return results

# -------------------------------------------------------------------------------------
# SAVE RESULTS TO YAML
# -------------------------------------------------------------------------------------

def save_results_to_yaml(
    results: dict,
    results_yaml: str = "results.yaml",
) -> None:
    """
    Save the results to a YAML file.
    """
    # Custom YAML dumper.
    class CustomDumper(yaml.SafeDumper):
        pass
    def float_representer(dumper, value):
        return dumper.represent_scalar("tag:yaml.org,2002:float", f"{value:+10.8E}")
    CustomDumper.add_representer(float, float_representer)
    # Write the results.
    with open(results_yaml, "w") as fileobj:
        yaml.dump(
            data=results,
            stream=fileobj,
            default_flow_style=None,
            width=1000,
            sort_keys=False,
            Dumper=CustomDumper,
        )

# -------------------------------------------------------------------------------------
# IF NAME MAIN
# -------------------------------------------------------------------------------------

if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print(f"\nExecution time = {stop - start:6.3f} [s]\n")

# -------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------
