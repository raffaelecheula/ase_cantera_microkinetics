# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

import cantera as ct
import numpy as np
from . import units

# -----------------------------------------------------------------------------
# GET X DICT
# -----------------------------------------------------------------------------

def get_X_dict(phase):

    X_dict = {}
    for ii in range(phase.n_species):
        name = phase.species_names[ii]
        X_dict[name] = phase.X[ii]

    return X_dict

# -----------------------------------------------------------------------------
# GET Y DICT
# -----------------------------------------------------------------------------

def get_Y_dict(phase):

    Y_dict = {}
    for ii in range(phase.n_species):
        name = phase.species_names[ii]
        Y_dict[name] = phase.Y[ii]

    return Y_dict

# -----------------------------------------------------------------------------
# GET ENTHALPIES DICT
# -----------------------------------------------------------------------------

def get_enthalpies_dict(phase, units_energy = units.J/units.kmol):

    h_dict = {}
    for ii in range(phase.n_species):
        name = phase.species_names[ii]
        h_dict[name] = (
            phase.standard_enthalpies_RT[ii]*ct.gas_constant*phase.T
        )
        h_dict[name] /= units_energy

    return h_dict

# -----------------------------------------------------------------------------
# GET STD ENTROPIES DICT
# -----------------------------------------------------------------------------

def get_std_entropies_dict(phase, units_energy = units.J/units.kmol):

    s0_dict = {}
    # standard_entropies_R are calculated at phase.P, not at reference P.
    if isinstance(phase, ct.Interface):
        log_P = 0.
    else:
        log_P = np.log(phase.P/ct.one_atm)
    for ii in range(phase.n_species):
        name = phase.species_names[ii]
        s0_dict[name] = (
            (phase.standard_entropies_R[ii]+log_P)*ct.gas_constant
        )
        s0_dict[name] /= units_energy

    return s0_dict

# -----------------------------------------------------------------------------
# GET STD GIBBS DICT
# -----------------------------------------------------------------------------

def get_std_gibbs_dict(phase, units_energy = units.J/units.kmol):

    g0_dict = {}
    # standard_gibbs_RT are calculated at phase.P, not at reference P.
    if isinstance(phase, ct.Interface):
        log_P = 0.
    else:
        log_P = np.log(phase.P/ct.one_atm)
    for ii in range(phase.n_species):
        name = phase.species_names[ii]
        g0_dict[name] = (
            (phase.standard_gibbs_RT[ii]-log_P)*ct.gas_constant*phase.T
        )
        g0_dict[name] /= units_energy

    return g0_dict

# -----------------------------------------------------------------------------
# GET ENTROPIES DICT
# -----------------------------------------------------------------------------

def get_entropies_dict(phase, units_energy = units.J/units.kmol):

    s0_dict = {}
    for ii in range(phase.n_species):
        name = phase.species_names[ii]
        if isinstance(phase, ct.Interface):
            xi = phase[name].coverages[0]
        else:
            xi = phase[name].X[0]
        log_xi = np.log(xi) if xi > 0. else -np.inf
        s0_dict[name] = (
            (phase.standard_entropies_R[ii]+log_xi)*ct.gas_constant
        )
        s0_dict[name] /= units_energy

    return s0_dict

# -----------------------------------------------------------------------------
# GET GIBBS DICT
# -----------------------------------------------------------------------------

def get_gibbs_dict(phase, units_energy = units.J/units.kmol):

    g0_dict = {}
    for ii in range(phase.n_species):
        name = phase.species_names[ii]
        if isinstance(phase, ct.Interface):
            xi = phase[name].coverages[0]
        else:
            xi = phase[name].X[0]
        log_xi = np.log(xi) if xi > 0. else -np.inf
        g0_dict[name] = (
            (phase.standard_gibbs_RT[ii]+log_xi)*ct.gas_constant*phase.T
        )
        g0_dict[name] /= units_energy

    return g0_dict

# -----------------------------------------------------------------------------
# REACTIONS FROM CAT TS
# -----------------------------------------------------------------------------

def reactions_from_cat_ts(gas, cat, cat_ts, h0_dict = None, s0_dict = None):

    if h0_dict is None:
        h0_dict = get_enthalpies_dict(gas)
        h0_dict.update(get_enthalpies_dict(cat))
        h0_dict.update(get_enthalpies_dict(cat_ts))
    if s0_dict is None:
        s0_dict = get_std_entropies_dict(gas)
        s0_dict.update(get_std_entropies_dict(cat))
        s0_dict.update(get_std_entropies_dict(cat_ts))

    for ii in [ii for ii in range(cat.n_reactions)]:

        reaction = cat.reaction(ii)
        name_ts = cat_ts.species_names[ii]

        reactants_gas = []
        reactants_ads = []
        for name in reaction.reactants:
            if name in gas.species_names:
                reactants_gas += [name]*int(reaction.reactants[name])
            elif name in cat.species_names:
                reactants_ads += [name]*int(reaction.reactants[name])

        h0_act = h0_dict[name_ts]
        s0_act = s0_dict[name_ts]
        for name in reaction.reactants:
            h0_act -= h0_dict[name]*reaction.reactants[name]
            s0_act -= s0_dict[name]*reaction.reactants[name]

        if isinstance(cat.reaction(ii).rate, ct.InterfaceArrheniusRate):

            A_pre = units.kB/units.hP*np.exp(s0_act/units.Rgas)
            A_pre *= (units.Rgas*cat.T/ct.one_atm)**(len(reactants_gas))
            A_pre *= cat.site_density**(1-len(reactants_ads))
            b_temp = 1.0

            reaction.rate = ct.InterfaceArrheniusRate(
                A = A_pre,
                b = b_temp,
                Ea = h0_act,
            )

        #elif isinstance(cat.reaction(ii).rate, ct.StickingArrheniusRate):
        #
        #    gas_spec = reactants_gas[0]
        #    mol_weight = gas[gas_spec].molecular_weights[0]
        #    A_pre = np.sqrt(units.Rgas/(2*np.pi*mol_weight))
        #    A_pre *= cat.site_density**(-len(reactants_ads))
        #    b_temp = 0.5
        #    k_for = A_pre*cat.T**b_temp
        #    s0_act = -units.Rgas*np.log(k_for/(units.kB*cat.T/units.hP))
        #    
        #    reaction.rate = ct.StickingArrheniusRate(Ea = h0_act)

        cat.modify_reaction(ii, reaction)

# -----------------------------------------------------------------------------
# GET DELTA MOLAR FLUXES
# -----------------------------------------------------------------------------

def get_delta_molar_fluxes(gas, sim, TDY, mdot, upstream):

    gas.TDY = TDY
    if upstream is not None:
        upstream.syncState()

    n_dot_zero = mdot/gas.mean_molecular_weight
    ni_dot_zero = gas.X*n_dot_zero

    advance_sim_to_steady_state(sim=sim, n_try_max=1000)

    n_dot = mdot/gas.mean_molecular_weight
    ni_dot = gas.X*n_dot
    
    dni_dot = ni_dot-ni_dot_zero
    dni_dot = np.array([x if abs(x) > 0. else 1e-20 for x in dni_dot])

    return dni_dot

# -----------------------------------------------------------------------------
# MODIFY ENTHALPY
# -----------------------------------------------------------------------------

def modify_enthalpy(
    spec,
    coeffs_zero = None,
    e_form = None,
    delta_e = None,
    units_energy = units.eV/units.molecule,
):

    if coeffs_zero is not None:
        coeffs = coeffs_zero.copy()
    else:
        coeffs = spec.thermo.coeffs.copy()

    if isinstance(spec.thermo, ct.ConstantCp):
        if e_form is not None:
            coeffs[1] = e_form*units_energy
        elif delta_e is not None:
            coeffs[1] += delta_e*units_energy
        spec.thermo = ct.ConstantCp(
            T_low = spec.thermo.min_temp,
            T_high = spec.thermo.max_temp,
            P_ref = spec.thermo.reference_pressure,
            coeffs = coeffs,
        )
    
    elif isinstance(spec.thermo, ct.NasaPoly2):
        if e_form is not None:
            coeffs[13] += e_form*units_energy/ct.gas_constant-coeffs[6]
            coeffs[6] = e_form*units_energy/ct.gas_constant
        elif delta_e is not None:
            coeffs[13] += delta_e*units_energy/ct.gas_constant
            coeffs[6] += delta_e*units_energy/ct.gas_constant
        spec.thermo = ct.NasaPoly2(
            T_low = spec.thermo.min_temp,
            T_high = spec.thermo.max_temp,
            P_ref = spec.thermo.reference_pressure,
            coeffs = coeffs,
        )

    elif isinstance(spec.thermo, ct.ShomatePoly2):
        if e_form is not None:
            coeffs[13] += e_form*units_energy/(units.kiloJoule/units.mole)-coeffs[6]
            coeffs[6] = e_form*units_energy/(units.kiloJoule/units.mole)
        elif delta_e is not None:
            coeffs[13] += delta_e*units_energy/(units.kiloJoule/units.mole)
            coeffs[6] += delta_e*units_energy/(units.kiloJoule/units.mole)
        spec.thermo = ct.ShomatePoly2(
            T_low = spec.thermo.min_temp,
            T_high = spec.thermo.max_temp,
            P_ref = spec.thermo.reference_pressure,
            coeffs = coeffs,
        )

    else:
        raise NotImplementedError("thermo class not implemented.")

    return spec

# -----------------------------------------------------------------------------
# ADVANCE SIM TO STEADY STATE
# -----------------------------------------------------------------------------

def advance_sim_to_steady_state(sim, n_try_max = 1000, try_except = True):

    if try_except is True:
        n_try = 0
        finished = False
        max_time_step = sim.max_time_step
        while finished is False and n_try < n_try_max:
            try:
                n_try += 1
                sim.reinitialize()
                sim.advance_to_steady_state()
                finished = True
            except KeyboardInterrupt:
                raise
            except:
                if sim.max_time_step == 0.0:
                    sim.max_time_step = 1.0
                else:
                    sim.max_time_step *= 0.1
        if finished is False:
            raise RuntimeError("Calculation Failed!")
        sim.max_time_step = max_time_step
    else:
        sim.reinitialize()
        sim.advance_to_steady_state()

# -----------------------------------------------------------------------------
# DEGREE RATE CONTROL
# -----------------------------------------------------------------------------

def degree_rate_control(
    gas,
    cat,
    sim,
    mdot,
    upstream,
    gas_spec_target,
    multip_value = 1.05,
    return_dict = False,
):

    TDY = gas.TDY
    index_gas_spec = gas.species_names.index(gas_spec_target)

    dni_dot_orig = get_delta_molar_fluxes(
        gas = gas,
        sim = sim,
        TDY = TDY,
        mdot = mdot,
        upstream = upstream,
    )

    dni_dot_orig = get_delta_molar_fluxes(
        gas = gas,
        sim = sim,
        TDY = TDY,
        mdot = mdot,
        upstream = upstream,
    )
    dni_dot_orig = get_delta_molar_fluxes(
        gas = gas,
        sim = sim,
        TDY = TDY,
        mdot = mdot,
        upstream = upstream,
    )

    DRC_vect = []
    DRC_dict = {}
    for ii in range(cat.n_reactions):
        cat.set_multiplier(value = multip_value, i_reaction = ii)
        dni_dot_mod = get_delta_molar_fluxes(
            gas = gas,
            sim = sim,
            TDY = TDY,
            mdot = mdot,
            upstream = upstream,
        )
        cat.set_multiplier(value = 1.0, i_reaction = ii)
        DRC_species = (
            (dni_dot_mod-dni_dot_orig)/dni_dot_orig/(multip_value-1.)
        )
        DRC_vect.append(DRC_species[index_gas_spec])
        if 'name' in cat.reaction(ii).input_data:
            name = cat.reaction(ii).input_data['name']
        else:
            name = cat.reaction(ii).equation
        DRC_dict[name] = DRC_species[index_gas_spec]

    if return_dict:
        return DRC_dict
    else:
        return DRC_vect

# -----------------------------------------------------------------------------
# DEGREE RATE CONTROL
# -----------------------------------------------------------------------------

def generalized_degree_rate_control(
    gas,
    cat,
    cat_ts,
    sim,
    mdot,
    upstream,
    gas_spec_target,
    delta_e = 0.001,
    units_energy = units.eV/units.molecule,
    reactions_DRC = True,
    adsorbates_DRC = True,
    exclude_species = [],
    return_dict = False,
):

    TDY = gas.TDY
    index_gas_spec = gas.species_names.index(gas_spec_target)
    reactions_from_cat_ts(
        gas = gas,
        cat = cat,
        cat_ts = cat_ts,
    )

    dni_dot_orig = get_delta_molar_fluxes(
        gas = gas,
        sim = sim,
        TDY = TDY,
        mdot = mdot,
        upstream = upstream,
    )

    DRC_vect = []
    DRC_dict = {}
    
    if reactions_DRC is True:
        for ii in range(cat_ts.n_species):
            spec = cat_ts.species(ii)
            if spec.name in exclude_species:
                continue
            coeffs_zero = spec.thermo.coeffs.copy()
            spec = modify_enthalpy(
                spec = spec,
                coeffs_zero = coeffs_zero,
                delta_e = delta_e,
                units_energy = units_energy,
            )
            cat_ts.modify_species(ii, spec)
            reactions_from_cat_ts(
                gas = gas,
                cat = cat,
                cat_ts = cat_ts,
            )
            dni_dot_mod = get_delta_molar_fluxes(
                gas = gas,
                sim = sim,
                TDY = TDY,
                mdot = mdot,
                upstream = upstream,
            )
            spec = modify_enthalpy(
                spec = spec,
                coeffs_zero = coeffs_zero,
                delta_e = 0.,
                units_energy = units_energy,
            )
            cat_ts.modify_species(ii, spec)
            reactions_from_cat_ts(
                gas = gas,
                cat = cat,
                cat_ts = cat_ts,
            )
            DRC_species = (
                -(dni_dot_mod-dni_dot_orig)/dni_dot_orig / 
                (delta_e*(units.eV/units.molecule)/units.Rgas/cat.T)
            )
            DRC_vect.append(DRC_species[index_gas_spec])
            DRC_dict[spec.name] = DRC_species[index_gas_spec]

    if adsorbates_DRC is True:
        for ii in range(cat.n_species):
            spec = cat.species(ii)
            if spec.name in exclude_species:
                continue
            coeffs_zero = spec.thermo.coeffs.copy()
            spec = modify_enthalpy(
                spec = spec,
                coeffs_zero = coeffs_zero,
                delta_e = -delta_e,
                units_energy = units_energy,
            )
            cat.modify_species(ii, spec)
            reactions_from_cat_ts(
                gas = gas,
                cat = cat,
                cat_ts = cat_ts,
            )
            dni_dot_mod = get_delta_molar_fluxes(
                gas = gas,
                sim = sim,
                TDY = TDY,
                mdot = mdot,
                upstream = upstream,
            )
            spec = modify_enthalpy(
                spec = spec,
                coeffs_zero = coeffs_zero,
                delta_e  = 0.,
                units_energy = units_energy,
            )
            cat.modify_species(ii, spec)
            reactions_from_cat_ts(
                gas = gas,
                cat = cat,
                cat_ts = cat_ts,
            )
            DRC_species = (
                -(dni_dot_mod-dni_dot_orig)/dni_dot_orig / 
                (delta_e*(units.eV/units.molecule)/units.Rgas/cat.T)
            )
            DRC_vect.append(DRC_species[index_gas_spec])
            DRC_dict[spec.name] = DRC_species[index_gas_spec]

    if return_dict:
        return DRC_dict
    else:
        return DRC_vect

# -----------------------------------------------------------------------------
# REACTION PATH ANALYSIS
# -----------------------------------------------------------------------------

def reaction_path_analysis(
    gas, cat, surf, filename = None, mode = 'w+', perc_thold = 0.,
):

    RPA_text = ""

    names_dict = {}
    r_net_dict = {}
    r_for_dict = {}
    r_rev_dict = {}
    phi_r_dict = {}
    num_r_dict = {}
    signs_dict = {}

    for ii in range(gas.n_reactions):
        react = gas.reaction(ii)
        r_net = gas.net_rates_of_progress[ii]
        r_for = gas.forward_rates_of_progress[ii]
        r_rev = gas.reverse_rates_of_progress[ii]
        if r_for < r_rev:
            r_for, r_rev = r_rev, r_for
            signs_dict[react] = -1
        else:
            signs_dict[react] = +1
        num_r_dict[react] = ii
        if 'name' in gas.reaction(ii).input_data:
            names_dict[react] = gas.reaction(ii).input_data['name']
        else:
            names_dict[react] = gas.reaction(ii).equation
        r_net_dict[react] = r_net
        r_for_dict[react] = r_for
        r_rev_dict[react] = r_rev
        phi_r_dict[react] = r_for/(r_for+r_rev+1e-50)

    for ii in range(cat.n_reactions):
        react = cat.reaction(ii)
        r_net = cat.net_rates_of_progress[ii]*surf.area
        r_for = cat.forward_rates_of_progress[ii]*surf.area
        r_rev = cat.reverse_rates_of_progress[ii]*surf.area
        if r_for < r_rev:
            r_for, r_rev = r_rev, r_for
            signs_dict[react] = -1
        else:
            signs_dict[react] = +1
        num_r_dict[react] = ii
        if 'name' in cat.reaction(ii).input_data:
            names_dict[react] = cat.reaction(ii).input_data['name']
        else:
            names_dict[react] = cat.reaction(ii).equation
        r_net_dict[react] = r_net
        r_for_dict[react] = r_for
        r_rev_dict[react] = r_rev
        phi_r_dict[react] = r_for/(r_for+r_rev+1e-50)

    rates_dict = {}
    for spec in gas.species_names:
        rates_dict[spec] = {}
    for spec in cat.species_names:
        rates_dict[spec] = {}

    for react in r_net_dict:
        for spec in react.reactants:
            rates_dict[spec][react] = -react.reactants[spec]*r_net_dict[react]
        for spec in react.products:
            rates_dict[spec][react] = +react.products[spec]*r_net_dict[react]

    for spec in rates_dict:

        prod_sum = sum([
            rates_dict[spec][react] for react in rates_dict[spec]
            if rates_dict[spec][react] > 0.
        ])
        cons_sum = sum([
            rates_dict[spec][react] for react in rates_dict[spec]
            if rates_dict[spec][react] < 0.
        ])

        RPA_block = '\n'+'-'*140+'\n'
        RPA_block += spec.ljust(72+3+6)
        RPA_block += 'r_net'.ljust(11+3)
        RPA_block += 'r_perc'.ljust(7+3)
        RPA_block += 'phi'.ljust(6+3)
        RPA_block += 'r_for'.ljust(11+3)
        RPA_block += 'r_rev'.ljust(11)
        RPA_block += '\n'+'-'*140+'\n'

        for react in sorted(
            rates_dict[spec], key = lambda react: rates_dict[spec][react]
        ):

            if rates_dict[spec][react] >= 0.:
                prod_perc = rates_dict[spec][react]/(prod_sum+1e-50)*100
            else:
                prod_perc = -rates_dict[spec][react]/(cons_sum+1e-50)*100

            if abs(prod_perc) > perc_thold:

                num_r = num_r_dict[react]*signs_dict[react]
                name_r = names_dict[react]
                if signs_dict[react] == -1:
                    if ' <=> ' in name_r:
                        name_r = ' <=> '.join(reversed(name_r.split(' <=> ')))
                    elif ' => ' in name_r:
                        name_r = ' <= '.join(reversed(name_r.split(' => ')))
                    elif ' = ' in name_r:
                        name_r = ' = '.join(reversed(name_r.split(' = ')))

                RPA_block += f' {num_r:+4d} '
                RPA_block += f'{name_r:72s}   '
                RPA_block += f'{rates_dict[spec][react]:+11.4e}   '
                RPA_block += f'{prod_perc:+7.2f}   '
                RPA_block += f'{phi_r_dict[react]:6.4f}   '
                RPA_block += f'{r_for_dict[react]:+11.4e}   '
                RPA_block += f'{r_rev_dict[react]:+11.4e}'
                RPA_block += '\n'

        RPA_text += RPA_block

    if filename is not None:
        with open(filename, mode) as fileobj:
            fileobj.write(RPA_text)

    return RPA_text

# -----------------------------------------------------------------------------
# GET NONEQUILIBRIUM RATIO
# -----------------------------------------------------------------------------

def get_nonequilirium_ratio(
    molfracs_dict,
    molfracs_eq_dict,
    reactants,
    products,
):
    
    eta_react = 1.
    for spec in products:
        eta_react *= molfracs_dict[spec]/molfracs_eq_dict[spec]
    for spec in reactants:
        eta_react *= molfracs_eq_dict[spec]/molfracs_dict[spec]

    return eta_react

# -----------------------------------------------------------------------------
# GET MOLFRACS FROM LAMBDA
# -----------------------------------------------------------------------------

def get_molfracs_from_lambda(
    molfracs_zero_dict,
    lambda_list,
    reactants_dict,
    products_dict,
):
    
    molfracs_dict = molfracs_zero_dict.copy()
    for spec in molfracs_dict:
        for ii, reaction in enumerate(reactants_dict):
            for spec_react in [
                spec_react for spec_react in reactants_dict[reaction]
                if spec_react == spec
            ]:
                molfracs_dict[spec] -= lambda_list[ii]
            for spec_prod in [
                spec_prod for spec_prod in products_dict[reaction]
                if spec_prod == spec
            ]:
                molfracs_dict[spec] += lambda_list[ii]

    molfracs_tot = sum([molfracs_dict[spec] for spec in molfracs_dict])
    for spec in molfracs_dict:
        molfracs_dict[spec] /= molfracs_tot

    return molfracs_dict

# -----------------------------------------------------------------------------
# ERROR ETA DICT
# -----------------------------------------------------------------------------

def error_eta_dict(
    x,
    eta_dict,
    molfracs_zero_dict,
    molfracs_eq_dict,
    reactants_dict,
    products_dict,
):
    
    molfracs_dict = get_molfracs_from_lambda(
        molfracs_zero_dict = molfracs_zero_dict,
        lambda_list = x,
        reactants_dict = reactants_dict,
        products_dict = products_dict,
    )

    error = []
    for react_name in eta_dict:
        reactants = reactants_dict[react_name]
        products = products_dict[react_name]
        eta_react = get_nonequilirium_ratio(
            molfracs_dict = molfracs_dict,
            molfracs_eq_dict = molfracs_eq_dict,
            reactants = reactants,
            products = products,
        )
        error.append(eta_react-eta_dict[react_name])

    for spec in molfracs_dict:
        if molfracs_dict[spec] < 0.:
            error += (0.-molfracs_dict[spec])*1e-3
        if molfracs_dict[spec] > 1.:
            error += (molfracs_dict[spec]-1.)*1e-3

    return error

# -----------------------------------------------------------------------------
# GET MOLFRACS FROM ETA DICT
# -----------------------------------------------------------------------------

def get_molfracs_from_eta_dict(
    gas,
    eta_dict,
    molfracs_zero_dict,
    x0 = None,
    method = 'hybr',
    options = None,
):

    from .reaction_mechanism import NameAnalyzer
    from scipy.optimize import root

    gas.TPX = gas.T, gas.P, molfracs_zero_dict
    molfracs_zero_dict = get_X_dict(gas)
    gas.equilibrate('TP')
    
    molfracs_eq_dict = {}
    for spec in gas.species_names:
        molfracs_eq_dict[spec] = gas[spec].X[0]
    
    name_analyzer = NameAnalyzer()
    reactants_dict = {}
    products_dict = {}
    for react_name in eta_dict:
        reactants_dict[react_name] = name_analyzer.get_reactants(react_name)
        products_dict[react_name] = name_analyzer.get_products(react_name)

    if x0 is None:
        x0 = [0.10]*len(eta_dict)
    
    results = root(
        fun = error_eta_dict,
        x0 = x0,
        method = method,
        options = options,
        args = (
            eta_dict,
            molfracs_zero_dict,
            molfracs_eq_dict,
            reactants_dict,
            products_dict,
        ),
    )
    if results.success is False:
        raise RuntimeError("root calculation not converged!")
    
    lambda_list = results.x
    
    molfracs_dict = get_molfracs_from_lambda(
        molfracs_zero_dict = molfracs_zero_dict,
        lambda_list = lambda_list,
        reactants_dict = reactants_dict,
        products_dict = products_dict,
    )

    gas.TPX = gas.T, gas.P, molfracs_dict

    return molfracs_dict

# -----------------------------------------------------------------------------
# GET DELTA G REACTION
# -----------------------------------------------------------------------------

def get_deltaG_reaction(
    reaction_name,
    gibbs_dict,
):
    
    from .reaction_mechanism import NameAnalyzer
    
    name_analyzer = NameAnalyzer()
    reactants = name_analyzer.get_reactants(reaction_name)
    products = name_analyzer.get_products(reaction_name)
    deltaG = 0.
    for spec in reactants:
        deltaG -= gibbs_dict[spec]
    for spec in products:
        deltaG += gibbs_dict[spec]

    return deltaG

# -----------------------------------------------------------------------------
# END
# -----------------------------------------------------------------------------
