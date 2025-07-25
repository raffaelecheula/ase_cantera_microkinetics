# -------------------------------------------------------------------------------------
# IMPORTS
# -------------------------------------------------------------------------------------

import numpy as np
from scipy.optimize import least_squares

from ase_cantera_microkinetics import units

# -------------------------------------------------------------------------------------
# READ VIB ENERGIES
# -------------------------------------------------------------------------------------

def read_vib_energies(
    filename: str = "vib.log",
    imaginary: bool = False,
) -> list:
    """Read vibrational energies from an ASE Vibrations log file."""
    vib_energies = []
    with open(filename, "rU") as fileobj:
        lines = fileobj.readlines()
    for i in range(3, len(lines)):
        if lines[i][0] == "-":
            break
        string = lines[i].split()[1]
        if string[-1] == "i":
            if imaginary is True:
                vib_energies.append(complex(0., float(string[:-1])*1e-3))
        else:
            vib_energies.append(complex(float(string)*1e-3))
    return vib_energies

# -------------------------------------------------------------------------------------
# H0 FROM NASA
# -------------------------------------------------------------------------------------

def H0_from_NASA(
    T: float,
    coeffs: np.ndarray,
) -> float:
    """
    Calculate enthalpy from NASA polynomial.
    """
    H0 = (
        coeffs[0] * T + 
        coeffs[1]/2 * T**2 + 
        coeffs[2]/3 * T**3 + 
        coeffs[3]/4 * T**4 + 
        coeffs[4]/5 * T**5 + 
        coeffs[5]
    )
    return H0

# -------------------------------------------------------------------------------------
# S0 FROM NASA
# -------------------------------------------------------------------------------------

def S0_from_NASA(
    T: float,
    coeffs: np.ndarray,
) -> float:
    """
    Calculate entropy from NASA polynomial.
    """
    S0 = (
        coeffs[0] * np.log(T) + 
        coeffs[1] * T + 
        coeffs[2]/2 * T**2 + 
        coeffs[3]/3 * T**3 + 
        coeffs[4]/4 * T**4 + 
        coeffs[6]
    )
    return S0

# -------------------------------------------------------------------------------------
# NASA COEFFS RESIDUALS FUN
# -------------------------------------------------------------------------------------

def NASA_coeffs_residuals_fun(
    coeffs: np.ndarray,
    T_array: np.ndarray,
    H0_array: np.ndarray,
    S0_array: np.ndarray,
):
    """
    Calculate residuals for NASA polynomial coefficients.
    """
    H0_fit = H0_from_NASA(T_array, coeffs)
    S0_fit = S0_from_NASA(T_array, coeffs)
    return np.append(H0_array-H0_fit, S0_array-S0_fit)

# -------------------------------------------------------------------------------------
# FIT NASA COEFFS
# -------------------------------------------------------------------------------------

def fit_NASA_coeffs(T_array, H0_array, S0_array):
    """
    Fit NASA polynomial coefficients.
    """
    opt = least_squares(
        fun  = NASA_coeffs_residuals_fun,
        x0 = np.ones(7),
        args = (T_array, H0_array, S0_array)
    )
    coeffs = opt.x
    return coeffs

# -------------------------------------------------------------------------------------
# PRINT NASA COEFFS
# -------------------------------------------------------------------------------------

def print_NASA_coeffs(
    coeffs: np.ndarray,
    name: str = None,
    ljust: int = 71,
    fileobj: object = None,
):
    """
    Print NASA polynomial coefficients.
    """
    if name is not None:
        print(f'"{name}"'.ljust(ljust), end=", ", file=fileobj)
    for coeff in coeffs:
        print(f"{coeff:+15.8E}", end=", ", file=fileobj)
    print("", file=fileobj)

# -------------------------------------------------------------------------------------
# PRINT NASA COEFFS
# -------------------------------------------------------------------------------------

def plot_NASA_coeffs(
    T_array: np.ndarray,
    H0_array: np.ndarray,
    S0_array: np.ndarray,
    coeffs: np.ndarray,
):
    """
    Plot NASA polynomial coefficients.
    """
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.subplot(211)
    plt.plot(T_array, H0_array, label="H0")
    plt.plot(T_array, H0_from_NASA(T_array, coeffs), label="H0 fit")
    plt.xlabel("temperature [K]")
    plt.ylabel("enthalpy/Rgas [1/K]")
    plt.legend()
    plt.subplot(212)
    plt.plot(T_array, S0_array, label="S0")
    plt.plot(T_array, S0_from_NASA(T_array, coeffs), label="S0 fit")
    plt.xlabel("temperature [K]")
    plt.ylabel("entropy/Rgas [-]")
    plt.legend()
    plt.show()

# -------------------------------------------------------------------------------------
# ASE THERMO TO NASA COEFFS
# -------------------------------------------------------------------------------------

def ase_thermo_to_NASA_coeffs(
    thermo: object,
    n_points: int = 100,
    t_low: float = 200,
    t_max: float = 1000,
    coeffs_ref: np.ndarray = None,
    subtract_ZPE = True,
    print_coeffs: bool = False,
    name: str = None,
    fileobj: object = None,
    plot_fit: bool = False,
):
    """
    Convert ASE thermo class to NASA polynomial coefficients.
    """
    T_array = np.zeros(n_points)
    H0_array = np.zeros(n_points)
    S0_array = np.zeros(n_points)
    # Get zero point energy.
    if hasattr(thermo, "get_zero_point_energy"):
        ZPE = thermo.get_zero_point_energy(verbose=False)
    else:
        ZPE = thermo.get_ZPE_correction()
    # Calculate enthalpy and entropy.
    for ii in range(n_points):
        T_array[ii] = t_low+ii*(t_max-t_low)/n_points
        if hasattr(thermo, "get_enthalpy"):
            H0 = thermo.get_enthalpy(
                temperature=T_array[ii],
                verbose=False,
            )
        if hasattr(thermo, "get_internal_energy"):
            H0 = thermo.get_internal_energy(
                temperature=T_array[ii],
                verbose=False,
            )
        else:
            raise RuntimeError("thermo class cannot calculate H0.")
        if subtract_ZPE is True:
            H0_array[ii] = (H0-ZPE) * units.eV/units.molecule/units.Rgas
        else:
            H0_array[ii] = H0 * units.eV/units.molecule/units.Rgas
        S0 = thermo.get_entropy(
            temperature=T_array[ii],
            verbose=False,
        )
        S0_array[ii] = S0 * units.eV/units.molecule/units.Rgas
        if coeffs_ref is not None:
            H0_array[ii] -= H0_from_NASA(T_array[ii], coeffs_ref)
            S0_array[ii] -= S0_from_NASA(T_array[ii], coeffs_ref)
    # Fit NASA polynomial coefficients.
    coeffs = fit_NASA_coeffs(
        T_array=T_array,
        H0_array=H0_array,
        S0_array=S0_array,
    )
    # Print NASA polynomial coefficients.
    if print_coeffs is True:
        print_NASA_coeffs(coeffs, name=name, fileobj=fileobj)
    # Plot NASA polynomial coefficients.
    if plot_fit is True:
        plot_NASA_coeffs(T_array, H0_array, S0_array, coeffs)
    # Return NASA polynomial coefficients.
    return coeffs

# -------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------
