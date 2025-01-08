# -------------------------------------------------------------------------------------
# IMPORTS
# -------------------------------------------------------------------------------------

import ase.units

# -------------------------------------------------------------------------------------
# PREFIX
# -------------------------------------------------------------------------------------

exa = 1e18
peta = 1e15
tera = 1e12
giga = 1e09
mega = 1e06
kilo = 1e03
hecto = 1e02
deca = 1e01
deci = 1e-01
centi = 1e-02
milli = 1e-03
micro = 1e-06
nano = 1e-09
pico = 1e-12
femto = 1e-15
atto = 1e-18

# -------------------------------------------------------------------------------------
# BASIC UNITS
# -------------------------------------------------------------------------------------

meter = mt = 1.0  # [mt]
second = sec = 1.0  # [sec]
kilogram = kg = 1.0  # [kg]
kilomole = kmol = 1.0  # [kmol]
Kelvin = K = 1.0  # [K]
Pascal = Pa = 1.0  # [kg/mt/sec^2]
Newton = N = 1.0  # [kg*mt/sec^2]
Joule = J = 1.0  # [kg*mt^2/sec^2]
Watt = W = 1.0  # [kg*mt^2/sec^3]
gram = kilogram / kilo
mol = mole = kilomole / kilo

# -------------------------------------------------------------------------------------
# DERIVED UNITS
# -------------------------------------------------------------------------------------

Angstrom = Ang = 1e-10 * meter
inch = 0.0254 * meter
litre = liter = 1e-03 * meter ** 3
minute = 60 * second
hour = 3600 * second
Hertz = Hz = 1 / second
Dalton = amu = 1.660539040e-27 * kg
atmosphere = atm = 101325 * Pascal
bar = 1e05 * Pascal
calorie = cal = 4.184 * Joule
eV = 1e03 / ase.units.kJ * Joule
molecule = 1.0 / ase.units.mol * mole

# -------------------------------------------------------------------------------------
# CONSTANTS
# -------------------------------------------------------------------------------------

k_Boltzmann = kB = ase.units.kB * eV  # [J/K]
h_Planck = hP = 6.62607e-34 * Joule * sec  # [J*sec]
N_Avogadro = Navo = 1.0 / molecule  # [1/kmol]
R_ideal_gas = Rgas = kB * Navo  # [J/kmol/K]
c_light = c_l = ase.units._c  # [m/s]

# -------------------------------------------------------------------------------------
# CONSTRUCTED UNITS
# -------------------------------------------------------------------------------------

decimeter = deci * meter
centimeter = centi * meter
millimeter = milli * meter
micrometer = micro * meter
nanometer = nano * meter

hectogram = hecto * gram
decagram = deca * gram
decigram = deci * gram
centigram = centi * gram
milligram = milli * gram
microgram = micro * gram
nanogram = nano * gram

millimole = milli * mole
micromole = micro * mole
nanomole = nano * mole

kiloPascal = kilo * Pascal
megaPascal = mega * Pascal

kiloNewton = kilo * Newton
megaNewton = mega * Newton

kiloJoule = kilo * Joule
megaJoule = mega * Joule

# -------------------------------------------------------------------------------------
# CELSIUS TO KELVIN
# -------------------------------------------------------------------------------------

def Celsius_to_Kelvin(temperature):
    return temperature + 273.15  # [K]

# -------------------------------------------------------------------------------------
# KELVIN TO CELSIUS
# -------------------------------------------------------------------------------------

def Kelvin_to_Celsius(temperature):
    return temperature - 273.15  # [C]

# -------------------------------------------------------------------------------------
# NORMAL LITER
# -------------------------------------------------------------------------------------

def NormalLiter(temperature, pressure):
    return 1e-3 * (101325 / pressure) * (temperature / 273.15)  # [m^3]

# -------------------------------------------------------------------------------------
# NORMAL CUBIC METER
# -------------------------------------------------------------------------------------

def NormalCubicMeter(temperature, pressure):
    return (101325 / pressure) * (temperature / 273.15)  # [m^3]

# -------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------
