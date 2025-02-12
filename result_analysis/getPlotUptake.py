
from ase.io import read
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os


# Constants
R = 8.314  # Gas constant in J/(mol·K)
N = 6.022e+23  # molecules per mole
BAR2PASCAL = 1e5


# Peng-Robinson EOS parameters for CO2
def calculate_pr_eos_constants(T):

    a = 0.45724 * (R**2 * T_c**2) / P_c
    b = 0.07780 * (R * T_c) / P_c
    alpha = (1 + (0.37464 + 1.54226 * omega - 0.26992 * omega**2) * (1 - np.sqrt(T / T_c)))**2
    return a, b, alpha


# Function to calculate molar volume V_m using the Peng-Robinson EOS
def calculate_molar_volume(P, T):
    a, b, alpha = calculate_pr_eos_constants(T)
    A = a * alpha * P / (R**2 * T**2)
    B = b * P / (R * T)

    # Coefficients of the cubic equation (Z^3 - Z^2 + Z - ...)
    coeffs = [1, -(1 - B), A - 2*B - 3*B**2, -(A*B - B**2 - B**3)]

    # Solve for compressibility factor Z (real roots only)
    Z_roots = np.roots(coeffs)
    Z = max(np.real(Z_roots[np.isreal(Z_roots)]))  # Use largest real root for gas phase

    # Calculate molar volume Vm = Z * R * T / P
    V_m = Z * R * T / P
    return V_m


def calcRhoIdealGas(pressure, temperature):
    """

    Parameters:
    - pressure (float): Pressure (in Pa).
    - temperature (float): Temperature (in K).
    - mol_weight (float): Molar mass of the gas (in kg/mol).

    Returns:
    - rho_bulk_gas (float): (in g/m^3).
    """
    return pressure / (R * temperature)


# Function to calculate density of CO2 using molar volume
def calcRho(P, T):
    V_m = calculate_molar_volume(P, T)
    rho_g = M_CO2 / V_m  # density in kg/m^3
    return rho_g


def getAvgExcess(abs_avg_nAds, void_volume, rho_bulk_gas):

    # Calculate number of molecules of CO2
    moles_co2 = rho_bulk_gas * void_volume
    n_molecules_co2 = moles_co2 * N
    print("# excess CO2", n_molecules_co2)

    avg_nAdsEx = abs_avg_nAds - n_molecules_co2
    return avg_nAdsEx


def calcRhoWithCoolProp(pressure, temperature, fluid="CO2"):
    import CoolProp.CoolProp as CP
    return CP.PropsSI("D", "T", temperature, "P", pressure, fluid)


def getAvgExcessWithCoolProp(abs_avg_nAds,pressure, temperature, volume, fluid="CO2"):

    import CoolProp.CoolProp as CP

    molar_mass = CP.PropsSI("M", fluid)
    density = CP.PropsSI("D", "T", temperature, "P", pressure, fluid)

    # Calculate the number of moles
    n_moles = (density * volume) / molar_mass
    n_molecules = n_moles * N

    return abs_avg_nAds - n_molecules


def getUptake(file_path, plot=False):

    if file_path.endswith(".exyxyz"):
        # from gas loaded exyxyz file
        fl_base = file_path.split("/")[-1].replace(".extxyz", "")
        atoms_list = read(extxyz_path, index=":")
        n_atoms_frame = len(atoms_frame)
        nAds_list = [(len(atoms)-n_atoms_frame)/3 for atoms in atoms_list]

    elif file_path.endswith(".csv"):
        # from status file
        df = pd.read_csv(csv_path)
        succ_insert = df[" succ_insertion"].tolist()
        succ_del = df["succ_deletion"].tolist()
        nAds_list = [n_insert-n_del for n_insert, n_del in zip(succ_insert,succ_del)]

    elif file_path.endswith(".npy"):
        # from npy file
        fl_base = results_dir
        nAds_list = np.load(file_path).tolist()

    abs_avg_nAds = np.array(nAds_list[int(len(nAds_list)/2):]).mean()
    #  excess_avg_nAds = getAvgExcess(abs_avg_nAds, void_volume, rho_bulk_gas)
    excess_avg_nAds = getAvgExcessWithCoolProp(abs_avg_nAds, pressure, temperature, void_volume, fluid="CO2")

    if plot:
        plt.plot(np.array(range(len(nAds_list))), nAds_list)
        plt.xlabel(r"Steps")
        plt.ylabel(r"Number of Molecules")
        plt.savefig(f"{fl_base}.png")
        plt.clf()
        #  plt.show()
    abs_uptake = (abs_avg_nAds)/sum(atoms_frame.get_masses())* 1000 # in mmol/g
    excess_uptake = (excess_avg_nAds)/sum(atoms_frame.get_masses())* 1000 # in mmol/g

    return abs_uptake, excess_uptake



# CO2 critical properties
T_c = 304.2  # K
P_c = 7.38e6  # Pa
omega = 0.225  # Acentric factor for CO2
M_CO2 = 0.04401  # kg/mol for CO2
#  M_CO2 = 44.01  # g/mol for CO2


# NOTE volume problem
#  mol_weight_co2 = 0.04401 # M_co2 in kg/mol
void_volume = 10* 2.41603e-27 # in m^3


results_dir_list = [it for it in os.listdir("./") if os.path.isdir(it) and "results" in it]
atoms_frame = read("./frame0.extxyz")


#  atoms_frame = read("./MgMOF74_clean_frame0.extxyz")
fl = open("uptakes.csv", "w")
print("FileName,Pressure(Bar),AbsUptake(mmol/g),ExcesUptake(mmol/g)", file=fl)
for results_dir in results_dir_list:
    pressure = float(results_dir.split("bar")[0].split("_")[-1]) # in bar
    temperature = float(results_dir.split("K")[0].split("_")[-1])
    print(pressure, temperature)
    #  extxyz_path = f"{results_dir}/trajectory_{pressure}bar.extxyz"
    #  csv_path = f"{results_dir}/status.csv"
    npy_path = f"{results_dir}/uptake_{pressure}bar.npy"
    #  uptake = getUptakeTraj(extxyz_path, plot=True)
    pressure *= BAR2PASCAL # in pascal
    #  rho_bulk_gas = calcRhoIdealGas(pressure, temperature)
    #  rho_bulk_gas = calcRho(pressure, temperature)
    rho_bulk_gas = calcRhoWithCoolProp(pressure, temperature)
    abs_uptake, excess_uptake = getUptake(npy_path, plot=False)
    print(f"{results_dir},{pressure},{abs_uptake},{excess_uptake}", file=fl)

