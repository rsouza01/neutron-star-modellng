import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import json
import math
import scipy.constants as phys_const


GAMMA = 3./2.
K = 1

mu = 1 # Pure hidrogen
mH = 1.67e-27 # Mass of hidrogen
k_B = 1.38e-23 # Boltzmann constant


def gradients(m, P, r):
    dPdr = -m * phys_const.G * rho(P)/(r**2. + 1e-20)
    dmdr = 4. * math.pi * r**2.0 * rho(P)
    return [dPdr, dmdr]

def rho(p):
    """
    Energy Density of a neutron star at a given pressure.

    Args:
        p: Pressure at a given value of r.

    Returns:
        The energy density at the given pressure.
    """
    # return 1 # constant
    return energy_density(p)


def energy_density(pressure):
    return ideal_gas_eos(pressure, GAMMA)
    # return polytropic_eos(pressure, GAMMA)

def ideal_gas_eos(pressure, temperature):
    return (mu * mH * pressure) / (k_B * temperature)


def polytropic_eos(pressure, gamma):
    if (pressure < 0):
        return 0
    
    return math.pow(pressure / K, 1./gamma)
