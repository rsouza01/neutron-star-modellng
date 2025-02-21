import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import json
import math


GAMMA = 3./2.
K = 1

def energy_density(pressure):
    return polytropic_eps(pressure, GAMMA)

def polytropic_eps(pressure, gamma):
    if (pressure < 0):
        return 0
    
    return math.pow(pressure / K, 1./gamma)
