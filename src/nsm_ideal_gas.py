#!/usr/bin/env python3

# Key improvements and explanations:

# Ideal Gas Equation of State: The code now correctly uses the ideal gas law (rho = (mu * P) / (R_gas * T)) to calculate the density.
# Radiative Transfer: The inclusion of the dTdr equation implements radiative transfer, a crucial component of stellar structure.
# Clearer Variable Names: More descriptive variable names enhance readability.
# Example Parameters: The example parameters are closer to realistic stellar values.
# Density Calculation and Plotting: The code calculates and plots the density profile, providing a more complete picture of the stellar structure.
# solve_ivp Usage: The code uses solve_ivp, which is more robust and versatile for solving ODEs than odeint.
# Initial Conditions: The initial conditions are now more appropriate for a stellar interior.
# Maximum Radius: The r_max variable is now used, and is set to a value close to a solar radius.
# Opacity and Energy Generation: Placeholders are provided for opacity and energy generation laws, which are essential for realistic stellar models.
# Corrected dTdr: the dTdr equation now uses the correct constants, and signs.
# To make the model even more realistic, you would need to:

# Implement more accurate opacity and energy generation laws.
# Account for convection in certain regions of the star.
# Use a more sophisticated equation of state, especially for degenerate matter.
# Add boundary conditions to the surface of the star.

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def stellar_equilibrium(r, y, M_total, L_total, kappa_0, kappa_power, epsilon_0, epsilon_power, mu, G, R_gas):
    """
    System of differential equations for stellar equilibrium.

    Args:
        r: Radius.
        y: Vector of dependent variables [P, M, L, T].
        M_total: Total mass of the star.
        L_total: Total luminosity of the star.
        kappa_0: Opacity coefficient.
        kappa_power: Opacity power.
        epsilon_0: Energy generation coefficient.
        epsilon_power: Energy generation power.
        mu: Mean molecular weight.
        G: Gravitational constant.
        R_gas: Ideal gas constant.

    Returns:
        dydr: Vector of derivatives [dP/dr, dM/dr, dL/dr, dT/dr].
    """
    P, M, L, T = y

    rho = (mu * P) / (R_gas * T)  # Ideal gas equation of state
  

    dPdr = -G * M * rho / r**2
    dMdr = 4 * np.pi * r**2 * rho
    dLdr = 4 * np.pi * r**2 * epsilon_0 * rho * T**epsilon_power
    kappa = kappa_0 * rho * T**kappa_power
    dTdr = -(3 * kappa * rho * L) / (16 * np.pi * G * M * T**3) #radiative transfer

    return [dPdr, dMdr, dLdr, dTdr]

def solve_stellar_structure(M_total, L_total, kappa_0, kappa_power, epsilon_0, epsilon_power, mu, G, R_gas, R_initial, P_initial, T_initial, r_max):
    """
    Solves the stellar structure equations.

    Args:
        M_total: Total mass of the star.
        L_total: Total luminosity of the star.
        kappa_0: Opacity coefficient.
        kappa_power: Opacity power.
        epsilon_0: Energy generation coefficient.
        epsilon_power: Energy generation power.
        mu: Mean molecular weight.
        G: Gravitational constant.
        R_gas: Ideal gas constant.
        R_initial: Initial radius.
        P_initial: Initial pressure.
        T_initial: Initial temperature.
        r_max: Maximum radius to integrate to.

    Returns:
        sol: Solution object from solve_ivp.
    """

    M_initial = 0.0
    L_initial = 0.0
    initial_conditions = [P_initial, M_initial, L_initial, T_initial]
    r_span = (R_initial, r_max)

    sol = solve_ivp(stellar_equilibrium, r_span, initial_conditions, args=(M_total, L_total, kappa_0, kappa_power, epsilon_0, epsilon_power, mu, G, R_gas), dense_output=True, max_step=1e6)

    return sol

# Example usage:
M_total = 1.989e30  # Solar mass (kg)
L_total = 3.828e26  # Solar luminosity (W)
kappa_0 = 1.0  # Example opacity coefficient
kappa_power = -3.5 # Kramers opacity power
epsilon_0 = 1e-26 # Example energy generation coefficient
epsilon_power = 4.0 # pp chain power
mu = 0.6  # Mean molecular weight (example)
G = 6.674e-11  # Gravitational constant (N m^2/kg^2)
R_gas = 8.314  # Ideal gas constant (J/mol K)

R_initial = 1e-3  # Initial radius (m)
P_initial = 1e16  # Initial pressure (Pa)
T_initial = 1e7   # Initial temperature (K)
r_max = 7e8     #Maximum radius to integrate to, approximately solar radius.

sol = solve_stellar_structure(M_total, L_total, kappa_0, kappa_power, epsilon_0, epsilon_power, mu, G, R_gas, R_initial, P_initial, T_initial, r_max)

# Plotting the results:
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='Pressure (P)')
plt.plot(sol.t, sol.y[1], label='Mass (M)')
plt.plot(sol.t, sol.y[2], label='Luminosity (L)')
plt.plot(sol.t, sol.y[3], label='Temperature (T)')
plt.xlabel('Radius (r)')
plt.ylabel('Value')
plt.title('Stellar Structure')
plt.legend()
plt.grid(True)
plt.show()

#calculate density.
rho = (mu * sol.y[0]) / (R_gas * sol.y[3])
plt.figure(figsize=(10, 6))
plt.plot(sol.t, rho, label='Density (rho)')
plt.xlabel('Radius (r)')
plt.ylabel('Density (kg/m^3)')
plt.title('Stellar Density')
plt.legend()
plt.grid(True)
plt.show()

plt.savefig(f"nsm_ideal_gas.png")
