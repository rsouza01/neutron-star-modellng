#!/usr/bin/env python3

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import json
import math
import scipy.constants as phys_const

import eos as eos

OUTPUT_FOLDER="../output"

# Define the system of differential equations
def system(r, y):
    """
    Defines the system of two first-order differential equations.

    Args:
        t: The independent variable (time).
        y: A NumPy array containing the two dependent variables (y[0] and y[1]).

    Returns:
        A NumPy array containing the derivatives of the two dependent variables (dy[0]/dt and dy[1]/dt).
    """
 
    # y1 = y[0] # dPdr
    # y2 = y[1] # dMdr
    # Example system: 
    # dy1/dt = -y1 + y2 
    # dy2/dt = y1 - 2*y2
    # dydt = [-y1 + y2, y1 - 2*y2]  # Replace with your actual equations

    P = y[0]
    m = y[1]

    dPdr = -m * phys_const.G * rho(P)/(r**2. + 1e-20)
    dmdr = 4. * math.PI * r**2.0 * rho(P)

    dydt = [
        dPdr,
        dmdr
    ]
    return dydt

def rho(p):
    """
    Energy Density of a neutron star at a given pressure.

    Args:
        p: Pressure at a given value of r.

    Returns:
        The energy density at the given pressure.
    """
    # return 1 # constant
    return eos.energy_density(p)


# Event function for stop condition
def event(r, y):
    """
    Defines the event that triggers the integration to stop.

    Args:
        r: The current time.
        y: The current values of the dependent variables.

    Returns:
        A value that, when crosses zero, triggers the event.
    """
    P = y[0]

    # Example stop condition: Stop when P(t) reaches 0
    return P  # Replace with your actual stop condition


def plot_results(r_span, sol):
    """
    Plots the results of the integration.

    Args:
        sol: The solution returned by solve_ivp.
    """
    # Generate points for plotting (up to the event, if it occurred)
    r_eval = np.linspace(r_span[0], sol.t[-1] if sol.status == 1 else r_span[1], 200)
    y = sol.sol(r_eval)

    # Plot the results
    plt.figure(figsize=(8, 6))
    plt.plot(r_eval, y[0], label='P(r)')
    plt.plot(r_eval, y[1], label='m(r)')
    # add a vertical line to plt at y=5
    
    plt.xlabel('Radius (r)')
    plt.ylabel('Solutions (P, m)')
    plt.title('Stellar toy model')
    plt.legend()
    plt.grid(True)
    # plt.show()
    plt.savefig(f"{OUTPUT_FOLDER}/rk4.png")


# add main function
if __name__ == "__main__":
    event.terminal = True  # Important: This tells solve_ivp to stop at the event
    event.direction = -1  # Optional:  1 for crossing from negative to positive, -1 for positive to negative, 0 for either.

    P_0 = 1
    M_0 = 0

    # Initial conditions
    y0 = [P_0, M_0]

    # Time span for the solution (a large enough span to potentially reach the event)
    r_span = (0, 2000)  # radius span

    # Solve the system of differential equations with the event
    sol = solve_ivp(system, r_span, y0, events=event, dense_output=True)

    # Check if the event occurred
    if sol.status == 1:
        print(f"Event occurred at t = {sol.t[-1]}")
    else:
        print("Event did not occur within the given time span.")

    result = {
        "radius": sol.t[-1], 
        "pressure": sol.y[0][-1],
        "mass": sol.y[1][-1]
        }

    '''Printing Inital Values'''
    print(f"Radius: {result['radius']}")
    print(f"Mass: {result['mass']}")

 
    plot_results(r_span, sol)

with open(f"{OUTPUT_FOLDER}/result.json", "w") as f:
    json.dump(result, f, indent=4)    

