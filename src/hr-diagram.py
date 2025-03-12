#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k, pi, sigma
from scipy.integrate import quad
from astropy import constants as cte 

def get_cross_product_stars(radiuses, temperatures):
    stars = []
    for radius in radiuses:
        for temperature in temperatures:

            luminosity = 4 * pi * (radius ** 2) * sigma * (temperature ** 4)
            stars.append(star(radius, temperature, luminosity))
    return stars

def plot_hr_diagram(stars):
    plt.figure(figsize=(10, 6))

    plt.xlabel("Effective surface temperature (T_eff/K)")
    plt.ylabel("Luminosity (L/L_sun)")
    plt.title("HR Diagram")
    plt.legend()
    plt.grid(True)
    # show more values in x axis
    plt.locator_params(axis="x", nbins=100)
    # invert x axis in plt
    plt.gca().invert_xaxis()
    # extract x and y values from stars, x = temperature, y = luminosity
    x = [star["temperature"] for star in stars]
    y = [star["luminosity"] for star in stars]

    plt.yscale("log") # logarithmic y scale is very useful.
    plt.xscale("log") # logarithmic x scale is also very useful.

    plt.scatter(x, y)

    plt.tight_layout() #prevents labels from being cut off.
    plt.savefig("../output/hr_diagram.png")

def star(radius, temperature, luminosity):
    return {"radius": radius, "temperature": temperature, "luminosity": luminosity}

# generate main
if __name__ == '__main__':

    star_radiuses_in_R_sun = [0.1, 1, 10, 100]
    star_radiuses_in_km = [radius * cte.R_sun.value for radius in star_radiuses_in_R_sun]

    effective_surface_temperatures_kelvin = [2000, 4000, 6000, 10000, 20000, 40000]
    stars = get_cross_product_stars(star_radiuses_in_km, effective_surface_temperatures_kelvin)

    plot_hr_diagram(stars)