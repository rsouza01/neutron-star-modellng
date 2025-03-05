#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k, pi, sigma
from scipy.integrate import quad

stefan_boltzmann_radiance = lambda temperature:  sigma * temperature**4.

def planck_law_wavelength(wavelength, temperature):
    """
    Calculates the spectral radiance of a blackbody at a given wavelength and temperature.

    Args:
        wavelength (float or numpy.ndarray): Wavelength in meters.
        temperature (float): Temperature in Kelvin.

    Returns:
        float or numpy.ndarray: Spectral radiance in W/sr/m^3.
    """
    a = 2 * h * c**2 / wavelength**5
    b = h * c / (wavelength * k * temperature)
    intensity = a / (np.exp(b) - 1)
    return intensity

def plot_blackbody_radiation(temperatures, log_scale=False, wavelengths=None):
    """
    Plots the blackbody radiation curves for the given temperatures.

    Args:
        temperatures (list or numpy.ndarray): List of temperatures in Kelvin.
        wavelengths (numpy.ndarray, optional): Wavelength range in meters. Defaults to a reasonable range.
    """
    # create string with all elements of temperatures
    temperatures_str = '_'.join(map(str, temperatures))

    if wavelengths is None:
        wavelengths = np.linspace(1e-9, 3e-6, 1000)  # Wavelength range from 1 nm to 3 µm

    plt.figure(figsize=(10, 6))
    for temp in temperatures:
        intensity = integrate_planck_wavelength(wavelengths, temp)
        plt.plot(wavelengths * 1e9, intensity, label=f'{temp} K') # plot in nm.

    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Spectral Radiance (W/sr/m$^3$)')
    plt.title('Blackbody Radiation')
    plt.legend()
    plt.grid(True)
    if log_scale:
        plt.yscale('log') # logarithmic y scale is very useful.
        plt.xscale('log') # logarithmic x scale is also very useful.
    plt.tight_layout() #prevents labels from being cut off.
    # plt.show()
    plt.savefig(f"../output/{temperatures_str}.png")

def integrate_planck_wavelength(temperature):
    """
    Integrates Planck's law over all wavelengths at a given temperature.

    Args:
        temperature (float): Temperature in Kelvin.

    Returns:
        float: Integrated spectral radiance (total power radiated per unit area) in W/m^2.
    """
    result, error = quad(planck_law_wavelength, 1e-9, np.inf, args=(temperature,)) #start from 1nm to avoid singularity at 0
    integrated_radiance = pi * result #multiply by pi to account for the integral over all solid angles.
    return integrated_radiance

def verify_stefan_boltzmann_wavelength(temperature: float) -> None:
    """
    Verifies the Stefan-Boltzmann law using the integrated Planck's law (wavelength).

    Args:
        temperature (float): Temperature in Kelvin.
    """
    integrated_radiance = integrate_planck_wavelength(temperature)
    sb_radiance = stefan_boltzmann_radiance(temperature)

    print(f"Temperature: {temperature} K")
    print(f"Integrated Planck Radiance (wavelength): {integrated_radiance:.8e} W/m^2")
    print(f"Stefan-Boltzmann Radiance: {sb_radiance:.8e} W/m^2")
    print(f">>> Difference: {abs(integrated_radiance - sb_radiance):.2e} W/m^2")

# # Example usage:
# temperatures = [3000, 4000, 5000, 6000]  # Example temperatures in Kelvin
# plot_blackbody_radiation(temperatures, log_scale=False)

# #Example using a custom wavelength range.
# custom_wavelengths = np.linspace(100e-9, 10e-6, 1000) # 100nm to 10um
# plot_blackbody_radiation([2500, 5500, 8000], log_scale=True, wavelengths=custom_wavelengths)


# # Example using a custom wavelength range.
# custom_wavelengths = np.linspace(100e-9, 10e-6, 1000) # 100nm to 10um
# plot_blackbody_radiation([2000, 6000, 12000], log_scale=True, wavelengths=custom_wavelengths)



# Example usage:
print(80*"=")
temperature3 = 4 #very cold temperature.
verify_stefan_boltzmann_wavelength(temperature3)

print(80*"=")
temperature = 3000  # Example temperature in Kelvin
verify_stefan_boltzmann_wavelength(temperature)

print(80*"=")
temperature2 = 5778 #temperature of the sun.
verify_stefan_boltzmann_wavelength(temperature2)

