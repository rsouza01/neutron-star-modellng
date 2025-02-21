#!/usr/bin/env python3


# Import necessary Libraries to get started
import numpy as np               # we will need a lot of arrays in this project
import matplotlib.pyplot as plt  # for visualising the results

def initial_n():

    '''Initial number of neutrons per unit volume at r = 0'''

    n = 1                                   # Initial value of n
    err = 1                                 # Inital Error
    tol = 1e-15                             # Tolerance Value
    count = 0                               # To track total iterations for convergence

    # Newton-Raphson method
    while err > tol :                       # Untill convergence of n
        count += 1                          # Increment the count of the iteration
        fn = 236*n**(2.54) + n*mn - rho_s   # Find f(n_old)
        dfn = 236*2.54*n**(1.54) + mn       # Find f'(n_old)
        temp = n - fn/dfn                   # Calculate n_new using Newton Raphson Method
        err = np.abs(n-temp)                # Calculate the error
        n = temp                            # Use this n value in next iteration for f(n_old), f'(n_old)
    print("Newton-Raphson Converged after ", count, "iterations")
    return n

def rho(p):

    '''
    Energy Density of a neutron star at a given pressure
    rho_s - Central Density at r = 0
    mn - mass of a neutron
    n - number density at given pressure
    '''

    n = (p*rho_s/363.44)**(1/2.54)
    return (236. * n**2.54 + n *mn)/rho_s


def dp_dr(r,m,p,flag):

    '''
    Pressure Gradient
    r - distance from the center of neutron star
    m - mass enclosed within the value of r
    p - pressure at that value of r
    flag - choose between classical and relativistic models
    '''

    if flag == 0:
        # Classical Model
        y = -m*rho(p)/(r**2 + 1e-20)

    else:
        # Relativistic Model
        y = -(p+rho(p))*(m + p*r**3)/(r**2 - 2*m*r + 1e-20)

    return y

def dm_dr(r,m,p):

    '''
    Mass Gradient
    Same r, m, p as we used for pressure gradient are to be used for mass gradient
    '''

    return rho(p)*r**2

def RK4Solver(r,m,p,h,flag):

    '''
    Runge - Kutta 4th order Method
    Calculates 4 different pressure and mass gradients
    K1_ - Mass Gradients
    K2_ - Pressure Gradients
    '''
    y = np.zeros(2)                                               # To store the next value of m and p

    # Gradient 1 at start point
    k11 = dm_dr(r, m, p)
    k21 = dp_dr(r, m, p, flag)

    # Gradient 2 at mid of start and end point
    k12 = dm_dr(r + 0.5*h, m + 0.5*k11*h, p + 0.5*k21*h)
    k22 = dp_dr(r + 0.5*h, m + 0.5*k11*h, p + 0.5*k21*h, flag)

    # Gradient 3 at mid of start and end point
    k13 = dm_dr(r + 0.5*h, m + 0.5*k12*h, p + 0.5*k22*h)
    k23 = dp_dr(r + 0.5*h, m + 0.5*k12*h, p + 0.5*k22*h, flag)

    # Gradient 4 at end point
    k14 = dm_dr(r + h, m + k13*h, p + k23*h)
    k24 = dp_dr(r + h, m + k13*h, p + k23*h, flag)

    # Updation of m and p
    y[0] = m + h*(k11 + 2*k12 + 2*k13 + k14)/6
    y[1] = p + h*(k21 + 2*k22 + 2*k23 + k24)/6

    return y


def plot_data(color, label):

    # Mass Profile
    plt.subplot(1, 2, 1)
    plt.plot(r*R0*1e-18, m*M0/Ms, color = color, linewidth = 1.2, label = label)
    plt.xlabel('Distance, $r$ (km)', fontsize = 13)
    plt.ylabel('Mass, $M/M_{sun}$', fontsize = 13)
    plt.title('Mass Profile of a Neutron Star', color = 'tab:red', weight = 'bold', fontsize = 15)
    plt.xlim(left = 0)
    plt.ylim(bottom = 0)
    plt.legend(fontsize = 13, frameon = False)

    # Pressure Profile
    plt.subplot(1, 2, 2)
    plt.plot(r*R0*1e-18, p*rho_s, color = color, linewidth = 1.2, label = label)
    plt.xlabel('Distance, $r$ (km)', fontsize = 13)
    plt.ylabel('Pressure, $P$ $(MeV/fm^{3})$', fontsize = 13)
    plt.title('Pressure Profile of a Neutron Star', color = 'tab:red', weight = 'bold', fontsize = 15)
    plt.xlim(left = 0)
    plt.ylim(bottom = 0)
    plt.legend(fontsize = 13, frameon = False)


# Simulation parameters and Conversion factors
hc = 197.327                                # conversion factor in MeV fm (h bar * c)
G = hc * 6.67259e-45                        # gravitational constant in MeV^-1 fm^3 kg^-1
Ms = 1.1157467e60                           # mass of Sun in MeV
rho_s = 1665.3                              # central density of neutron star (density at r = 0) in MeV/fm^3
M0 = (4*3.14159265*(G**3)*rho_s)**(-0.5)
R0 = G*M0
mn = 938.926                                # mass of a neutron in MeV c^-2


#Initialising Values and Arrays

'''Values for radius, step size and tolerance value'''
N = 1501                              # Total number of data points
r = np.linspace(0,15,N)               # values of radius to compute enclosed mass within it, pressure and density profiles
h = r[1]-r[0]                         # step size for RK4Solver
tol = 9e-5                            # tolerance for RK4Solver

'''Arrays to store the updation values in RK4Solver'''
m = np.zeros(N)                       # mass
p = np.zeros(N)                       # pressure

'''Initial Number Density of Neutrons at r = 0'''
ni = initial_n()

'''Setting intial values in the array - values at the center of the neutron star'''
r[0] = 0                              # radius
m[0] = 0                              # mass
p[0] = 363.44 * (ni**2.54)/rho_s      # pressure

'''Setting flags for choosing between classical and relativistic model'''
flag_set = [0,1]

'''Printing Inital Values'''
print("Initial number density, ni =", ni)
print("Initial Pressure, P[0] = ", p[0]*rho_s, "MeV/fm^3")
print("Simulation range, R = 0 to R =", r[-1]*R0*1e-18, "km") # 1e-18 is the conversion factor from femtometers to kilometers
print("Step size for RK4 Solver:", h)


'''Using RK4 Numerical Method for modeling a neutron star'''
plt.figure(figsize = (18, 5))                                   # setting up the figure size for subplots
for k in range(0,2):                                            # k will take two values in this run --> 0 and 1
    flag = flag_set[k]                                          # flag_set[0] = 0, flag_set[1] = 1
    for i in range(0,N-1):                                      # for 1500 iterations
        if flag == 0:                                           # if the flag is set to 0
            [m[i+1], p[i+1]] = RK4Solver(r[i],m[i],p[i],h,flag) # use classical model
        else:                                                   # if not 0
            [m[i+1], p[i+1]] = RK4Solver(r[i],m[i],p[i],h,flag) # use relativistic model
        if p[i+1] < tol:                                        # if the value of the pressure goes below tolerance (close to 0 as we know)
            break                                               # the simulation is completed
    print()                                                     # for a new line between outputs
    if i == N-2:                                                # if the pressure value was not able to go below the tolerance value
        lbl1 = "Program didn't converge to P = 0, extend the maximum value of r" # print to increase the itetarations of simulation --> values of r
    else:                                                       # if simulation stopped before this
        lbl1 = f"P < {tol} found after {i} runs"                # print the conclusion statement

    '''Keep only the used indices of array discarding the remaining ones'''
    m = m[:i+2] # m[0:i+2 --> Exclusive index]
    p = p[:i+2]
    r = r[:i+2]

    '''Time to visualise and print the results'''
    if flag == 0:
        lbl = "Classical Model"
        plot_data('tab:orange', "Classical Model")
    else:
        lbl = "Relativistic Model"
        plot_data('tab:cyan', "Relativistic Model")

    '''Printing the overall output'''
    print ("==================================================================")
    print (lbl, "Results:", lbl1)
    print ("==================================================================")
    print ("Initial density, rho_s = ", rho_s, "MeV/fm^3")
    print ("Total mass = ", m[-1]*M0/Ms, "times Solar mass")
    print ("Radius of the Neutron star = ", r[-1]*R0*1e-18, "km")

print()
plt.subplots_adjust(wspace = 0.15)
plt.savefig('nsm.png')

