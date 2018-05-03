# FUNCTIONS TO OBTAIN THE SOLUTIONS TO THE LINEAR WAVE PROBLEM
# ============================================================

import numpy as np
import sympy as sp
import scipy as sc
from scipy import optimize

# eta(x,t)
def eta_lin(x, t):
    return a * np.cos(k * x - omega * t + psi)

# phi(x,z,t)
def phi_lin(x, z, t):
    return a * omega / k * f_lin(z) * np.sin(k * x - omega * t + psi)

# u(x,z,t)
def u_lin(x, z, t):
    return a * omega * f_lin(z) * np.cos(k * x - omega * t + psi)

# w(x,z,t)
def w_lin(x, z, t):
    return a * omega * f1_lin(z) * np.sin(k * x - omega * t + psi)

# f(x,z,t)
def f_lin(z):
    return np.cosh(k * (z + h)) / np.sinh(k * h)

# f1(x,z,t)
def f1_lin(z):
    return np.sinh(k * (z + h)) / np.sinh(k * h)

# Function to obtain omega from h and k by dispersion relation
def dispersion_omega(h, k):
    return np.sqrt(g * k * np.tanh(k * h))

# Function to obtain k from omega and k by dispersion relation
def dispersion_h(omega, k):
    return np.arctanh((omega ** 2) / (g * k)) / k

# Function to obtain k from h and omega by dispersion relation
def dispersion_k(h, omega):
    def equation(k):
        return g * k * np.tanh(k * h) - (omega ** 2)

    return sc.optimize.fsolve(equation, 1)[0]

# Function to compute the parameters of the linear free surface wave problem
def compute_parameters(h=None, H=None, a=None, lamda=None,
                       k=None, T=None, omega=None, psi=0, print_results=False):
    OK = True

    # Converting inputs into floats
    try:
        h = float(h)
    except:
        pass
    try:
        H = float(H)
    except:
        pass
    try:
        a = float(a)
    except:
        pass
    try:
        lamda = float(lamda)
    except:
        pass
    try:
        k = float(k)
    except:
        pass
    try:
        T = float(T)
    except:
        pass
    try:
        omega = float(omega)
    except:
        pass
    try:
        psi = float(psi)
    except:
        pass

    # CHecking/computing wave height
    if isinstance(H, float) and a == None:
        a = H / 2
    elif isinstance(a, float) and H == None:
        H = 2 * a
    elif isinstance(H, float) and isinstance(a, float):
        if not np.isclose(H, 2 * a):
            OK = False
            if print_results:
                print('ERROR: H is not 2*a')
    else:
        OK = False
        if print_results:
            print('ERROR: Wave height (H or a) not defined')

    # CHecking/computing wave length/number
    if isinstance(lamda, float) and k == None:
        k = 2 * np.pi / lamda
    elif isinstance(k, float) and lamda == None:
        lamda = 2 * np.pi / k

    # CHecking/computing wave period/frequency
    if isinstance(T, float) and omega == None:
        omega = 2 * np.pi / T
    elif isinstance(omega, float) and T == None:
        T = 2 * np.pi / omega

    # CHecking/computing dispersion
    if isinstance(h, float) and isinstance(k, float) and omega == None:
        omega = dispersion_omega(h, k)
        T = 2 * np.pi / omega
    elif isinstance(omega, float) and isinstance(k, float) and h == None:
        h = dispersion_h(omega, k)
    elif isinstance(h, float) and isinstance(omega, float) and k == None:
        k = dispersion_k(h, omega)
        lamda = 2 * np.pi / k
    elif isinstance(h, float) and isinstance(k, float) and isinstance(omega, float):
        if not np.isclose(omega, dispersion_omega(h, k)):
            OK = False
            if print_results:
                print('ERROR: Parameters do not satisfy Dispersion Relationship')
    else:
        OK = False
        if print_results:
            print('ERROR: Not enough parameters defined')

    if OK:
        # Addimensional parameters and Froude number
        omega_adim = omega * np.sqrt(h / g)
        k_adim = k * h
        Fr = (np.sqrt(k * np.tanh(k_adim)) - omega_adim) / k_adim

    if print_results:
        # Print parameters
        print('h = ' + str(h) + ' m')
        print('H = ' + str(H) + ' m')
        print('a = ' + str(a) + ' m')
        print('lamda = ' + str(lamda) + ' m')
        print('k = ' + str(k) + ' rad/m')
        print('T = ' + str(T) + ' s')
        print('omega = ' + str(omega) + ' rad/s')
        print('psi = ' + str(psi) + ' rad')
        if OK:
            print('Fr = ' + str(Fr))
            if H / lamda < (1 / 7):
                print('H/lamda = ' + str(H / lamda) + ' < 0.143. OK')
            else:
                print('H/lamda = ' + str(H / lamda) + ' > 0.143. ¿LINEARITY?')
        else:
            print('**ERRORS WERE FOUND**')

    return h, H, a, lamda, k, T, omega, psi, OK


# P_d(x,z,t): linear dynamic pressure
def Pd_lin(x, z, t):
    return -rho * g * eta_lin(x, t) * np.cosh(k * (z + h)) / np.cosh(k * h)  # -rho*g*z
