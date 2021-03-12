import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
rc('text', usetex=True)
rc('figure', figsize=(11.69, 8.27))

"""
Calculate the transmission
through a defect with two
linear leads and the defect
being a one level defect with
binding t' and e'
"""

"""
The following equation gets used:
T(E,V) = 4 Tr[Gamma_L * G_cc^r * Gamma_R * G_cc^a]

with:
Gamma_L/R = Im{Sigma_L/R^a}
G_cc^a = Dagger(G_cc^r) = [E + i0 - epsilon - Sigma_L - Sigma_R]
Sigmal_L/R = t_L/R^2 * G_LL/RR^r,a
"""


"""
Given parameters

t: Binding between lead sites
tt: Binding between leads and defect
eps: Energy of leads site
eps_dfct: Energy of defect site
E: Energy (incl. infinitesimal imaginary part)
"""
t = 1
tt = 0.5
eps = 5
eps_dfct = 2.5
E = np.arange(-10, 10, 0.2) + 0.01j

"""
Necessary functions to determine
the following quantities

g_LL/RR^r,a:
Sigma_L/R^a: Self-energy
Gamma_L,R: Scattering matrices
G_cc^a,r: retarded/avanced Green function of defect
"""
def calc_g(E, eps):
    return 1/(E-eps)

def calc_Sigma(t, g):
    return t**2 * g**2

def calc_Gamma(Sigma):
    return Sigma.imag

def calc_G(E, eps_dftc, Sigma_L, Sigma_R):
    return (E - eps_dftc - Sigma_L, Sigma_R)**(-1)

def calc_Transmission(Gamma_L, G_cc_r, Gamma_R, G_cc_a):
    tmp = Gamma_L * G_cc_r * Gamma_R * G_cc_a
    return 4*np.trace(tmp)

"""
Variables containing the previously
defined quantities
"""

g_LL_r = np.array([calc_g(energy, eps) for energy in E])
g_LL_a = np.conjugate(g_LL_r)
g_RR_r = g_LL_r # Assume LL and RR are same (change later)
g_RR_a = g_LL_a # Assume LL and RR are same (change later)

Sigma_L_r = np.array([calc_Sigma(t, value) for value in g_LL_r])
Sigma_L_a = np.array([calc_Sigma(t, value) for value in g_LL_a])
Sigma_R_r = Sigma_L_r # Assume L and R are same (change later)
Sigma_R_a = Sigma_L_a # Assume L and R are same (change later) 

Gamma_L = np.array([calc_Gamma(value) for value in Sigma_L_a])
Gamma_R = np.array([calc_Gamma(value) for value in Sigma_R_a])

for n in range(1):
    print(E[n])
    print(Sigma_L_r[n])
    G_cc_r = calc_G(E[n], eps_dfct, Sigma_L_r[n], Sigma_R_r[n])
    G_cc_a = np.conjugate(G_cc_r)

# T = calc_Transmission(Gamma_L, G_cc_r, Gamma_R, G_cc_a)
