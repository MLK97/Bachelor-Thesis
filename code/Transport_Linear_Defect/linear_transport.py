import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from numpy.linalg import multi_dot
from numpy.linalg import inv

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
    return tt**2 * g

def calc_Gamma(Sigma):
    return Sigma.imag

def calc_G(E, eps_dftc, Sigma_L, Sigma_R):
    return np.divide(1, (E - eps_dftc - Sigma_L, Sigma_R))  # maybe this gives wrong result?

def calc_Transmission(Gamma_L, G_cc_r, Gamma_R, G_cc_a):
    tmp = Gamma_L * G_cc_r * Gamma_R * G_cc_a
    return 4 * tmp # eigentlich Spur

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

G_cc_r = np.array([calc_G(energy, eps_dfct, sig_L, sig_R) for energy, sig_L, sig_R in zip(E, Sigma_L_r, Sigma_R_r)])
G_cc_a = np.conjugate(G_cc_r)

T = np.array([calc_Transmission(gam_L, green_r, gam_R, green_a) for gam_L, green_r, gam_R, green_a in zip(Gamma_L, G_cc_r, Gamma_R, G_cc_a)])

# print("g_LL_r:", len(g_LL_r))
# print("Sigma_L_r:", len(Sigma_L_r))
# print("Gamma_L:", len(Gamma_L))
# print("G_cc_r:", len(G_cc_r))
# print("T:", T)

plt.plot(E.real, T.real, label="Transmission Real")

plt.grid(True)
plt.title("Transmission plot")
plt.ylabel("T")
plt.xlabel("E")
plt.legend()
plt.show()
