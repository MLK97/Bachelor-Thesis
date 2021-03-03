import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
rc('text', usetex=True)
rc('figure', figsize=(11.69, 8.27))

"""
New approach for approximating G_11
via itertion and starting at g_11

Now the following eq gets used:
G_11 = g_11*[1 - g_11 * t^2 * G_11]^-1

g_11 is defined as:
g_11 = 1/(E-epsilon) and
E = E +- i eta
"""

t = 1
epsilon = 0.5

"""
Compute g_11 for a given Energy
"""
Green = np.complex128([]) # Stores values G_11(E)
Green_analytical = np.complex128([])
E = np.arange(-10, 10, 0.2) + 0.01j # Stores Energy values (x-axis)
g_11 = [] # Stores g_11(E) values


def calculate_g(E):
    return 1/(E-epsilon)

for value in E:
    g_11 = np.append(g_11, calculate_g(value))

def calculate_Green(g_11):
    G_11 = g_11
    print("start at: ", g_11)
    G_11_prev = G_11
    G_11 = g_11*(1 - g_11 * t**2 * G_11)**(-1)
    """ For 10**(-1) slowly looses perfect shape"""
    while np.abs(G_11-G_11_prev) > 10**(-2):
        G_11_prev = G_11
        G_11 = g_11*(1 - g_11 * t**2 * G_11)**(-1)
        print("G_11 Iteration: ", G_11)
    return G_11

for value in g_11:
    Green = np.append(Green, calculate_Green(value))

Green_real = Green.real
Green_imag = Green.imag

"""
Computing the analytical equations (3.21) and (3.22)
from Zilly.
"""

def analytical_Green(E):
    if (E-epsilon) <= -2*t:
        G_11 = (E-epsilon)/(2*t**2) + 1/(2*t**2) * np.sqrt((E-epsilon)**2-4*t**2)
    if np.abs(E-epsilon) < 2*t:
        G_11 = (E-epsilon)/(2*t**2) - 1j/(2*t**2) * np.sqrt(4*t**2-(E-epsilon)**2)
    if E-epsilon >= 2*t:
        G_11 = (E-epsilon)/(2*t**2) - 1/(2*t**2) * np.sqrt((E-epsilon)**2-4*t**2)
    return G_11


for value in E:
    Green_analytical = np.append(Green_analytical, analytical_Green(value))

Green_analytical_real = Green_analytical.real
Green_analytical_imag = Green_analytical.imag

plt.plot(E, Green_real, label="Re($G_{11}$)", color="blue")
plt.plot(E, Green_imag, label="Im($G_{11}$)", color="blue")

plt.plot(E, Green_analytical_real, label="Analytical Re($G_{11}$)", linestyle="dashed", color="red")
plt.plot(E, Green_analytical_imag, label="Analytical Im($G_{11}$)", linestyle="dashed", color="red")

plt.grid(True)
plt.title("Approximation of $G_{11}$ by Iteration")
plt.ylabel("$G_{11}$")
plt.xlabel("Energy")
plt.legend()
plt.show()
