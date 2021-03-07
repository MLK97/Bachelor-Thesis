import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from numpy.linalg import inv
from numpy.linalg import multi_dot

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
rc('text', usetex=True)
rc('figure', figsize=(11.69, 8.27))

"""
New approach for approximating G_11
via itertion and starting at g_11

G_11 = [1 - g_11 * t^+ * G_11 * t]^-1 * g_11,
Now the following eq gets used:

where g_11, t^+, G_11 are matrices

g_11 is defined as:
g_11 = [(E)1-H]^-1 and
E = E +- i eta
"""

t = 0.5
epsilon = 0

t_mat = np.array([[1, 0], [0, 1]])
t_dag = np.transpose(t_mat)

"""
Compute g_11 for a given Energy
"""
Green = np.complex128([]) # Stores values G_11(E)
Green_analytical = np.complex128([])
E = np.arange(-20, 20, 0.2) + 0.01j # Stores Energy values (x-axis)
g = [] # Stores g_11(E) values

def calculate_g(E):
    mat = np.array([[(E-epsilon), -t],
                     [-t, (E-epsilon)]])
    mat_inv = inv(mat)
    return mat_inv

g = np.array([calculate_g(value) for value in E])

def calculate_Green(g):
    G = g
    G_prev = G
    G = np.matmul(inv((np.identity(2) - multi_dot([g, t_mat, G, t_dag]))), g)
    while np.abs(G[1][0]-G_prev[1][0]) > 10**(-4):
        G_prev = G
        G = np.matmul(inv((np.identity(2) - multi_dot([g, t_mat, G, t_dag]))), g)
    return G


Green = np.array([calculate_Green(value) for value in g])

Green_real = np.array([Green[i][1][0].real for i in range(len(Green))])
Green_imag = np.array([Green[i][1][0].imag for i in range(len(Green))])

"""
Computing the analytical equations (3.21) and (3.22)
from Zilly.
"""

#Problem for E = [-1, 1]
def analytical_Green(E, t):
    G = 0
    if E < -2-t:
        G = (E/2)+(1/4)*(np.sqrt((E-t)**2-4)+np.sqrt((E+t)**2-4))
    if -2 - t <= E and E <= -t + 2:
        G = (E/2)+(1/4)*(np.sqrt((E-t)**2-4)-1j*np.sqrt(4-(E+t)**2))
    if -t + 2 < E and E < t - 2:
        G = (E/2)+(1/4)*(np.sqrt((E-t)**2-4)-np.sqrt((E+t)**2-4))
    if t - 2 <= E <= t + 2:
        G = (E/2)+(1/4)*(-np.sqrt((E+t)**2-4)-1j*np.sqrt(4-(E-t)**2))
    if 2 + t < E:
        G = (E/2)+(1/4)*(-np.sqrt((E+t)**2-4)-np.sqrt((E-t)**2-4))
    return G


for value in E:
    Green_analytical = np.append(Green_analytical, analytical_Green(value, t))

Green_analytical_real = Green_analytical.real
Green_analytical_imag = Green_analytical.imag

plt.plot(E, Green_real, label="Re($G_{11}$)", color="red")
plt.plot(E, Green_imag, label="Im($G_{11}$)", color="red")

#plt.plot(E, Green_analytical_real, label="Analytical Re($G_{12}$)", linestyle="dashed", color="blue")
#plt.plot(E, Green_analytical_imag, label="Analytical Im($G_{12}$)", linestyle="dashed", color="blue")

plt.grid(True)
plt.title("$G_{12}$ by Iteration for $t = 3 \in [2, \infty)$ and $\epsilon = 0$")
plt.ylabel("$G_{11}$")
plt.xlabel("Energy")
plt.legend()
plt.show()
