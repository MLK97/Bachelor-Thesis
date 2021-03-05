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

Now the following eq gets used:
G_11 = [1 - g_11 * t^+ * G_11 * t]^-1 * g_11,

where g_11, t^+, G_11 are matrices

g_11 is defined as:
g_11 = [(E)1-H]^-1 and
E = E +- i eta
"""

t = 1
epsilon = 0

t_mat = np.array([[t, 0], [0, t]])
t_dag = np.transpose(t_mat)

"""
Compute g_11 for a given Energy
"""
Green = np.complex128([]) # Stores values G_11(E)
Green_analytical = np.complex128([])
E = np.arange(-10, 10, 0.2) + 0.01j # Stores Energy values (x-axis)
g_11 = [] # Stores g_11(E) values

def calculate_g(E):
    mat = np.array([[(E-epsilon), -t],
                     [-t, (E-epsilon)]])
    mat_inv = inv(mat)
    return mat_inv

g_11 = np.array([calculate_g(value) for value in E])

def calculate_Green(g_11):
    G_11 = g_11
    G_11_prev = G_11
    G_11 = np.matmul(inv((np.identity(2) - multi_dot([g_11, t_mat, G_11, t_dag]))), g_11)
    while np.abs(G_11[0][0]-G_11_prev[0][0]) > 10**(-4):
        G_11_prev = G_11
        G_11 = np.matmul(inv((np.identity(2) - multi_dot([g_11, t_mat, G_11, t_dag]))), g_11)
    return G_11


Green = np.array([calculate_Green(value) for value in g_11])

Green_real = np.array([Green[i][0][0].real for i in range(len(Green))])
Green_imag = np.array([Green[i][0][0].imag for i in range(len(Green))])

"""
Computing the analytical equations (3.21) and (3.22)
from Zilly.
"""

#Problem for E = [-1, 1]
def analytical_Green(E, t):
    G_11 = 0
    if E < -2+t:
        G_11 = (E/2)+(1/4)*(np.sqrt((E-t)**2-4)+np.sqrt((E+t)**2-4))
    if -2 + t <= E and E <= 2 - t:
        # makes some dumb shit when rewriting the last sum without
        # explicitily telling that it is complex there
        G_11 = (E/2)+(1/4)*(-np.sqrt((E+t)**2-4)-1j*np.sqrt(4-(E-t)**2))
    if 2 - t <= E:
        G_11 = (E/2)+(1/4)*(-np.sqrt((E+t)**2-4)-np.sqrt((E-t)**2-4))
    return G_11


for value in E:
    Green_analytical = np.append(Green_analytical, analytical_Green(value, t))

Green_analytical_real = Green_analytical.real
Green_analytical_imag = Green_analytical.imag

plt.plot(E, Green_real, label="Re($G_{11}$)", color="red")
plt.plot(E, Green_imag, label="Im($G_{11}$)", color="red")

plt.plot(E, Green_analytical_real, label="Analytical Re($G_{11}$)", linestyle="dashed", color="blue")
plt.plot(E, Green_analytical_imag, label="Analytical Im($G_{11}$)", linestyle="dashed", color="blue")

plt.grid(True)
plt.title("$G_{11}$ by Iteration for $t = 2 \in [0, 2]$ and $\epsilon = 0$")
plt.ylabel("$G_{11}$")
plt.xlabel("Energy")
plt.legend()
plt.show()
