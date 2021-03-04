import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from numpy.linalg import inv
from numpy.linalg import multi_dot

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
rc('text', usetex=True)
rc('figure', figsize=(11.69, 8.27))

"""
Display if and how fast G_11 converges.

The following eq gets used:
G_11 = [1 - g_11 * t^+ * G_11 * t]^-1 * g_11,

where g_11, t^+, G_11 are matrices

g_11 is defined as:
g_11 = [(E)1-H]^-1 and
E = E +- i eta
"""

t = 1
epsilon = 0.5
t_mat = np.array([[0, 0], [t, 0]])
t_dag = np.transpose(t_mat)

"""
Compute g_11 for a given Energy
"""

E = np.arange(-10, 10, 0.2) + 0.01j # Stores Energy values (x-axis)
g_11 = [] # Stores g_11(E) values
conv = np.array([])

def calculate_g(E):
    mat = np.array([[(E-epsilon), -t],
                    [-t, (E-epsilon)]])
    mat_inv = inv(mat)
    return mat_inv

g_11 = np.array([calculate_g(value) for value in E])

def calculate_Green(g_11, conv):
    G_11 = g_11
    for n in range(300):
        G_11 = np.matmul(inv((np.identity(2) - multi_dot([g_11, t_mat, G_11, t_dag]))), g_11)
        conv = np.append(conv, G_11[0][0])
    return (G_11, conv)

Green = np.zeros(len(g_11)) # Stores values G_11(E)
Green = np.array([calculate_Green(value, conv) for value in g_11], dtype=object)

iteration = np.arange(0, 300, 1)
x_axis = np.array([energy*np.ones(300) for energy in E])

for i in range(len(Green)):
    deviation = np.array([np.abs(value.imag-Green[i][0][0][0].imag) for value in Green[i][1]])
    #print(deviation)
    plot = plt.scatter(iteration, x_axis[i], cmap='hot', c=deviation, alpha=0.8)
"""
pro energie gibt es 300 iterationen:
(energy, iteration)

for every energy point I need to make 300 energies
"""
# make labels
plt.title('Deviation Im($G_{11}$)')
plt.xlabel('Iteration')
plt.ylabel('Energy')
# plt.zlabel('Deviation')
plt.colorbar(plot)
plt.show()
