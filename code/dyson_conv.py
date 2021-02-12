import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

rc('font', **{'family': 'serif',
              'serif': ['Computer Modern'],
              'size': 20})
rc('text', usetex=True)
rc('figure', figsize=(11.69, 8.27))


"""
The Dyson equation is of the form
G = g + gVG,
where G is a matrice.

Eventually, for the components we arrive at
G_11 = g_11 + g_11 t^2 G_11^2,

where g_11 = 1/(E-epsilon) and
E = E +- i eta
"""

"""
Init Values:
t is transfer matrice element and known
g is the unperturbed Green function 
"""
t = 1
g_11 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

G = np.empty((7,40))

diff = np.empty(40)
diff[0] = 0
x = np.arange(40)

for a in range(len(g_11)):
    G[a][0] = g_11[a]

    print("t =", t, "\t", "energie =", g_11[a])
    print("Step\tG[n-1]\tG[n]\tDifference")

    for n in x:
        G[a][n] = g_11[a] + g_11[a] * t**2 * G[a][n-1]**2
        diff[n] = np.abs(G[a][n-1] - G[a][n])

        print(n, "\t", G[a][n-1], "\t", G[a][n], "\t", diff[n])
    
# plt.plot(x, diff, marker="o", label="Green Convergence")
plt.plot(x, G[0], marker=".", label="$G_{11}$ @ 1/(E - $\epsilon$) = 0.1")
plt.plot(x, G[1], marker="v", label="$G_{11}$ @ 1/(E - $\epsilon$) = 0.2")
plt.plot(x, G[2], marker="s", label="$G_{11}$ @ 1/(E - $\epsilon$) = 0.4")
plt.plot(x, G[4], marker="x", label="$G_{11}$ @ 1/(E - $\epsilon$) = 0.5")
plt.plot(x, G[5], marker="D", label="$G_{11}$ @ 1/(E - $\epsilon$) = 0.6")
plt.plot(x, G[6], marker="p", label="$G_{11}$ @ 1/(E - $\epsilon$) = 0.7")

plt.title("Numerical Evaluation of $G_{11}$ for fixed E")
plt.ylabel("$G_{11}$")
plt.xlabel("Step")

plt.ylim(0, 2)
plt.grid(True)
plt.legend()
plt.show()
