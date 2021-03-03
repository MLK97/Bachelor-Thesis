import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
rc('text', usetex=True)
rc('figure', figsize=(11.69, 8.27))
"""
The main problem of this method
is its very high instability
when approximating G_11 via
iteration and starting at g_11.

Since G_11 appears as square,
every value that is
- higher than |0.5|
- extremely small 10e-100
causes very big fluctuations of the iteration
"""


"""
The Dyson equation is of the form
G = g + gVG,
where G is a matrice.

Eventually, for the 11-component we arrive at
G_11 = g_11 + g_11 t^2 G_11^2,

where g_11 = 1/(E-epsilon) and
E = E +- i eta
"""

"""
Init Values:
t is transfer matrice element and assumed known
g is the unperturbed Green function 
"""

"""
This code is divided into 3 functions

g_values(): creates the matrix containing g_11 values used in approximate_G()
approximate_G(): approximates G_11 with the values t, g_11
plot_G(): takes the best approximation of approximate_G() and plots them as G(E)

plot_G()'s plot should be identical/similar to Fig 5.2 in Cuevas
"""

def approximate_G(t, g_11):
    G = np.empty((13,40)) # change first value according to len(E)
    x = np.arange(40)
    g_11 = g_11.real
    for m in range(len(g_11)):
        if abs(g_11[m]) > 0.5:
            """ checks if value is 0.5 and thus would diverge anyway"""
            print("Initial Value =", g_11[m], " got skipped.")
            continue
        G[m][0] = g_11[m]

        print("Current Value of m =", m)
        print("Initial Value =", g_11[m])
        print("Step\tG[n]")

        for n in x:
            G[m][n] = g_11[m] + g_11[m] * t**2 * G[m][n-1]**2
            print(n, "\t", G[m][n])

    print(G[0])
    #for l in range(0, len(g_11)-1):
    plt.plot(x, G[0], marker="v", label=g_11[1])

    #plt.ylim(-50000, 50000)
    plt.grid(True)
    #plt.legend()
    plt.title("Approximation of $G_{11}$ by Iteration")
    plt.ylabel("$G_{11}$")
    plt.xlabel("Step")
    plt.show()
    return G

def plot_G(G, E):
    """
    From approximate_G() we got
    approximated values for G in (2,2) form
    G[n][max-length] is the final approximated value
    for the n-th variation of g_11.

    approximate_G() currently varies g_11
    by changing the Energy.
    Thus we can plot G(E).
    """
    best_approx = len(G[0])-1
    Green = np.empty(len(G))

    for n in range(len(Green)-1):
        print(G[n][best_approx])
        Green[n] = G[n][best_approx]
    
    plt.plot(E, Green, marker="p", label="$G_{11}$ @ 1/(E - $\epsilon$) = 0.7")

    plt.title("Numerical Evaluation of $G$ over E")
    plt.ylabel("G")
    plt.xlabel("E")

    # plt.ylim(0, 2)
    plt.grid(True)
    #plt.legend()
    plt.show()

def g_values(E):
    """
    g_11 is the unperturbed function of G_11
    For the half infinite chain it is defined as
    g_11 = 1/(E-epsilon) and
    E = E +- i eta
    """

    """
    Create Array of Energy values with
    different etas
    """
    E_real = E
    #E_imag = 1j*np.array([-0.1, -0.01, -0.001, 0.001, 0.01, 0.1])
    E_imag = 1j*np.array([0.001])
    Energy = np.array([0])
    for n in E_real:
        Energy = np.append(Energy, n + 1j * E_imag)

    """
    Compute g_11 with the formula
    g_11 = 1/(E-epsilon).
    At the moment epsilon is fixed and
    epsilon = 0.5
    """
    epsilon = 0.5
    g_11 = 1/(Energy-epsilon)
    print(len(g_11))
    return g_11
    
if __name__ == "__main__":
    t = 1 # Value used in Zilly
    epsilon = [0.1, 0.2, 0.3, 0.4, 0.5] # Zilly uses 0.5  
    E = np.array([-3, -2.5, -2, 1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3])
    
    g_11 = g_values(E)
    G = approximate_G(t, g_11)
    E = np.array([-3, -2, -1, 1, 2, 3])
    #plot_G(G, E)
