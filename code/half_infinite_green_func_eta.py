import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import cmath

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
rc('text', usetex=True)
rc('figure', figsize=(11.69, 8.27))

"""
The numerical stable method for solving quadratic equations
"""

"""
The equation to solve is:
G_11 = g_11 + g_11*t^2*G_11^2
with g_11 = 1/(E-epsilon)

Thus (E-epsilon)G_11 = g_11 + t^2*G_11^2.

Thus:
a = t^2
b = (E-epsilon)
c = 1
"""
E = np.arange(-4, 4, 0.1)
solution = np.ndarray((len(E)-1), dtype=np.complex128)
solution_s = np.ndarray((len(E)-1), dtype=np.complex128)
solution_xs = np.ndarray((len(E)-1), dtype=np.complex128)
solution_xxs = np.ndarray((len(E)-1), dtype=np.complex128)
a = 1 # t, coupling energy
b = E - 5 # Epsilon like Zilly, Energy acts as x
c = 1

def quadSolver(a, b,c, tol = 1e-18):
    print('Equation: {0}x**2 + {1}x + {2}'.format(a,b,c))
    if a==b==0:
        if c!=0:
            print('Not a valid equation')
        else:
            print(' 0=0 is not an interesting equation')
        return

    if a==0:
        print('Single solution is x =', -c/b)
        return

    discriminant = b**2 - 4 * a * c
    if discriminant > 0:
        root1 = (-b + np.sqrt(discriminant))/ (2 * a)
        root2 = (-b - np.sqrt(discriminant))/ (2 * a)
        print('Has two roots:')
        print(np.complex128(root2))
        return np.complex128(root2)
    elif discriminant == 0:
        root1 = float(-b + np.sqrt(discriminant))/ (2 * a)
        print('Has a double root:')
        return np.complex128(root1)
    elif discriminant < 0:
        root1 = (-b + cmath.sqrt(discriminant))/ (2 * a)
        root2 = (-b - cmath.sqrt(discriminant))/ (2 * a)
        print('Has two complex roots:')
        print(root1)
        print(root2)
        return np.complex128(root2)

for n in range(len(E)-1):
    b = E[n] - 0.5
    solution[n] = quadSolver(a,b,c) # Eta = 0
    solution_s[n] = quadSolver(a,b-0.1j,c) # Eta = -0.1
    solution_xs[n] = quadSolver(a,b-0.01j,c) # Eta = -0.01
    solution_xxs[n] = quadSolver(a,b-0.001j,c) # Eta = -0.001

solution_real = solution.real
solution_imag = solution.imag

plt.plot(E[:79], solution_real, marker="v", label="Re($G_{11}), \eta=0$", markevery=2)
plt.plot(E[:79], solution_imag, marker="x", label="Im($G_{11}), \eta=0$", markevery=2)

plt.plot(E[:79], solution_s.real, marker="v", label="Re($G_{11}), \eta=-0.1$", markevery=2)
plt.plot(E[:79], solution_s.imag, marker="x", label="Im($G_{11}), \eta=-0.1$", markevery=2)

plt.plot(E[:79], solution_xs.real, marker="v", label="Re($G_{11}), \eta=-0.01$", markevery=2)
plt.plot(E[:79], solution_xs.imag, marker="x", label="Im($G_{11}), \eta=-0.01$", markevery=2)

plt.plot(E[:79], solution_xxs.real, marker="v", label="Re($G_{11}), \eta=-0.001$", markevery=2)
plt.plot(E[:79], solution_xxs.imag, marker="x", label="Im($G_{11}), \eta=-0.001$", markevery=2)


plt.grid(True)
plt.legend()

plt.title("Surface Green Function for different $\eta$")
plt.ylabel("$G_{11}$")
plt.xlabel("Energy")
plt.show()

