import numpy as np
import matplotlib.pyplot as plt

def second_derivative_matrix(N, h):
    diagonal = np.eye(N) * -2
    upper_off_diagonal = np.eye(N, k=1)
    lower_off_diagonal = np.eye(N, k=-1)
    D2 = (1/h**2) * (diagonal + upper_off_diagonal + lower_off_diagonal)
    return D2

def hamiltonian(N, h, V): #note to self, h = kg m^2 / t, hbar, kg, L set to 1
    D2 = second_derivative_matrix(N, h)
    H = -0.5 * D2
    V_matrix = np.diag(V)
    return H + V_matrix

def solve(N, x_min, x_max, V):
    L = x_max - x_min
    h = L / (N+1)
    x = np.linspace(x_min + h, x_max -h , N)
    H = hamiltonian(N, h, V)
    E, psi = np.linalg.eigh(H)
    return E, psi, h, x

def normalise(psi, h, n_states):
    return [psi[:, i] / np.sqrt(h) for i in range (n_states)]





if __name__ == "__main__":

    N = 3
    h = 0.1
    u2 = 0.9
    u3 = 1.0
    u4 = 1.08
    u = np.array([u2, u3, u4])
    D2U = second_derivative_matrix(N, h) @ u
    print(D2U)

    x = np.linspace(-1, 1, N)
    V_func = 0.5 * x ** 2
    H = hamiltonian(N, h, V_func)
    print(H)

    N = 100
    x_min = 0
    x_max = 1
    V = np.zeros(N)
    E, psi, h, x = solve(N, x_min, x_max, V)
    for i in range(5):
        print(E[i], np.pi**2 * (i+1)**2 / 2)

    L = x_max - x_min
    h = L / (N+1)
    x = np.linspace(h, L-h, N)
    states = normalise(psi, h, 3)

    for n,phi in enumerate(states):
        print(n+1, np.sum(phi**2)*h)
        if phi[0]<0:
            phi = -phi
        plt.plot(x, phi, label=f"n={n+1}")
        plt.axhline(y=0, color='k', linestyle='-')
    plt.title("square well eigenfunctions")
    plt.legend()
    plt.show()


    N = 200
    x_min = -5
    x_max = 5
    h = (x_max - x_min) / (N + 1)
    x = np.linspace(x_min + h, x_max - h, N)
    V = 0.5 * (x ** 2)
    E, psi, h, x = solve(N, x_min, x_max, V)
    states = normalise(psi, h, 3)


    for n in range(5):
        print(n, "expected E", n + 0.5, "E obtained", E[n])

    for n, phi in enumerate(states):
        if phi[0] < 0:
            phi = -phi
        plt.plot(x, phi, label=f"n={n}")
        plt.axhline(y=0, color='k', linestyle='-')
    plt.title("harmonic oscillator eigenfunctions")
    plt.legend()
    plt.show()

    x_min = -10
    x_max = 10
    N = 400
    h = (x_max - x_min) / (N + 1)
    x = np.linspace(x_min + h, x_max - h, N)
    V = []
    for xn in x:
        if -1 < xn < 1:
            V.append(-20)
        else:
            V. append(0)
    V = np.array(V)
    E, psi, h, x = solve(N, x_min, x_max, V)
    phi_bound = np.sum(E < 0)
    states = normalise(psi, h, phi_bound)
    bound = E[E < 0]
    print(len(bound), "bound states:", bound)
    for n, phi in enumerate(states):
        plt.plot(x, phi, label=f"n={n}")
        plt.axhline(y=0, color='k', linestyle='-')
    plt.title("finite square well eigenfunctions")
    plt.legend()
    plt.show()

    top = phi_bound - 1
    phi = psi[:, top] / np.sqrt(h)
    plt.plot(x, phi)
    plt.axvline(-1, color = 'k', linestyle = '--')
    plt.axvline(1, color = 'k', linestyle = '--')
    plt.axhline(0, color = 'r', linewidth = 0.8)
    plt.title("L, N")
    plt.show()

    x_min = -20
    x_max = 20
    N = 800
    h = (x_max - x_min) / (N + 1)
    x = np.linspace(x_min + h, x_max - h, N)
    V = []
    for xn in x:
        if -1 < xn < 1:
            V.append(-20)
        else:
            V. append(0)
    V = np.array(V)
    E, psi, h, x = solve(N, x_min, x_max, V)
    phi_bound = np.sum(E < 0)
    bound = E[E < 0]

    top = phi_bound - 1
    phi = psi[:, top] / np.sqrt(h)
    plt.plot(x, phi)
    plt.axvline(-1, color = 'k', linestyle = '--')
    plt.axvline(1, color = 'k', linestyle = '--')
    plt.axhline(0, color = 'r', linewidth = 0.8)
    plt.title("L x 2, N x 2 ")
    plt.show()
    print(len(bound), "bound states:", bound)

    x_min = -40
    x_max = 40
    N = 1600
    h = (x_max - x_min) / (N + 1)
    x = np.linspace(x_min + h, x_max - h, N)
    V = []
    for xn in x:
        if -1 < xn < 1:
            V.append(-20)
        else:
            V.append(0)
    V = np.array(V)
    E, psi, h, x = solve(N, x_min, x_max, V)
    phi_bound = np.sum(E < 0)
    bound = E[E < 0]

    top = phi_bound - 1
    phi = psi[:, top] / np.sqrt(h)
    plt.plot(x, phi)
    plt.axvline(-1, color='k', linestyle='--')
    plt.axvline(1, color='k', linestyle='--')
    plt.axhline(0, color='r', linewidth=0.8)
    plt.title("L x 4, N x 4 ")
    plt.show()
    print(len(bound), "bound states:", bound)


    x_min = -80
    x_max = 80
    N = 3200
    h = (x_max - x_min) / (N + 1)
    x = np.linspace(x_min + h, x_max - h, N)
    V = []
    for xn in x:
        if -1 < xn < 1:
            V.append(-20)
        else:
            V.append(0)
    V = np.array(V)
    E, psi, h, x = solve(N, x_min, x_max, V)
    phi_bound = np.sum(E < 0)
    bound = E[E < 0]

    top = phi_bound - 1
    phi = psi[:, top] / np.sqrt(h)
    plt.plot(x, phi)
    plt.axvline(-1, color='k', linestyle='--')
    plt.axvline(1, color='k', linestyle='--')
    plt.axhline(0, color='r', linewidth=0.8)
    plt.title("L x 8, N x 8 ")
    plt.show()
    print(len(bound), "bound states:", bound)