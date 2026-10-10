# schrodinger-solver

Finds the energy levels and wavefunctions of a particle in one dimension, numerically:

second derivative by finite differences -> written as a matrix -> put into the Hamiltonian -> for a given space and potential, the Hamiltonian is now a matrix -> its eigenvalues are the energies and its eigenvectors are the wavefunctions

## Why

I wanted to understand why energy levels are discrete.

## The physics

Energy is kinetic plus potential:

    E = p^2/2m + V(x)

In quantum mechanics momentum becomes an operator, p = -i hbar d/dx, so p^2 = -hbar^2 d^2/dx^2. The energy then becomes the Hamiltonian operator:

    H = -(hbar^2/2m) d^2/dx^2 + V(x)

The allowed energies are the solutions of H phi = E phi.

The container is built into H through the potential. Where V < E, phi oscillates. Where V > E the region is classically forbidden and phi has to decay, because the solution has to satisfy H phi = E phi and also be normalisable (total probability 1). Only certain energies manage both, so the levels are discrete.

On a grid with spacing h, the second derivative becomes a matrix:

    phi''(x_i) ~ (phi_(i-1) - 2*phi_i + phi_(i+1)) / h^2

Adding the potential on the diagonal makes H a matrix too, so H phi = E phi is a matrix eigenvalue problem. The eigenvalues are the energies and the eigenvectors are the wavefunctions.

Three potentials: hard walls (infinite square well), finite walls (finite square well) and a bowl (harmonic oscillator, V = x^2/2). Units are hbar = m = 1.

## Results

All three match the analytic results:

| Potential | Analytic | Numerical |
|---|---|---|
| Hard walls, n = 1, 2, 3 | 4.935, 19.739, 44.413 | 4.934, 19.733, 44.381 |
| Bowl, n = 0, 1, 2 | 0.5, 1.5, 2.5 | 0.4999, 1.4996, 2.4990 |
| Finite walls, all 5 bound states | -19.08, -16.35, -11.91, -6.00, -0.027 | -19.08, -16.35, -11.92, -6.03, -0.041 |

## Known limitations

- The second derivative is approximated on a grid, so the answers are only as good as the step size h. The error gets bigger for higher states: the 5th hard wall level is off by 0.2%.
- The clearest case is the top finite well state. It is only just bound, so it is the most sensitive: -0.041 against -0.027. Making the box bigger does not fix it, only a smaller step size would.
- The code is one script that runs each case in turn, not reusable functions.

## Running it

Needs Python with numpy and matplotlib:

    pip install numpy matplotlib
    python schrodinger-solver.py

It prints each energy next to the expected value, then opens the plots one at a time. Close one to see the next. In PyCharm, open the file and press Run.
