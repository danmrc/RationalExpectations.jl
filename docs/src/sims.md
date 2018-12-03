# Sims Solver

```@docs

sims(G0::AbstractArray,G1::AbstractArray,Pi::AbstractArray,Psi::AbstractArray)

irf(Theta1::AbstractArray,Theta2::AbstractArray,t::Int,shock)
```

## Issues

Although most implementations force to use the complex Schur, we don`t do it here (and discourage it). Here what happens when you use the schur decomposition with complex $\Gamma_0$ and $\Gamma_1$:

```@example complex_er

using Plots
using RationalExpectations

sigma = 1
phi_pi = 1.5
phi_y = 0.5/4
beta = 0.99
phi = 1
alph = 1/3
ep = 6
theta = 2/3
rho_v = 0.5

Theta = (1-alph)/(1-alph+alph*ep)
lambda = (1-theta)*((1-beta*theta)/theta)*Theta
kappa = lambda*(sigma+(phi+alph)/(1-alph))

G0 = [[1 0 0 0];[-1 1 0 0];[0 -1 sigma 1];[0 0 0 beta]]
G1 = [[rho_v 0 0 0];[0 0 phi_y phi_pi];[0 0 sigma 0];[0 0 -kappa 1]]
Psi = [1 0 0 0]'
Pi = [[0 phi_pi 0 1];[0 phi_y sigma -kappa]]'

sol_sims = sims(complex(G0),complex(G1),Pi,Psi)

resul = irf(real(sol_sims.Theta1),real(sol_sims.Theta2),12,0.25)

```

Let's plot the IRF of the shock to the monetary policy, that is an AR(1) with $\rho_v = 0.5$:

```@example complex_er

plot(resul[:,1], leg = false)

```
Solutions are welcomed. On the other hand, `klein` _does not_ suffer from this issue.
