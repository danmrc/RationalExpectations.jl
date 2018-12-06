We will implement the model from Galí(2008), chapter 3. We will need to write two separated models, as the sintax for the Sims solver is different of the sintaxe of the Klein solver. First, let's use the same calibration as Galí(2008):

```@example1

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

```

The equations are:

$\pi_t = \beta E_t(\pi_{t+1}) + \kappa\tilde{y}_t\\
\tilde{y}_t = -\frac{1}{\sigma}(i_t - E_t(\pi_{t+1})) + E_t(\tilde{y}_{t+1})\\
i_t = \rho + \phi_{\pi}\pi_t + \phi_{\tilde{y}_t}\tilde{y}_t + v_t\\
v_t = \rho_v v_{t-1} + \varepsilon_v\\$


We ignore the $r^n_t$ term as Galí(2008) does.

## Klein

The ordering of the variables for this model is:

$x_{t+1} =
\begin{pmatrix}
v_{t+1}\\
i_t\\
E_t(\pi_{t+1})\\
E_t(\tilde{y}_{t+1})\\
\end{pmatrix}$

Notice that we are using $\mathbf{v_{t+1}}$ on the left hand side. Using $v_{\mathbf{t}}$ will generate the wrong matrices. We will write this in a way that is consistent with the Klein method (equation 1), so the matrices are

$A =
\begin{pmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & \beta & 0\\
0 & -1 & 1 & \sigma\\
\end{pmatrix},

B = \begin{pmatrix}
\rho_v & 0 & 0 & 0\\
1 & 0 & \phi_\pi & \phi_y\\
0 & 0 & 1 & -\kappa\\
0 & 0 & 0 & \sigma\\
\end{pmatrix},

C = \begin{pmatrix}
1\\
0\\
0\\
0\\\end{pmatrix}$

As in Gali(2008), we will set a initial monetary poliocy shock of 0.25. This works as an initial condition for the model. As usual, we set that the expected value of the shocks are zero. We will receive the impulse response function automatically.

```@example1

using Plots, RationalExpectations

A = [[1 0 0 0];[0 1 0 0]; [0 0 beta 0]; [0 -1 1 sigma]]
B = [[rho_v 0 0 0];[1 0 phi_pi phi_y];[0 0 1 -kappa];[0 0 0 sigma]]

C = [1;0; 0; 0]

k0 = [0.25;0]

t=12

choque = [0.; 0;0; 0]

klein_sol = klein(A,B,C,t,k0,choque,[3 4])

```
We can plot the elements of `klein_sol` to see the irf

## Sims

Sims methods requires that we write expectations error, e.g. $\eta_t^{\pi} = \pi_t - E_{t-1}(\pi_t)$. We can work it to obtain the following matrices:

$\Gamma_0 = \begin{pmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & -1 & \sigma & 1\\
0 & 0 & 0 & \beta\\
\end{pmatrix},

\Gamma_1 = \begin{pmatrix}
 \rho_v & 0 & 0 & 0\\
1 & 0 & \phi_y & \phi_pi\\
0 & 0 & \sigma & 0\\
0 & 0 & -\kappa & 1\\
\end{pmatrix},

\Psi = \begin{pmatrix}
1\\
0\\
0\\
0\\
\end{pmatrix},

\Pi = \begin{pmatrix}
0 & 0\\
\phi_y & \phi_pi\\
\sigma & 0\\
-\kappa & 1\\
\end{pmatrix}$

And the variables are ordered as:

$x_{t+1} = \begin{pmatrix}
v_{t+1}\\
i_t\\
E_t(\tilde{y}_{t+1})\\
E_t(\pi_{t+1})\\
\end{pmatrix}$

See the end of this article for the whole maths of this transformation. Here is it, in Julia:

```@example1

G0 = [[1 0 0 0];[0 1 0 0];[0 -1 sigma 1];[0 0 0 beta]]
G1 = [[rho_v 0 0 0];[1 0 phi_y phi_pi];[0 0 sigma 0];[0 0 -kappa 1]]
Psi = [1 0 0 0]'
Pi = [[0 phi_pi 0 1];[0 phi_y sigma -kappa]]'

sol_sims = sims(G0,G1,Pi,Psi)

resul = irf(sol_sims,12,0.25)

```

## Analytical Solution and comparison

Last, but not least, Galí(2008) gives an analytical solution for $\tilde{y}_t$ and $\pi_t$ They are:

$\tilde{y}_t = -(1-\beta{}\rho_v) \Lambda_v v_t\\
\pi_t = -\kappa\Lambda_v v_t$

And $\Lambda_v = \frac{1}{(1-\beta{}\rho_v)[\sigma(1-\rho_v)+\phi_y]+\kappa(\phi_{\pi}-\rho_v)}$

Lets put the analytical solutions in Julia:

```@example1

Lambda_v = 1/((1-beta*rho_v)*(sigma*(1-rho_v)+phi_y)+kappa*(phi_pi - rho_v))
true_y(v) = -(1-beta*rho_v)*Lambda_v*v
true_pi(v) = -kappa*Lambda_v*v

true_path = zeros(13,3)
initial_shock = 0.25
shock = zeros(13)
shock[2] = initial_shock

for j = 2:13
    true_path[j,1] = rho_v*true_path[j-1] + shock[j]
    true_path[j,2] = true_y(true_path[j,1])
    true_path[j,3] = true_pi(true_path[j,1])
end

```

We are able to compare the analytical solution with the estimated solutions. First, the shock on $v_t$

```@example1

plot(true_path[2:13,1], lab = "Analytical Solution", linestyle = :dash, lw = 3)
plot!(resul[:,1], lab = "Gensys Answer", linestyle = :dash, lw = 3)
plot!(klein_sol[:,1], lab = "Klein Answer", linestyle = :dot, lw = 3)

```

Here is the shock in the output gap:

```@example1

plot(true_path[2:13,2], lab = "Analytical Solution", linestyle = :dash,lw = 3)
plot!(resul[:,3], lab = "Gensys Answer", linestyle = :dash, lw = 3)
plot!(klein_sol[:,4], lab = "Klein Answer", linestyle = :dot, lw = 3)

```

And the shock in the inflation:

```@example1

plot(4*true_path[2:13,3], lab = "Analytical Solution", linestyle = :dash, lw = 3)
plot!(4*resul[:,4], lab = "Gensys Answer", linestyle = :dash, lw = 3)
plot!(4*klein_sol[:,3], lab = "Klein Answer", linestyle = :dot, lw = 3)

```
