"""
```
sims(G0,G1,Pi,Psi)
```
This function solves a rational expectations model using Sims(2000) (aka gensys). The model should be written in the following way:

\$\\Gamma_0 \\mathbf{x_{t+1}} = \\Gamma_1 \\mathbf{x_{t}} + \\Psi \\mathbf{\\varepsilon_t} + \\Pi \\mathbf{\eta_t}\$

The arguments of the function are:

* G0 is \$\\Gamma_0\$
* G1 is \$\\Gamma_1\$
* Pi is \$\\Pi\$
* Psi is \$\\Psi\$

Where \$\\Gamma_0\$ can be singular. \$\\eta_t\$ is the expectation error and \$E_{t}(\\nu_{t+1})=0\$ by construction.

Sims`s solver does not require that we say what variables are jump variables.

## Value

The function returns theta1 and theta2, such that:

\$y_t = \\Theta_1 y_{t-1} + \\Theta_2 \\varepsilon_t\$

## Note

Although the original algorithm allows cases with multiple equilibrium (sunspots), this has not been implemented thus far.
"""
function sims(G0,G1,Pi,Psi)

    F = schur(G0,G1)
    sel = abs.(F.beta ./ F.alpha) .< 1
    ordschur!(F, sel)
    u_index = findfirst(abs.(F.beta ./F.alpha) .> 1)

    Z,Q = F.Z, F.Q

    n = size(F.S,2)
    S11 = F.S[1:(u_index-1),1:(u_index-1)]
    S12 = F.S[1:(u_index-1),(u_index):n]
    S22 = F.S[(u_index):n,(u_index):n]

    T11 = F.T[1:(u_index-1),1:(u_index-1)]
    T12 = F.T[1:(u_index-1),(u_index):n]
    T22 = F.T[(u_index):n,(u_index):n]

    Q1 = Q[1:(u_index-1),:]
    Q2 = Q[u_index:n,:]

    aux = Q2*Pi
    m = size(Pi,2)
    svd_dec = svd(aux)
    if m > size(svd_dec.S,1)
        error("No unique Stable Solution. Not implemented yet")
    elseif m == size(svd_dec.S,1)
        @info "Unique Solution Available"
    else
        error("No stable solution")
    end

    xi = Q1*Pi*pinv(aux)

    bb0= S12-xi*S22
    bb1 = [[S11 bb0];[zeros(size(S22,1),size(S11,2)) Matrix{Float64}(I,size(S22,1),size(S22,1))]]
    bb2= T12-xi*T22
    bb3 = [[T11 bb2];zeros(size(T22,1),(size(T11,2)+size(bb2,2)))]
    bb4 = [Q1 - xi*Q2;zeros(size(T22,1),(size(T11,2)+size(bb2,2)))]

    theta1 = Z*inv(bb1)*bb3*Z'
    theta2 = Z*inv(bb1)*bb4*Psi

    return (theta1 = theta1, theta2 = theta2)
end

function simulate_sys(Theta1,Theta2,t,burn)
    y = zeros((t+burn),size(Theta1,1))
    for j = 2:(t+burn)
        y[j,:] = Theta1*y[j+1,:] + Theta2*randn(size(Theta1,1))
    end
    y = y[burn+1:t,:]
    return y
end

"""
```
irf(Theta1,Theta2,t,shock)
```

Generates the IRF from the matrices calculate using the gensys.

* Theta1 is the theta1 matrix from `sims`
* Theta2 is the theta1 matrix from `sims`
* t is the number of periodos to be simulated
* shock is the size of the shock

See also [`sims(G0,G1,Pi,Psi)`] (@ref)
"""
function irf(Theta1,Theta2,t,shock)
    resp = zeros((t+1),size(Theta1,1))
    for j = 1:(t+1)
        resp[j,:] = Theta1^(j-1)*Theta2*shock
    end
    return resp
end
