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
function sims(G0::AbstractArray,G1::AbstractArray,Pi::AbstractArray,Psi::AbstractArray)

    F = schur(G0,G1)
    sel = abs.(F.beta ./ F.alpha) .< 1
    ordschur!(F, sel)
    gen_eigen = abs.(F.beta ./F.alpha)
    u_index = findfirst(gen_eigen .> 1)

    unst = sum(gen_eigen .> 1)
    stb = sum(gen_eigen .< 1)

    Z,Q,S,T = F.Z, F.Q',F.S,F.T

    n = size(F.S,2)

    Q1 = Q[1:(u_index-1),:]
    Q2 = Q[u_index:n,:]

    aux = Q2*Pi
    m = size(Pi,2)
    svd_dec = svd(aux)
    r = size(svd_dec.S,1)
    if m > r
        @info "No unique Stable Solution"
    elseif m == r
        @info "Unique Solution Available"
    else
        error("No stable solution")
    end

    xi = Q1*Pi*pinv(aux)

    ixi = [Matrix{Float64}(I,size(xi,1),size(xi,2)) -xi]

    bb1 = ixi*S
    bb2 = ixi*T
    bb3 = ixi*Q

    b1 = [bb1; zeros(unst,unst) Matrix{Float64}(I,size(bb1,1),size(bb1,1))]
    b2 = [bb2;zeros(unst,size(bb2,2))]
    b3 = [bb3; zeros(unst,size(bb3,2))]

    theta1 = Z*inv(b1)*b2*Z'
    theta2 = Z*inv(b1)*b3*Psi

    if m > r
        V1 = svd_dec.V[1:r,:]
        V2 = svd_dec.V[r+1:size(svd_dec.V,1),:]
        Vs = V1*V1'
        bb5 = Q1*Pi*(Matrix{Float64}(I,size(Vs,1),size(Vs,2))-Vs)
        bb6 = [bb5;zeros(unst,size(bb5,2))]
        theta3 = Z*inv(b1)*bb6*V2
        ret = (theta1,theta2,theta3)
    else
    theta3 = nothing
    ret = SimsSol(theta1,theta2,theta3)
    end

    return ret
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
irf(model::SimsSol,t::Int,shock)
irf(Theta1,Theta2,t,shock)
```

Generates the IRF from the matrices calculate using the gensys.

This function allows you to pass the matrices directly or to provide the whole model.

* model receives the full output of the function `sims`
* Theta1 is the theta1 matrix from `sims`
* Theta2 is the theta1 matrix from `sims`
* t is the number of periodos to be simulated
* shock is the size of the shock

See also [`sims(G0,G1,Pi,Psi)`] (@ref)

"""
function irf(Theta1::AbstractArray,Theta2::AbstractArray,t::Int,shock)
    resp = zeros((t+1),size(Theta1,1))
    for j = 1:(t+1)
        resp[j,:] = Theta1^(j-1)*Theta2*shock
    end
    return resp
end

function irf(model::SimsSol,t::Int,shock)
    irf(model.Theta1::AbstractArray,model.Theta2::AbstractArray,t::Int,shock)
end
