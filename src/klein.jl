"""
```
klein(A,B,C,t,k0,shock_exp,jumps)
```

* `A,B,C` see bellow
* `t` is the number of periods to be simulated
* `k0` is the initial condition for the stable part of the system.
* `shock_exp` is the expected value of the future shocks, usually zero.
* `jump` is a vector that tells the position of the jump variables

Solves a rational expectation models using Klein(2000) method. The model has to be written in the following way:

\$A E_{t}(\\mathbf{x_{t+1}}) = B \\mathbf{x_t} + C \\mathbf{\\varepsilon_t}\$

Where \$\\mathbf{x_t}\$ is the vector we are interested in and \$\\varepsilon_t\$ is a vector of random shocks. For now, the jumps variables must be the last variables in the system.

To guarantee that there is a stable solution, we need to check the Blanchard Khan conditions: the number of eigenvalues bigger than 1 is equal to the number of non predeterminated variables in the problem (e.g. variables that depend on the expectation of its future value), which are also called jump variables.

This method accepts a matrix A that is singular.
"""
function klein(A::AbstractArray,B::AbstractArray,C::AbstractArray,t::Int,k0::AbstractArray,shock_exp::AbstractArray,jumps::AbstractArray)

    F = schur(A,B)

    gen_eig = F.beta ./ F.alpha
    select = abs.(gen_eig) .< 1
    ordschur!(F,select)
    gen_eig = F.beta ./ F.alpha #Para acertar a ordenação

    check_bk = size(A,2) - length(jumps)
    greater1 = sum(abs.(gen_eig) .> 1)
    less1 = sum(abs.(gen_eig) .< 1)

    if(greater1!=length(jumps))
        error(string(length(jumps), " jump variables but ", greater1, " eigens > 1"))
    elseif(check_bk != less1)
        error(string(test2, " predeterminate, but ", less1, "eigens < 1"))
    else (@info "Blanchard Khan conditions satisfied!")
    end

    Q = F.Q
    Z = F.Z
    T = F.T
    S = F.S

    u_index = findfirst(abs.(gen_eig) .> 1)

    T11 = T[1:(u_index-1),1:(u_index-1)]
    T12 = T[1:(u_index-1),u_index:size(T,2)]
    T21 = T[u_index:size(T,1),1:u_index-1]
    T22 = T[u_index:size(T,1),u_index:size(T,2)]

    Q1 = Q[1:(u_index-1),:]
    Q2 = Q[u_index:size(Q,1),:]

    S11 = S[1:(u_index-1),1:(u_index-1)]
    S12 = S[1:u_index-1,u_index:size(T,2)]
    S21 = S[u_index:size(T,1),1:u_index-1]
    S22 = S[u_index:size(T,1),u_index:size(T,2)]

    Z11 = Z[1:(u_index-1),1:(u_index-1)]
    Z12 = Z[1:(u_index-1),u_index:size(T,2)]
    Z21 = Z[u_index:size(T,1),1:(u_index-1)]
    Z22 = Z[u_index:size(T,1),u_index:size(T,2)]

    u = -inv(T22)*inv(I-inv(T22)*S22)*Q2.*C'*shock_exp
    s0 = inv(Z11)*(k0-Z12*u)

    s = Array{Float64}(undef,t,check_bk)
    s[1,:] = s0

    for j = 1:(t-1)
        s[j+1,:] = inv(S11)*T11*s[j,:]+inv(S11)*T12*u-inv(S11)*S12*u+inv(S11)*Q1.*C'*shock_exp
    end

    u_aux = Array{Float64}(undef,t,length(jumps))
    for j = 1:t
        u_aux[j,:] = u
    end

    x_star = [s u_aux]
    x = similar(x_star)
    for j = 1:t
        x[j,:] = Z*x_star[j,:]
    end
    return x
end
