using LinearAlgebra

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

    weird_e = Q1*Pi*pinv(aux)

    bb0= S12-weird_e*S22
    bb1 = [[S11 bb0];[zeros(size(S22,1),size(S11,2)) Matrix{Float64}(I,size(S22,1),size(S22,1))]]
    bb2= T12-weird_e*T22
    bb3 = [[T11 bb2];zeros(size(T22,1),(size(T11,2)+size(bb2,2)))]
    bb4 = [Q1 - weird_e*Q2;zeros(size(T22,1),(size(T11,2)+size(bb2,2)))]

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

function irf(Theta1,Theta2,t,shock)
    resp = zeros((t+1),size(Theta1,1))
    for j = 1:(t+1)
        resp[j,:] = Theta1^(j-1)*Theta2*shock
    end
    return resp
end
