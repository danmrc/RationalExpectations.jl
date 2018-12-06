
function doH(lags,leads,h::Vararg{Array})

    tam = length(h)
    H = h[1]
    for i = 2:tam
        H = hcat(H,h[i])
    end

    ret = AIMH(H,lags,leads,size(h[1],2))

    return(ret)

end

function f1aux(H,lags,leads,nvars,kmax)
    hth = H[:,((nvars*(lags+leads))+1):size(H,2)]
    k = 0
    while det(hth) == 0 && k < kmax
        qr_dec = qr(hth)
        q = qr_dec.Q'
        un = q[1,:]'
        uz = q[2,:]'
        block1 = uz*H[:,1:nvars*(lags+leads)]
        block2 = fill!(similar(block1),0)
        block3 = un*H
        H = [[block2 block1];block3]
        Z = [Z;block1]
        hth = H[:,((nvars*(lags+leads))+1):size(H,2)]
        k = k+1
    end
    Z = Z[2:size(Z,1),:]
    return (H = H, Z = Z, k=k)
end

function f1(arg::AIMH,kmax)
    anh = f1aux(arg.H,arg.lags,arg.leads,arg.nvars,kmax)
    thh = arg.nvars*(arg.lags+arg.leads) #indexes the last column before H_\theta, the H for the largest lag
    hth = anh.H[:,(thh+1):size(arg.H,2)]
    hh = anh.H[:,1:(thh)]
    G = inv(hth)*hh
    b1 = zeros(arg.nvars*arg.lags,arg.nvars*arg.lags)
    b2 = Matrix(I,arg.nvars*arg.leads,arg.nvars*arg.leads)
    if size(b2) == (0,0)
        A = [b1;G]
    else
        A = [b1 b2;G]
    end
    AIM1(H = anh.H,A = A, Z = anh.Z)
end

function orthogonal_space(A,Z)
     eig=eigen(A)
     check = abs.(eig.values) .> 1
     broots = findall(check)
     left_eig = inv(eig.vectors)
     V = left_eig[broots,:]
     if length(Z) == 0
         Q = V
     else if length(V) == 0
         Q = Z
     else
         Q = [Z;V]
     end
     return Q
 end

function f3(Q,lags,leads,nvars)
    cnt = size(Q,1)
    if cnt != leads*nvars
        error("No unique solution")
    end
    @info Unique Solution Available

    Q_l = Q[:,1:nvars*lags]
    Q_r = Q[:,(nvars*lags+1):size(Q,2)]

    B = -inv(Q_r)*Q_l

    return B
end

function AIM(model::AIMH,kmax = 500)
        first = f1(model,kmax)
        secd = orthogonal_space(first.A,first.Z)
        soll = f3(secd,model.lags,model.leads,model.nvars)
        br = soll[:,(model.nvars*model.lags+1):size(soll,2)]
        




        return
end
