using LinearAlgebra
using DSP
using StatsBase
using NumericalIntegration

import Base.+

+(A::Array{Float64,1},B) = A .+ B
matmul(A,B) = A * B

mutable struct cache
    k
    w
    Hk
    vmat
    Gwk
end

function fermi(w,μ,T)
    kT = 8.617332478e-5 * T
    if T == 0
        return float(w .< μ)
    else
        x = (w - μ) / kT
        return inv(1.0 + exp(x))
    end
end

function calc_DOS(model::model)
    if model.dim_k == 0
        eigvals = solve_eig(model)
        dw = model.cache.w[2] - model.cache.w[1]
        wmesh = append!(model.cache.w,model.cache.w[end]+dw)

        hist = fit(Histogram,vec(eigvals),wmesh)
        dos  = hist.weights / ( size(eigvals,1) * size(eigvals,2) )

        return wmesh,dos
    else
        if model.cache == 0
            @error "cache has not been setup. Please, warm up first."
        end
        Gw = tr.(sum.(model.cache.Gwk)) / length(model.cache.k)
        dos = -imag(Gw) / π

        return model.cache.w, dos
    end
end

function warmup_calc!(model::model,kmesh,wmesh,w_window; η=0.05)
    #check kmesh dim
    if length(kmesh) != model.dim_k
        @error "kmesh size must be same size as dim_k"
    end

    #warmup params
    k_vec = k_uniform_mesh(model, kmesh)
    w = LinRange(w_window[1],w_window[2],wmesh)
    nk = length(k_vec)
    nst = model.nstates

    #collect Ham
    @info "Collecting Hamiltonian matrix"
    H(k) = hamiltonian(model,k)
    Hk = H.(k_vec)

    #collect vk
    @info "Collectting Velocity matrix"
    vx = []
    vy = []
    vz = []

    dk = 1e-5
    if model.dim_k == 1
        for i in 1:nk
            dkx = [k_vec[i][1] + dk]
            v = (H(dkx) - Hk[i]) / dk
            push!(vx,real(v))
        end
    elseif model.dim_k == 2
        for i in 1:nk
            dkx = [k_vec[i][1] + dk, k_vec[i][2]]
            v = (H(dkx) - Hk[i]) ./ dk
            push!(vx,real(v))
            dky = [k_vec[i][1], k_vec[i][2] + dk]
            v = (H(dky) - Hk[i]) ./ dk
            push!(vy,real(v))
        end
    elseif model.dim_k == 3
        for i in 1:nk
            dkx = [k_vec[i][1] + dk, k_vec[i][2], k_vec[i][3]]
            v = (H(dkx) - Hk[i]) ./ dk
            push!(vx,real(v))
            dky = [k_vec[i][1], k_vec[i][2] + dk, k_vec[i][3]]
            v = (H(dky) - Hk[i]) ./ dk
            push!(vy,real(v))
            dkz = [k_vec[i][1], k_vec[i][2], k_vec[i][3] + dk]
            v = (H(dkz) - Hk[i]) ./ dk
            push!(vz,real(v))
        end
    end


    #collect Gwk
    @info "Collecting Gwk matrix"
    Gwk = []
    Id = Matrix(I,nst,nst)
    for (i,wx) in enumerate(w)
        Iw = Id * (wx + 1im*η)
        Gk = []
        for (ik,Hkx) in enumerate(Hk)
            push!(Gk,inv(Iw .- Hkx))
        end
        push!(Gwk,float.(Gk))
    end

    model.cache = cache(k_vec,w,Hk,(float.(vx),float.(vy),float.(vz)),Gwk)
end

function calc_OptCond(model::model; T,dir=[1,2])
    if model.cache == 0
        @error "cache has not been setup. Please, warm up first."
    end

    function KK(ww,s)
        isigma = zeros(Float64,length(s))
        for (iw,w1) in enumerate(ww)
            wsq = ww.^2 .- w1^2
            intg = s./ wsq
            intg[iw] = 0
            isigma[iw] = integrate(ww,intg,TrapezoidalFast())
        end
        return -2ww .* isigma ./ π
    end

    μ = model.mu
    nw = length(model.cache.w)
    nk = length(model.cache.k)

    @info "Preparing conductivity calculation"
    dw = model.cache.w[2] - model.cache.w[1]
    if nw % 2 != 0
        half = floor(Int,nw/2)
    else
        half = Int(nw/2) + 1
    end
    kT = 8.617332478e-5 * T

    sigma = []

    @info "Calculating optical conductivity..."

    for (iw,w) in enumerate(model.cache.w[half:end])
        tmp =  zeros(Float64,nw)
        sumTr = 0.0
        for (inu,nu) in enumerate(model.cache.w)
            inup = floor(Int64,(nu+w)/dw) + half
            if (inup > 0) && (inup <= nw)
                Awk = -imag(model.cache.Gwk[inup]) ./ π
                v2A = matmul.(model.cache.vmat[dir[2]],Awk)
                Awk = -imag(model.cache.Gwk[inu]) ./ π
                v1A = matmul.(model.cache.vmat[dir[1]],Awk)
                vAvA = matmul.(v1A,v2A)
                sumTr = sum(tr.(vAvA)) / nk
            end

            if w == 0 && T != 0
                fm = fermi(w,μ,T)
                fm = (fm * (1 - fm)) / kT
            elseif w == 0 && T == 0
                fm = 1.0
            elseif w != 0
                fm = (fermi(nu,μ,T) - fermi(nu+w,μ,T)) / w
            end

            tmp[inu] = fm * sumTr
        end
        sig = integrate(model.cache.w,tmp,TrapezoidalFast())
        push!(sigma,sig)
    end

    rw = length(sigma)
    hw = model.cache.w[(nw-rw+1):end]

    isigma = KK(hw,sigma)

    return hw,sigma .+ 1im*isigma
end

function calc_Epsilon(hw,sigma;eps0,eps_infty)
    return eps_infty .+ 1im .* sigma ./ (eps0 .* hw)
end

function calc_Reflectance(epsilon)
    sqeps = sqrt.(epsilon)
    sqeps = (sqeps .- 1) ./ (sqeps .+ 1)
    refle = abs.(sqeps).^2
    return refle
end
