using LinearAlgebra
using DSP
using StatsBase

####################################
# Base of the calc.jl
###################################
const HALF = 1/2

function integrate(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    return (x[2] - x[1]) * (HALF * (y[1] + y[end]) + sum(y[2:end-1]))
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

####################################
# Calc.jl function is started here
###################################

"""
    calc_DOS(model, nw, kmesh, wrange, eta)

Calculate density of state of the model.

# Optional Arguments
- `eta` : Lorentzian Width for approximating δ(E)

# Returns
- `w` : energy mesh
- `DOS` : density of state along the energy mesh

# Examples
```julia
# setup the model
tb = tb_model(2,2,lat,orb)
set_onsite! ...
set_hop! ...
# calculate dos
w, dos = cacl_DOS(tb, nw=100, kmesh=[100,100], wrange=(-4,4))
```
"""
function calc_DOS(model;nw,kmesh,wrange,eta=0.01)
    if length(kmesh) != model.dim_k
        error("number of kmesh must same as number of k direction")
    end
    if typeof(wrange) != Tuple(Int64,Int64)
        error("please specify right frequency range in Tuple")
    end

    k = k_uniform_mesh(model,kmesh)
    H(k) = hamiltonian(model,k)

    w = LinRange(wrange[1],wrange[2],nw)

    dos = zeros(Float64,size(w)...)

    @info "calculating dos"
    @inbounds for ik in 1:length(k)
        kp = k[ik]
        for (ie,ep) in enumerate(eigvals(H(kp)))
            dos += eta ./ ( (w .- ep).^2 .+ eta^2 )
        end
    end

    return w,dos / (kmesh[1]*kmesh[2])
end

"""
    calc_OptCond(model, hwmesh,hnw, T=10,dir=[1,2])

Calculate optical conductivity of the model in given direction by `dir`.
The optical conductivity is defined by Kubo formula

```math
\\text{Re}\\sigma_{\\alpha\\beta}(\\omega) = \\int_{-\\infty}^{\\infty} d\\nu \\frac{f(\\omega,\\mu,T) - f(\\nu+\\omega,\\mu,T)}{\\omega} \\sum_{\\mathbf{k}} \\text{Tr} \\; \\boldsymbol{\\mathsf{v_\\alpha}}(\\mathbf{k}) \\boldsymbol{\\mathsf{A(\\omega,\\mathbf{k})}}\\boldsymbol{\\mathsf{v_\\beta}}(\\mathbf{k})\\boldsymbol{\\mathsf{A(\\omega+\\nu,\\mathbf{k})}}
```

and the imaginary part have kramers-kronig relation of the real part

```math
\\text{Im}\\sigma_{\\alpha\\beta}(\\omega) = -\\frac{2\\omega}{\\pi}\\int_0^\\infty d\\omega' \\frac{\\text{Re}\\sigma_{\\alpha\\beta}(\\omega')}{\\omega^2 - \\omega'^2}
```

# Parameters
- `hwmesh` : Mesh of the optical frequency
- `hnw` : Mesh size of the optical frequency
- `T` : Temperature in Kelvin
- `dir` : Direction of the optical response

# Optional Arguments
- `report` : Boolean argument. It will printout the current task that this function
            calculate. Default is `false`

# Returns
- `hw` : optical frequency mesh
- `sigma` : complex optical conductivity

# Examples
```julia
# setup the model
tb = tb_model(2,2,lat,orb)
set_onsite! ...
set_hop! ...
# warmup the calculation
warmup_calc!(tb,[200,200],150,[-4,4])
# calculate optical conductivity
hw, sigma_yx = calc_OptCond(tb, (0,3), 100, T=50,dir=[2,1],report=true)
```

"""
function calc_OptCond(model::model, hwmesh,hnw; T,dir,report=false)
    if model.cache == 0
        @error "cache has not been setup. Please, warm up first."
    end

    function KK(ww,s)
        isigma = zeros(Float64,length(s))
        for (iw,w1) in enumerate(ww)
            wsq = ww.^2 .- w1^2
            intg = s./ wsq
            intg[iw] = 0
            isigma[iw] = integrate(ww,intg)
        end
        return -2ww .* isigma ./ π
    end

    μ = model.mu
    nw = length(model.cache.w)
    nk = length(model.cache.k)

    if report == true
        @info "Preparing conductivity calculation"
    end
    hw = LinRange(hwmesh...,hnw)
    dw = model.cache.w[2] - model.cache.w[1]

    function check(nw,w,calb)
        dw = w[2] - w[1]
        i1 = floor(Int,nw/2) + 20
        i2 = floor(Int,nw/2) - 30
        gt = w[i1] + w[i2]
        half = floor(Int,nw/2) + calb
        id = floor(Int,(w[i1]+w[i2])/dw) + half
        diff = (w[id] - gt) / dw
        return diff
    end

    half = 0
    for calb in -2:2
        if abs(check(nw,model.cache.w,calb)) < 0.51
            half = floor(Int,nw/2) + calb
            break
        end
    end

    kT = 8.617332478e-5 * T

    sigma = []

    if report == true
        @info "Calculating optical conductivity..."
    end
    for (iw,w) in enumerate(hw)
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
        sig = integrate(model.cache.w,tmp)
        push!(sigma,sig)
    end

    isigma = KK(hw,float.(sigma))

    return hw,sigma .+ 1im*isigma
end

"""
    calc_Epsilon(hw, sigma; eps0,eps_infty)

Calculate epsilon for given optical conductivity.

# Examples
```julia
# setup the model
tb = tb_model(2,2,lat,orb)
set_onsite! ...
set_hop! ...
# warmup the calculation
warmup_calc!(tb,[200,200],150,[-4,4])
# calculate optical conductivity
hw, sigma_xx = calc_OptCond(tb, (0,3), 100, T=300,dir=[1,1],report=true)
eps = calc_Epsilon(hw, sigma_xx, 1.0, 2.0)
```

"""
function calc_Epsilon(hw,sigma;eps0,eps_infty)
    return eps_infty .+ 1im .* sigma ./ (eps0 .* hw)
end

"""
    calc_Reflectance(epsilon)

Calculate reflectance for given epsilon.

# Examples
```julia
# setup the model
tb = tb_model(2,2,lat,orb)
set_onsite! ...
set_hop! ...
# warmup the calculation
warmup_calc!(tb,[200,200],150,[-4,4])
# calculate optical conductivity
hw, sigma_xx = calc_OptCond(tb, (0,3), 100, T=300,dir=[1,1],report=true)
eps = calc_Epsilon(hw, sigma_xx, 1.0, 2.0)
rw = calc_Reflectance(eps)
```

"""
function calc_Reflectance(epsilon)
    sqeps = sqrt.(epsilon)
    sqeps = (sqeps .- 1) ./ (sqeps .+ 1)
    refle = abs.(sqeps).^2
    return refle
end

"""
    calc_LossFunction(epsilon)

Calculate loss function for given epsilon.

```math
\\text{LF}(\\omega) = \\frac{\\epsilon_2^2(\\omega)}{\\epsilon_1^2(\\omega) + \\epsilon_2^2(\\omega)}
```

# Examples
```julia
# setup the model
tb = tb_model(2,2,lat,orb)
set_onsite! ...
set_hop! ...
# warmup the calculation
warmup_calc!(tb,[200,200],150,[-4,4])
# calculate optical conductivity
hw, sigma_xx = calc_OptCond(tb, (0,3), 100, T=300,dir=[1,1],report=true)
eps = calc_Epsilon(hw, sigma_xx, 1.0, 2.0)
LF = calc_LossFunction(eps)
```
"""
function calc_LossFunction(epsilon)
    imeps = imag(epsilon)
    return imeps ./ (real(epsilon).^2 .+ imeps)
end
