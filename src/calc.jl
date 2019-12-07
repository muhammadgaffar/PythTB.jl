using StatsBase

function calc_DOS(model::model,kmesh=75,wmesh=100;window=[-4,4],finite=false)

    meshEnergy = LinRange(window[1],window[2],wmesh+1)
    w = 0.5 * (meshEnergy[1:end-1] .+ meshEnergy[2:end])

    if finite == true
        eigvals = solve_eig(model)
    else

        kpts = []
        if model.dim_k == 1
            for i in 1:kmesh
                k = [i/kmesh]
                push!(kpts,k)
            end
        elseif model.dim_k == 2
            for i in 1:kmesh, j in 1:kmesh
                k = [i/kmesh, j/kmesh]
                push!(kpts,k)
            end
        elseif model.dim_k == 3
            for i in 1:kmesh, j in 1:kmesh, k in 1:kmesh
                k = [i/kmesh, j/kmesh, k/kmesh]
                push!(kpts,k)
            end
        end

        eigvals = solve_eig(model,kpts)
    end

    hist = fit(Histogram,vec(eigvals),meshEnergy)
    dos  = hist.weights / ( size(eigvals,1) * size(eigvals,2) )

    return w,dos
end
