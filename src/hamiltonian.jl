using LinearAlgebra

function gen_ham(tb::PythTB.model, kpoint::AbstractVector)
    if length(kpoint) != tb.dim_k
        throw("Dimension Mismatch, dimension must same as dim_k")
    end

    # initialize hamiltonian array
    if tb.nspin == 1
        ham = zeros(ComplexF64, tb.norb,tb.norb)
    else
        ham = zeros(ComplexF64, tb.norb,2,tb.norb,2)
    end

    # fill diagonal matrices
    @inbounds for i in 1:tb.norb
        if tb.nspin == 1
            ham[i,i] = view(tb.site_energies,i,:)
        else
            ham[i,:,i,:] = view(tb.site_energies, i,:,:)
        end
    end

    # fill hoppings
    @inbounds Threads.@threads for ih in axes(tb.hoppings,1)
        h = view(tb.hoppings, ih,:)
        i,j = h[2]+1, h[3]+1
        if tb.dim_k > 0
            # rv = R_j - R_i
            rv = tb.orb[j,tb.per] + h[4] - tb.orb[i,tb.per]
            # calculate phase of hoppings
            amp = h[1] .* exp(2im*π*(kpoint' * rv))
        end
        # important caveat!
        # this section below take most time and allocation
        # must be optimized with some low-level trick
        # on array operations.
        if tb.nspin == 1
            ham[i,j] .+= amp
            ham[j,i] .+= conj(amp)
        else
            ham[i,:,j,:] .+= amp
            ham[j,:,i,:] .+= conj(transpose(amp))
        end
    end

    return Hermitian(reshape(ham,tb.nstates,tb.nstates))
end

function solve(tb, kpoints; eig_vec=false)
    if ndims(kpoints) == 1
        kpoints .= kpoints'
    end
    nk = size(kpoints,1)
    # initialize array
    eval = zeros(Float64, nk,tb.nstates)
    if eig_vec == true
        if tb.nspin == 1
            evec = zeros(ComplexF64, nk,tb.nstates,tb.norb)
        else
            evec = zeros(ComplexF64, nk,tb.nstates,tb.norb,2)
        end
    end

    # solve!
    @inbounds for ik in 1:nk
        if eig_vec == false
            eval[ik,:] .= eigvals(gen_ham(tb,kpoints[ik,:]))
        else
            F = eigen(gen_ham(tb,kpoints[ik,:]))
            # the corresponding eigvec is v[:,i] for eigval w[i]
            # transpose the eigvec so the corresponding index is
            # v[i,:]
            F.vectors .= transpose(F.vectors)
            eval[ik,:] .= F.values
            if tb.nspin == 1
                evec[ik,:,:] .= F.vectors
            else
                evec[ik,:,:,:] .= reshape(F.vectors, tb.nstates,tb.norb,2)
            end
        end
    end

    # return
    if eig_vec == false
        return eval
    else
        return eval,evec
    end
end
