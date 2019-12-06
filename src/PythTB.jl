module PythTB
    using PyCall
    using StatsBase

    export TB,tb_model,set_hop!,set_onsite!,k_path,solve_eig, get_DOS
    export setup_wf_array!, berry_phase, berry_flux, impose_pbc!, cut_piece
    export position_expectation, solve_on_grid!, position_hwf, make_supercell
    export w90_model,dist_hop,bands_consistency
    export show,visualize

    const TB = PyNULL()

    function __init__()
        copy!(TB, pyimport_conda("pythtb", "pythtb"))
    end

    mutable struct model
        model
        array
        w90
        dim_k
        dim_r
        nspin
        per
        norb
        nstates
        latvec
        orb
        site_energies
        hoppings
    end

    include("visualize.jl")

    function tb_model(dim_k,dim_r,lat,orb; per = nothing,nspin = 1)
        if per == nothing
            tb = TB.tb_model(dim_k,dim_r,lat,orb,per = per,nspin = nspin)
        else
            tb = TB.tb_model(dim_k,dim_r,lat,orb,per = per .- 1,nspin = nspin)
        end
        tb = model(tb,0,false,0,0,0,0,0,0,0,0,0,0)
        tb.dim_k = tb.model._dim_k
        tb.dim_r = tb.model._dim_r
        tb.nspin = tb.model._nspin
        tb.per   = tb.model._per .+ 1
        tb.norb  = tb.model._norb
        tb.nstates = tb.model._nsta
        tb.latvec = lat
        tb.orb   = orb
        tb.site_energies = tb.model._site_energies
        tb.hoppings = tb.model._hoppings
        return tb
    end

    function set_hop!(model::model,hop_amp,i,j,R = nothing; mode = "set",
                    allow_conjugate_pair = true)

        model.model.set_hop(hop_amp,i-1,j-1,R .- 1,mode,allow_conjugate_pair)
        model.hoppings = model.model._hoppings
    end

    function set_onsite!(model::model,en,ind_i = nothing; mode = "set")
        model.model.set_onsite(en,ind_i,mode = mode)
        model.site_energies = model.model._site_energies
    end

    function k_path(model::model,kpts,nk; report = true)
        return model.model.k_path(kpts,nk,report)
    end

    function solve_eig(model::model,k_list = nothing; eig_vec = false)
        if k_list != nothing
            if length(k_list) <= 3
                model.model.solve_one(k_list, eig_vectors = eig_vec)
            else
                model.model.solve_all(k_list,eig_vectors = eig_vec)
            end
        else
            model.model.solve_all(k_list,eig_vectors = eig_vec)
        end
    end

    function get_DOS(model::model,kmesh=75,wmesh=100;window=[-4,4],finite=false)
        meshEnergy = LinRange(window[1],window[2],wmesh)
        w = 0.5 * (meshEnergy[1:end-1] .+ meshEnergy[2:end])

        kpts = []
        for i in 1:kmesh, j in 1:kmesh
            k = [i/kmesh, j/kmesh]
            push!(kpts,k)
        end

        if finite == true
            eigvals = solve_eig(model)
        else
            eigvals = solve_eig(model,kpts)
        end

        hist = fit(Histogram,vec(eigvals),meshEnergy)
        dos  = hist.weights / ( size(eigvals,1) * size(eigvals,2) )

        return w,dos
    end

    function setup_wf_array!(model::model,mesh)
        model.array = TB.wf_array(model.model,mesh)
    end

    function berry_phase(model::model,occ,dir=nothing;
                        contin=true,berry_evals=false)

        return model.array.berry_phase(occ,dir,contin,berry_evals)
    end

    function berry_flux(model::model,occ,dirs=nothing;
                        individual_phases=false)

        return model.array.berry_flux(occ,dirs,individual_phases = individual_phases)
    end

    function impose_pbc!(model::model,mesh_dir,k_dir)
        model.array.impose_pbc(mesh_dir,k_dir)
    end

    function cut_piece(mod::model,num,finite_dir;glue_edgs=false)
        finite_tb = mod.model.cut_piece(num,finite_dir-1,glue_edgs)
        finite_tb = model(finite_tb,0,false,0,0,0,0,0,0,0,0,0,0)
        finite_tb.dim_k = finite_tb.model._dim_k
        finite_tb.dim_r = finite_tb.model._dim_r
        finite_tb.nspin = finite_tb.model._nspin
        finite_tb.per   = finite_tb.model._per .+ 1
        finite_tb.norb  = finite_tb.model._norb
        finite_tb.nstates = finite_tb.model._nsta
        finite_tb.latvec = mod.latvec
        finite_tb.orb   = mod.orb
        finite_tb.site_energies = finite_tb.model._site_energies
        finite_tb.hoppings = finite_tb.model._hoppings
        return finite_tb
    end

    function position_expectation(model::model,eig_vectors,dir)
        return model.model.position_expectation(eig_vectors,dir-1)
    end

    function solve_on_grid!(model::model,start_k)
        model.array.solve_on_grid(start_k)
    end

    function position_hwf(model::model,evec,dir;hwf_evec=false,basis="orbital")
        return model.model.position_hwf(evec,dir,hwf_evec,basis)
    end

    function make_supercell(mod::model,sc_red_lat;return_sc_vectors=false,to_home=true)
        sc_tb = mod.model.make_supercell(sc_red_lat,return_sc_vectors,to_home)
        sc_tb = model(sc_tb,0,false,0,0,0,0,0,0,0,0,0,0)
        sc_tb.dim_k = sc_tb.model._dim_k
        sc_tb.dim_r = sc_tb.model._dim_r
        sc_tb.nspin = sc_tb.model._nspin
        sc_tb.per   = sc_tb.model._per .+ 1
        sc_tb.norb  = sc_tb.model._norb
        sc_tb.nstates = sc_tb.model._nsta
        sc_tb.latvec = mod.latvec
        sc_tb.orb   = mod.orb
        sc_tb.site_energies = sc_tb.model._site_energies
        sc_tb.hoppings = sc_tb.model._hoppings
        return sc_tb
    end

    function w90_model(path,prefix; zero_energy = 0.0, min_hopping_norm = nothing,
                    max_distance=nothing, ignorable_imaginary_part=nothing)

        w90_init = TB.w90(path,prefix)
        w90_tb = w90_init.model(zero_energy,min_hopping_norm,
                                max_distance,ignorable_imaginary_part)
        w90_tb = model(w90_tb,0,w90_init,0,0,0,0,0,0,0,0,0,0)
        return w90_tb
    end

    function dist_hop(model::model)
        dist, ham = model.w90.dist_hop()
        return dist, ham
    end

    function bands_consistency(model::model)
        kpts, energy = model.w90.w90_bands_consistency()
        return kpts, energy
    end

end # module PythTB
