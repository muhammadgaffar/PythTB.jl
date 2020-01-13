module PythTB
    using PyCall
    using Requires

    export TB,tb_model,set_hop!,set_onsite!,k_path,solve_eig, hamiltonian
    export val_to_block, remove_orb, k_uniform_mesh, position_matrix
    export setup_wf_array!, berry_phase, berry_flux, impose_pbc!, cut_piece
    export position_expectation, solve_on_grid!, position_hwf, make_supercell
    export impose_loop!, w90_model,dist_hop,bands_consistency
    export warmup_calc!, calc_DOS, calc_OptCond, calc_Epsilon, calc_Reflectance, calc_LossFunction
    export show
    export σ

    export KITE, configuration, calc_dos, calc_opticalConductivity, calc_dcConductivity

    ####################################
    # PythTB Wrapper and Plots dependencies
    ###################################

    const TB = PyNULL()
    const pb = PyNULL()

    function __init__()
        copy!(TB, pyimport_conda("pythtb", "pythtb"))
        copy!(pb, pyimport_conda("pybinding", "pybinding"))
        @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("visualize.jl")
    end

    ####################################
    # Base of the module
    ###################################

    import Base: show

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
        mu
        cache

        model(model) = new(model,0,false,
                            0,0,0,0,0,0,
                            0,0,0,0,0,0)
    end

    struct pauli
        x
        y
        z
    end
    const σx = [0,1,0,0]
    const σy = [0,0,1,0]
    const σz = [0,0,0,1]

    const σ = pauli(σx,σy,σz)

    function show(io::IO, model::model)
        println(io,"Tight-binding model summary")
        println(io,"===========================================================")
        println(io,"k-space dimension = ", model.dim_k)
        println(io,"r-space dimension = ", model.dim_r)
        println(io,"periodic direction = ", model.per)
        println(io,"number of orbitals = ", model.norb)
        println(io,"number of electronic states = ", model.nstates)
        println(io,"------------------------------------------")
        println(io,"lattice vectors:")
        for i in 1:size(model.latvec,1)
            o = model.latvec[i,:]
            println(io,"# â[$i] ==> ", round.(o,digits=4))
        end
        println(io,"------------------------------------------")
        println(io,"position of orbitals:")
        for i in 1:size(model.orb,1)
            o = model.orb[i,:]
            println(io,"# r[$i] ==> ", round.(o,digits=4))
        end
        println(io,"------------------------------------------")
        println(io,"site energies:")
        if model.nspin == 2
            for i in 1:size(model.site_energies,1)
                o = model.site_energies[i,:,:]
                println(io,"# ε[$i] : ")
                show(IOContext(io, :limit=>false), MIME"text/plain"(), o)
                println(io,"")
            end
        else
            for i in 1:size(model.site_energies,1)
                o = model.site_energies[i]
                println(io,"# ε[$i] ==> ", round.(o,digits=4))
            end
        end
        println(io,"------------------------------------------")
        println(io,"hoppings:")
        if model.hoppings == []
            println(io,"no hopping")
        elseif model.nspin == 2
            for i in 1:size(model.hoppings,1)
                o = model.hoppings[i,:]
                if length(o) == 4
                    println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)+$(o[4])⟩ : ")
                    show(IOContext(io, :limit=>false), MIME"text/plain"(), o[1])
                    println(io,"")
                    println(io,"")
                else
                    println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)⟩ : ")
                    show(IOContext(io, :limit=>false), MIME"text/plain"(), o[1])
                    println(io,"")
                    println(io,"")
                end
            end
        else
            for i in 1:size(model.hoppings,1)
                o = model.hoppings[i,:]
                if length(o) == 4
                    println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)+$(o[4])⟩ = ", round.(complex(o[1]),digits=4))
                else
                    println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)⟩ = ", round.(complex(o[1]),digits=4) )
                end
            end
        end

    end

    function show(io::IO,pauli::pauli)
        println(io,"Pauli matrices (x,y,z)")
    end

    ####################################
    # Real module functions begin here
    ###################################

    """
        tb_model(dim_k,dim_r,lat,orb,per,nspin)

    construct object that (will) contains all tightbinding informations of the model.
    parameter `dim_r` can be larger than `dim_k`

    # Optional Arguments
    - `per` : List of lattice vectors which are be considered to be periodic.
              By default, all lattice vectors are assumed to be periodic.
              If `dim_k` is smaller than `dim_r`, then by default the first
              dim_k vectors are considered to be periodic.
    - `nspin` : Number of explicit spin components assumed for each orbital in orb.
                Allowed values of `nspin` are 1 and 2. If `nspin` is 1 then the model is spinless,
                if `nspin` is 2 then it is explicitly a spinfull model and each orbital is
                assumed to have two spin components. Default value of this parameter is 1.

    # Examples
    ```julia
    # Creates model that is two-dimensional in real space but only
    # one-dimensional in reciprocal space. Second lattice vector is
    # chosen to be periodic (since per=[2]). Three orbital
    # coordinates are specified.

    lat = [1.0 0.5; 0.0 2.0]
    orb = [0.2 0.3; 0.1 0.1; 0.2 0.2]
    model  = tb_model(1, 2, lat, orb, per=[2])

    # read the summary of the model
    show(model)
    ```
    """
    function tb_model(dim_k,dim_r,lat,orb; per = nothing,nspin = 1)

        #spin must be 1 or 2
        if nspin != 1 && nspin != 2
            error("nspin must either 1 or 2")
        end

        if per == nothing
            tb = TB.tb_model(dim_k,dim_r,lat,orb,per = per,nspin = nspin)
        else
            tb = TB.tb_model(dim_k,dim_r,lat,orb,per = per .- 1,nspin = nspin)
        end
        tb = model(tb)
        tb.dim_k = tb.model._dim_k
        tb.dim_r = tb.model._dim_r
        tb.nspin = tb.model._nspin
        tb.per   = tb.model._per .+ 1
        tb.norb  = tb.model._norb
        tb.nstates = tb.model._nsta
        tb.latvec = tb.model._lat
        tb.orb   = tb.model._orb
        tb.site_energies = tb.model._site_energies
        tb.hoppings = tb.model._hoppings
        return tb
    end

    """
        set_hop!(model, hop_amp,ind_i,ind_j,ind_R,
                    mode = "set",allow_conjugate_pair = true)

    Defines hopping parameters between tight-binding orbitals of the TightBinding object `model`.

    This function specifies the following object
    ```math
    H_{ij}(\\mathbf{R}) = \\left\\langle \\phi_{\\mathbf{0},i}| H | \\phi_{\\mathbf{R},j} \\right\\rangle
    ```

    # Parameters
    - `hop_amp` :   Hopping amplitude; can be real or complex number.
                    If nspin is 2 then hopping amplitude can be given either as a single number,
                    or as an array of four numbers, or as 2x2 matrix. If a single number is given,
                    it is interpreted as hopping amplitude for both up and down spin component.
                    If an array of four numbers is given, these are the coefficients of I, sigma_x,
                    sigma_y, and sigma_z (that is, the 2x2 identity and the three Pauli spin matrices) respectively.
                    Finally, full 2x2 matrix can be given as well.
    - `ind_i` :     Index of bra orbital. This orbital is assumed to be in the home unit cell.
    - `ind_j` :     Index of ket orbital. This orbital does not have to be in the home unit cell;
                    its unit cell position is determined by parameter ind_R.
    - `ind_R` :     Lattice vector (integer array, in reduced coordinates) pointing to the unit cell where the ket orbital is located.

    # Optional arguments
    - `mode` :      Speficies way in which parameter hop_amp is used. Default is `set`.
                    It can either set value of hopping term from scratch (`set`),
                    reset it (`reset`), or add to it (`add`).
    - `allow_conjugate_pair` : Default is `True`. If set to `True` code will allow user to specify hopping i -> j+R
                             even if conjugate-pair hopping j -> i-R has been specified

    # Examples
    ```julia
    # Specifies complex hopping amplitude between first orbital in home
    # unit cell and third orbital in neigbouring unit cell.
    set_hop!(model, 0.3+0.4j, 1, 3, [0, 1])
    # change value of this hopping
    set_hop!(model, 0.1+0.2j, 1, 3, [0, 1], mode="reset")
    # add to previous value (after this function call below,
    # hopping term amplitude is 100.1+0.2j)
    set_hop!(model, 100.0, 1, 3, [0, 1], mode="add")
    ```
    """
    function set_hop!(model::model,hop_amp,i,j,R = nothing; mode = "set",
                    allow_conjugate_pair = true)

        model.model.set_hop(hop_amp,i-1,j-1,R,mode,allow_conjugate_pair)
        model.hoppings = model.model._hoppings
    end

    """
        set_onsite!(model, onsite_en,ind_i,mode = "set")

    Defines on-site energies for tight-binding orbitals. One can either set
    energy for one tight-binding orbital, or all at once.

    # Parameters
    - `onsite_en` : Either a list of on-site energies (in arbitrary units)
                    for each orbital, or a single on-site energy. In the case
                    when `nspin` is 1 (spinless) then each on-site energy is a
                    single number. If `nspin` is 2 then on-site energy can be given
                    either as a single number, or as an array of four numbers, or
                    2x2 matrix.
    - `ind_i` : Index of tight-binding orbital whose on-site energy you wish
                to change. This parameter should be specified only when `onsite_en`
                is a single number (not a list).

    # Optional arguments
    - `mode` :  Speficies way in which parameter onsite_en is used. Default is `set`.
                It can either set value of on-site energy from scratch, reset it,
                or add to it.

    # Examples
    ```julia
    # Defines on-site energy of first orbital to be 0.0, second 1.0, and third 2.0
    set_onsite!(model, [0.0, 1.0, 2.0])
    # Increases value of on-site energy for second orbital
    set_onsite!(model, 100.0, 1, mode="add")
    # Changes on-site energy of second orbital to zero
    set_onsite!(model, 0.0, 1, mode="reset")
    # Sets all three on-site energies at once
    set_onsite!(model, [2.0, 3.0, 4.0], mode="reset")
    ```
    """
    function set_onsite!(model::model,en,ind_i = nothing; mode = "set")
        if ind_i == nothing
            model.model.set_onsite(en,ind_i,mode = mode)
        else
            model.model.set_onsite(en,ind_i - 1,mode = mode)
        end
        model.site_energies = model.model._site_energies
    end

    function val_to_block(model::model,val)
        return model.model._val_to_block(val)
    end

    function remove_orb(mod::model,to_remove)
        tb = mod.model.remove_orb(to_remove .- 1)
        tb = model(tb)
        tb.dim_k = tb.model._dim_k
        tb.dim_r = tb.model._dim_r
        tb.nspin = tb.model._nspin
        tb.per   = tb.model._per .+ 1
        tb.norb  = tb.model._norb
        tb.nstates = tb.model._nsta
        tb.latvec = tb.model._lat
        tb.orb   = tb.model._orb
        tb.site_energies = tb.model._site_energies
        tb.hoppings = tb.model._hoppings
        return tb
    end

    function k_uniform_mesh(model::model,mesh_size)
        dim = model.dim_k
        k = []
        if dim == 1
            for i in 0:(mesh_size[1]-1)
                push!(k,[i/(mesh_size[1]-1)])
            end
        elseif dim == 2
            for i in 0:(mesh_size[1]-1), j in 0:(mesh_size[2]-1)
                push!(k,[i/(mesh_size[1]-1),j/(mesh_size[2]-1)])
            end
        elseif dim == 3
            for i in 0:(mesh_size[1]-1), j in 0:(mesh_size[2]-1), l in 0:(mesh_size[3]-1)
                push!(k,[i/(mesh_size[1]-1),j/(mesh_size[2]-1),l/(mesh_size[3]-1)])
            end
        end
        return k
    end

    """
        position_matrix(model, eigvec,dir)

    Returns matrix elements of the position operator along direction `dir` for
    eigenvectors `evec` at a single k-point. Position operator is defined in
    reduced coordinates.

    The returned object X is

    ```math
    X^\\alpha_{mn\\mathbf{k}} = \\left\\langle u_{m\\mathbf{k}} | r^\\alpha | u_{n\\mathbf{k}} \\right\\rangle
    ```

    # Examples
    ```julia
    # diagonalizes Hamiltonian at some k-points
    evals, evecs = solve_eig(model, k_vec,eig_vec=true)
    # computes position operator matrix elements for 3-rd kpoint
    # and bottom five bands along first coordinate
    pos_mat = position_matrix(model, evecs[1:5,2],1)
    ```
    """
    function position_matrix(model::model,eigvec,dir)
        return model.model.position_matrix(eigvec,dir .- 1)
    end

    """
        k_path(model, kpts,nk,report=true)

    Interpolates a path in reciprocal space between specified k-points. In 2D or
    3D the k-path can consist of several straight segments connecting high-symmetry
    points (“nodes”), and the results can be used to plot the bands along this path.

    # Optional arguments
    - `report` : Optional parameter specifying whether printout is desired
                (default is `false`).

    # Returns
    - `k_vec` : Array of (nearly) equidistant interpolated k-points. The distance between
                the points is calculated in the Cartesian frame, however coordinates
                themselves are given in dimensionless reduced coordinates! This
                is done so that this array can be directly passed to function `solve_eig`.
    - `k_dist` : Array giving accumulated k-distance to each k-point in the path.
                Unlike array k_vec this one has dimensions! (Units are defined here
                so that for an one-dimensional crystal with lattice constant equal
                to for example 10 the length of the Brillouin zone would equal 1/10=0.1.
                In other words factors of 2π are absorbed into k.) This array can be used
                to plot path in the k-space so that the distances between the
                k-points in the plot are exact.
    - `k_node` : Array giving accumulated k-distance to each node on the path in
                Cartesian coordinates. This array is typically used to plot nodes
                (typically special points) on the path in k-space.

    # Examples
    ```julia
    # Construct a path connecting four nodal points in k-space
    # Path will contain 401 k-points, roughly equally spaced
    path = [0.0 0.0; 0.0 0.5; 0.5 0.5; 0.0 0.0]
    k_vec,k_dist,k_node = k_path(model, path,401)
    # solve for eigenvalues on that path
    evals = solve_eig(model, k_vec)
    ```
    """
    function k_path(model::model,kpts,nk; report = false)
        return model.model.k_path(kpts,nk,report)
    end

    """
        hamiltonian(model, k)

    Generate hamiltonian for a model at certain of k-point. k-point is given in
    reduced coordinate.

    # Examples
    ```julia
    # for 3D k-space
    k_point = [1 π 1/3π]
    Hmat = hamiltonian(model, k_point)
    ```
    """
    function hamiltonian(model::model,k)
        return model.model._gen_ham(k)
    end

    """
        solve_eig(model, k_list=nothing,eig_vec=false)

    Solves for eigenvalues and (optionally) eigenvectors of the tight-binding
    model on a given one-dimensional list of k-vectors.

    # Examples
    ```julia
    # Returns eigenvalues for three k-vectors
    eval = solve_eig(model, [0.0 0.0; 0.0 0.2; 0.0 0.5])
    # Returns eigenvalues and eigenvectors for two k-vectors
    eval, evec = solve_eig(model, [0.0 0.0; 0.0 0.2],eig_vec=true)
    ```
    """
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

    """
        setup_wf_array!(model, mesh)

    This function is pre-process to construct grid in k-space. It will be used
    to perform calculation such Berry Phase, Berry curvature, 1st chern number, etc.

    If `model.array` is used for closed paths, either in a reciprocal-space or parametric direction,
    then one needs to include both the starting and ending eigenfunctions even though
    they are physically equivalent.  If the array dimension in question is a k-vector
    direction and the path traverses the Brillouin zone in a primitive reciprocal-lattice
    direction, `impose_pbc!` can be used to associate the starting and ending points
    with each other; if it is a non-winding loop in k-space or a loop in parameter space,
    then `impose_loop!` can be used instead.

    # Examples
    ```julia
    # Construct wf_array capable of storing an 11x21 array of wavefunctions
    setup_wf_array!(model, [11,21])
    # populate this wf_array with regular grid of points in Brillouin zone
    solve_on_grid!(model, [0.0, 0.0])
    # Compute set of eigenvectors at one k-point
    eval, evec = solve_eig([kx, ky], eig_vec = true)
    # Store it manually into a specified location in the array
    model.array[3, 4] = evec
    ```
    """
    function setup_wf_array!(model::model,mesh)
        model.array = TB.wf_array(model.model,mesh)
    end

    """
        berry_phase(model, occ,dir,contin=true,berry_evals=false)

    Computes the Berry phase along a given array direction and
    for a given set of occupied states.  This assumes that the
    occupied bands are well separated in energy from unoccupied
    bands. It is the responsibility of the user to check that
    this is satisfied.  By default, the Berry phase traced over
    occupied bands is returned, but optionally the individual
    phases of the eigenvalues of the global unitary rotation
    matrix (corresponding to "maximally localized Wannier
    centers" or "Wilson loop eigenvalues") can be requested
    (see parameter `berry_evals` for more details).

    For an array of size *N* in direction `dir`, the Berry phase
    is computed from the *N-1* inner products of neighboring
    eigenfunctions.  This corresponds to an "open-path Berry
    phase" if the first and last points have no special
    relation.  If they correspond to the same physical
    Hamiltonian, and have been properly aligned in phase using
    `impose_pbc!` or `impose_loop!`, then a closed-path Berry phase
    will be computed.

    For a one-dimensional `model.array` (i.e., a single string), the
    computed Berry phases are always chosen to be between -π and π.
    For a higher dimensional wf_array, the Berry phase is computed
    for each one-dimensional string of points, and an array of
    Berry phases is returned. The Berry phase for the first string
    (with lowest index) is always constrained to be between -π and
    π. The range of the remaining phases depends on the value of
    the input parameter `contin`.

    # Parameters
    - `occ` : Array of indices of energy bands which are considered to be occupied.

    # Optional arguments
    - `dir` : Index of `model.array` direction along which Berry phase is computed.
            This parameters needs not be specified for a one-dimensional `model.array`.
    - `contin` : Optional boolean parameter. If `true` then the branch choice of
                the Berry phase (which is indeterminate modulo 2π) is made so
                that neighboring strings (in the direction of increasing index value)
                have as close as possible phases. The phase of the first string
                (with lowest index) is always constrained to be between -π and π.
                If `false`, the Berry phase for every string is constrained to be
                between -π and π. The default value is `true`.
    - `berry_evals` : Optional boolean parameter. If `true` then will compute and
                return the phases of the eigenvalues of the product of overlap matrices.
                (These numbers correspond also to hybrid Wannier function centers.)
                These phases are either forced to be between -π and π (if `contin`* is
                `false`) or they are made to be continuous (if `contin` is `true`).

    # Return
    - `phase` : If `berry_evals` is `false` (default value) then
            returns the Berry phase for each string. For a
            one-dimensional `model.array` this is just one number. For a
            higher-dimensional `model.array` `phase` contains one phase for
            each one-dimensional string in the following format. For
            example, if `model.array` contains k-points on mesh with
            indices [i,j,k] and if direction along which Berry phase
            is computed is `dir=1` then `phase` will be two dimensional
            array with indices [i,k], since Berry phase is computed
            along second direction. If `berry_evals` is `true` then for
            each string returns phases of all eigenvalues of the
            product of overlap matrices.

    # Examples
    ```julia
    # Computes Berry phases along second direction for three lowest
    # occupied states. For example, if wf is threedimensional, then
    # phase[2,3] would correspond to Berry phase of string of states
    # along model.array[2,:,3]
    phase = berry_phase!(model, [1, 2, 3], 2)
    ```
    """
    function berry_phase(model::model,occ,dir=nothing;
                        contin=true,berry_evals=false)

        occ = Int.(occ)
        if dir == nothing
            return model.array.berry_phase(occ .-1,dir,contin,berry_evals)
        else
            dir = Int.(dir)
            return model.array.berry_phase(occ .-1,dir .-1,contin,berry_evals)
        end
    end

    """
        berry_flux(model, occ,dir,individual_phases=false)

    In the case of a 2-dimensional `model.array` array calculates the
    integral of Berry curvature over the entire plane.  In higher
    dimensional case (3 or 4) it will compute integrated curvature
    over all 2-dimensional slices of a higher-dimensional
    `model.array`

    # Parameters
    - `occ` : Array of indices of energy bands which are considered to be occupied.

    # Optional arguments
    - `dirs` : Array of indices of two `model.array` directions on which
            the Berry flux is computed. This parameter needs not be
            specified for a two-dimensional wf_array.  By default `dirs` takes
            first two directions in the array.
    - `individual_phases` : If `true` then returns Berry phase
                for each plaquette (small square) in the array. Default
                value is `false`.

    # Returns
    - `flux` : - In a 2-dimensional case returns and integral
            of Berry curvature (if `individual_phases` is `true` then
            returns integral of Berry phase around each plaquette).
            In higher dimensional case returns integral of Berry
            curvature over all slices defined with directions `dirs`.
            Returned value is an array over the remaining indices of
            `model.array`.  (If `individual_phases` is `true` then it
            returns again phases around each plaquette for each
            slice. First indices define the slice, last two indices
            index the plaquette.)

    # Examples
    ```julia
    flux  = berry_flux(model, [1,2,3])
    ```
    """
    function berry_flux(model::model,occ,dir=nothing;
                        individual_phases=false)
        occ = Int.(occ)
        if dir == nothing
            return model.array.berry_flux(occ .-1,dir,individual_phases = individual_phases)
        else
            dir = Int.(dir)
            return model.array.berry_flux(occ .-1,dir .-1,individual_phases = individual_phases)
        end

    end

    """
        impose_pbc!(model, mesh_dir,k_dir)

    This function will impose these periodic boundary conditions along
    one direction of the array. We are assuming that the k-point
    mesh increases by exactly one reciprocal lattice vector along
    this direction. This is currently not checked by the code;
    it is the responsibility of the user. Currently `model.array`
    does not store the k-vectors on which the model was solved;
    it only stores the eigenvectors (wavefunctions).

    # Examples
    ```julia
    # Imposes periodic boundary conditions along the mesh_dir=1
    # direction of the wf_array object, assuming that along that
    # direction the k_dir=2 component of the k-vector is increased
    # by one reciprocal lattice vector.  This could happen, for
    # example, if the underlying tb_model is two dimensional but
    # wf_array is a one-dimensional path along k_y direction.
    impose_pbc!(model, 1,2)
    ```
    """
    function impose_pbc!(model::model,mesh_dir,k_dir)
        model.array.impose_pbc(mesh_dir .-1,k_dir .-1)
    end

    """
        impose_loop!(model, mesh_dir)

    If the user knows that the first and last points along the
    `mesh_dir` direction correspond to the same Hamiltonian (this
    is not checked), then this routine can be used to set the
    eigenvectors equal (with equal phase), by replacing the last
    eigenvector with the first one (for each band, and for each
    other mesh direction, if any).

    This routine should not be used if the first and last points
    are related by a reciprocal lattice vector; in that case,
    `impose_pbc!` should be used instead.

    # Examples
    ```julia
    # Suppose the wf_array object is three-dimensional
    # corresponding to (kx,ky,lambda) where (kx,ky) are
    # wavevectors of a 2D insulator and lambda is an
    # adiabatic parameter that goes around a closed loop.
    # Then to insure that the states at the ends of the lambda
    # path are equal (with equal phase) in preparation for
    # computing Berry phases in lambda for given (kx,ky),
    impose_loop!(model, mesh_dir=3)
    ```
    """
    function impose_loop!(model::model, mesh_dir)
        model.array.impose_loop(mesh_dir .- 1)
    end

    """
        finite_model = cut_piece(model, num,finite_dir,glue_edgs=false)

    Constructs a (d-1)-dimensional tight-binding model out of a
    d-dimensional one by repeating the unit cell a given number of
    times along one of the periodic lattice vectors. The real-space
    lattice vectors of the returned model are the same as those of
    the original model; only the dimensionality of reciprocal space
    is reduced.

    # Examples
    ```julia
    A = tb_model(3, 3, ...)
    # Construct TB model out of model A
    # by repeating model along second lattice vector ten times
    B = cut_piece(A, 10, 2)
    # Further cut two-dimensional model B into one-dimensional model
    # A by repeating unit cell twenty times along third lattice
    # vector and allow hoppings from one edge to the other
    C = cut_piece(B, 20, 2, glue_edgs=true)
    ```
    """
    function cut_piece(mod::model,num,finite_dir;glue_edgs=false)
        finite_tb = mod.model.cut_piece(num,finite_dir-1,glue_edgs)
        finite_tb = model(finite_tb)
        finite_tb.dim_k = finite_tb.model._dim_k
        finite_tb.dim_r = finite_tb.model._dim_r
        finite_tb.nspin = finite_tb.model._nspin
        finite_tb.per   = finite_tb.model._per .+ 1
        finite_tb.norb  = finite_tb.model._norb
        finite_tb.nstates = finite_tb.model._nsta
        finite_tb.latvec = finite_tb.model._lat
        finite_tb.orb   = finite_tb.model._orb
        finite_tb.site_energies = finite_tb.model._site_energies
        finite_tb.hoppings = finite_tb.model._hoppings
        return finite_tb
    end

    """
        position_expectation(model, eigvec,dir)

    Returns diagonal matrix elements of the position operator.
    These elements can be interpreted as an average position of
    n-th Bloch state `eigvec[n]` along direction `dir`.

    # Examples
    ```julia
    # diagonalizes Hamiltonian at some k-points
    evals, evecs = solve_eig(mode, k_vec,eig_vec=true)
    # computes average position for 3-rd kpoint
    # and bottom five bands along first coordinate
    pos_exp = my_model.position_expectation(evecs[1:5,2], 1)
    ```
    """
    function position_expectation(model::model,eig_vectors,dir)
        return model.model.position_expectation(eig_vectors,dir-1)
    end

    """
        solve_on_grid!(model, start_k)

    Solve a tight-binding model on a regular mesh of k-points covering
    the entire reciprocal-space unit cell. Both points at the opposite
    sides of reciprocal-space unit cell are included in the array.

    # Examples
    ```julia
    solve_on_grid(model, [-0.5, -0.5])
    ```
    """
    function solve_on_grid!(model::model,start_k)
        model.array.solve_on_grid(start_k)
    end

    """
        position_hwf(model, eigvec,dir,hwf_evec=false,basis="orbital")

    Returns eigenvalues and optionally eigenvectors of the
    position operator matrix in either Bloch or orbital basis.
    These eigenvectors can be interpreted as linear combinations of
    Bloch states `eigvec` that have minimal extent along direction `dir`.
    The eigenvalues are average positions of these localized states.

    Note that these eigenvectors are not maximally localized
    Wannier functions in the usual sense because they are
    localized only along one direction.  They are also not the
    average positions of the Bloch states `eigvec`, which are
    instead computed by `position expectation`.

    # Optional arguments
    - `hwf_evec` : Boolean variable. If set to `true`
                this function will return not only eigenvalues but also
                eigenvectors of position operator. Default value is `false`.
    - `basis` : If `basis="bloch"` then hybrid Wannier function `hwf_evec` is written
            in the Bloch basis.  I.e. hwf[i,j] correspond to the weight of j-th
            Bloch state from `eigvec` in the i-th hybrid Wannier function.
            If `basis="orbital"` and `nspin=1` then hwf[i,orb] correspond to the
            weight of orb-th orbital in the i-th hybrid Wannier function.
            If `basis="orbital"` and `nspin=2` then hwf[i,orb,spin] correspond to the
            weight of orb-th orbital, spin-th spin component in the i-th hybrid
            Wannier function.  Default value is "orbital".

    # Examples
    ```julia
    # diagonalizes Hamiltonian at some k-points
    evals, evecs = solve_eig(model, k_vec,eig_vec=true)
    # computes hybrid Wannier centers (and functions) for 3-rd kpoint
    # and bottom five bands along first coordinate
    hwfc, hwf = position_hwf(model, evecs[1:5,2], 1, hwf_evec=true)
    ```
    """
    function position_hwf(model::model,evec,dir;hwf_evec=false,basis="orbital")
        return model.model.position_hwf(evec,dir -1,hwf_evec,basis)
    end

    """
        make_supercell(model, sc_red_lat, return_sc_vectors=false, to_home=true)

    Returns tight-binding model. Representing a super-cell of a current object.
    This function can be used together with `cut_piece` in order to create slabs
    with arbitrary surfaces.

    # Optional arguments
    - `return_sc_vectors` : Default value is `false`. If `true` returns
                        also lattice vectors inside the super-cell.
    - `to_home` : If `true` will shift all orbitals to the home cell.
                Default value is `true`.

    # Examples
    ```julia
    # Creates super-cell out of 2d tight-binding model tb
    sc_tb = make_supercell(model, [2 1; -1 2])
    ```
    """
    function make_supercell(mod::model,sc_red_lat;return_sc_vectors=false,to_home=true)
        sc_tb = mod.model.make_supercell(sc_red_lat,return_sc_vectors,to_home)
        sc_tb = model(sc_tb)
        sc_tb.dim_k = sc_tb.model._dim_k
        sc_tb.dim_r = sc_tb.model._dim_r
        sc_tb.nspin = sc_tb.model._nspin
        sc_tb.per   = sc_tb.model._per .+ 1
        sc_tb.norb  = sc_tb.model._norb
        sc_tb.nstates = sc_tb.model._nsta
        sc_tb.latvec = sc_tb.model._lat
        sc_tb.orb   = sc_tb.model._orb
        sc_tb.site_energies = sc_tb.model._site_energies
        sc_tb.hoppings = sc_tb.model._hoppings
        return sc_tb
    end

    """
        w90_model(path,prefix,zero_energy,min_hopping_norm,
                max_distance,ignorable_imaginary_part)

    This function imports tight-binding model parameters from an output
    of a Wannier90 code.

    The interface from Wannier90 to PythTB will use only the following
    files created by Wannier90:
    - `prefix.win`
    - `prefix_hr.dat`
    - `prefix_centres.xyz`
    - `prefix_band.kpt` (optional)
    - `prefix_band.dat` (optional)

    Warning : So far PythTB assumes that the position operator is diagonal in the tight-binding basis.

    # Optional arguments
    - `zero_energy` : Sets the zero of the energy in the band
                structure.  This value is typically set to the Fermi level
                computed by the density-functional code (or to the top of the
                valence band).  Units are electron-volts.
    - `min_hopping_norm` : Hopping terms read from Wannier90 with
                complex norm less than `min_hopping_norm` will not be included
                in the returned tight-binding model.
    - `max_distance` : Hopping terms from site i to site j+R will
                be ignored if the distance from orbital i to j+R is larger
                than `max_distance`.  This parameter is given in Angstroms.
    - `ignorable_imaginary_part` : The hopping term will be assumed to
            be exactly real if the absolute value of the imaginary part as
            computed by Wannier90 is less than `ignorable_imaginary_part`.

    # Examples
    ```julia
    # reads Wannier90 from folder called `example_a` folder
    # it assumes that that folder contains files "silicon.win" and so on
    w90model = w90("example_a", "silicon")
    ```
    """
    function w90_model(path,prefix; zero_energy = 0.0, min_hopping_norm = nothing,
                    max_distance=nothing, ignorable_imaginary_part=nothing)

        w90_init = TB.w90(path,prefix)
        w90_tb = w90_init.model(zero_energy,min_hopping_norm,
                                max_distance,ignorable_imaginary_part)
        w90_tb = model(w90_tb)
        w90_tb.w90 = w90_init
        w90_tb.dim_k = w90_tb.model._dim_k
        w90_tb.dim_r = w90_tb.model._dim_r
        w90_tb.nspin = w90_tb.model._nspin
        w90_tb.per   = w90_tb.model._per .+ 1
        w90_tb.norb  = w90_tb.model._norb
        w90_tb.nstates = w90_tb.model._nsta
        w90_tb.latvec = w90_tb.model._lat
        w90_tb.orb   = w90_tb.model._orb
        w90_tb.site_energies = w90_tb.model._site_energies
        w90_tb.hoppings = w90_tb.model._hoppings
        return w90_tb
    end

    """
        dist_hop(model)

    This is one of the diagnostic tools that can be used to help in
    determining `min_hopping_norm` and `max_distance` parameter in `w90_model` function.

    # Examples
    ```julia
    # get distances and hopping terms
    dist, ham = dist_hop(w90model)
    ```
    """
    function dist_hop(model::model)
        dist, ham = model.w90.dist_hop()
        return dist, ham
    end

    """
        bands_consistency(model)

    This function reads in band structure as interpolated by
    Wannier90.  Please note that this is not the same as the band
    structure calculated by the underlying DFT code.  The two will
    agree only on the coarse set of k-points that were used in
    Wannier90 generation.

    The purpose of this function is to compare the interpolation
    in Wannier90 with that in PythTB.

    The code assumes that the following files were generated by Wannier90,
    - `prefix_band.kpt`
    - `prefix_band.dat`

    # Examples
    ```julia
    # get band structure from wannier90
    w90_kpt,w90_evals = bands_consistency(w90model)
    # solve simplified model on the same k-path as in wannier90
    evals = solve_eig(w90model, w90_kpt)
    ```
    """
    function bands_consistency(model::model)
        kpts, energy = model.w90.w90_bands_consistency()
        return kpts, energy
    end

    include("calc.jl")
    include("kite.jl")

end # module PythTB
