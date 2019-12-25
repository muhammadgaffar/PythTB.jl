using PythTB
using Plots

# Two dimensional tight-binding 2D Kane-Mele model
# C.L. Kane and E.J. Mele, PRL 95, 146802 (2005) using Eq. 1

function kane_mele(topological)
    #define hexagonal lattice vector
    lat = [[1,0],[0.5,sqrt(3/2)]]
    #put orbital
    orb = [[1/3,1/3],[2/3,2/3]]
    #tb model
    model = tb_model(2,2,lat,orb,nspin=2)
    #topological phase
    onsite = 0
    if topological == "even"
        onsite = 2.5
    elseif topological == "odd"
        onsite = 1.0
    end
    #parameters
    t = 1
    spin_orb = 0.6*t*0.5
    rashba = 0.25*t
    #set onsite
    set_onsite!(model,[onsite,-1*onsite])
    #set hopping for spin dependent hop
    set_hop!(model,t,1,2,[0,0])
    set_hop!(model,t,1,2,[0,-1])
    set_hop!(model,t,1,2,[-1,0])
    #σz effect
    set_hop!(model,-1im*spin_orb*σ.z,1,1,[0,1])
    set_hop!(model, 1im*spin_orb*σ.z,1,1,[1,0])
    set_hop!(model,-1im*spin_orb*σ.z,1,1,[1,-1])
    set_hop!(model, 1im*spin_orb*σ.z,2,2,[0,1])
    set_hop!(model,-1im*spin_orb*σ.z,2,2,[1,0])
    set_hop!(model, 1im*spin_orb*σ.z,2,2,[1,-1])
    #rashba first neighbour hop σxdy - σydzx
    sq32 = sqrt(3)/2
    set_hop!(model,1im*rashba*( 0.5*σ.x-sq32*σ.y),1,2,[0,0],mode="add")
    set_hop!(model,1im*rashba*(-1.0*σ.x         ),1,2,[0,-1],mode="add")
    set_hop!(model,1im*rashba*( 0.5*σ.x+sq32*σ.y),1,2,[-1,0],mode="add")
    return model
end

for topo_index in ["even","odd"]
    global fig1_odd,fig2_odd,fig1_even,fig2_even
    #get the model
    model = kane_mele(topo_index)
    #list of symmetry point in k-space
    path = [[0,0],[2/3,1/3],[0.5,0.5],[1/3,2/3],[0,0]]
    #labels
    labels = ["\\Gamma","K","M","K'","\\Gamma"]
    #genereate kpoints
    nk = 201
    k_vec,k_dist,k_node = k_path(model,path,nk)

    #solve the model
    eigvals = solve_eig(model,k_vec)

    #plot
    fig1 = plot(framestyle=:box,legend=false)
    plot!(fig1,k_dist,eigvals',c="black")
    title!(fig1,"Kane-Mele: "*topo_index*" phase")
    xticks!(fig1,k_node,labels)
    xlims!(fig1,(k_node[1],k_node[end]))
    vline!(fig1,k_node,linewidth=0.5,alpha=0.5,c="black")
    xlabel!(fig1,"k-space")
    ylabel!(fig1,"Band Energy")

    #initialize wf array
    nk = 41
    setup_wf_array!(model,[nk,nk])
    #solve wf automatically by on grid
    #starting kpoint = -1/2,-1/2
    solve_on_grid!(model,[-0.5,-0.5])

    # calculate Berry phases around the BZ in the k_x direction
    # (which can be interpreted as the 1D hybrid Wannier centers
    # in the x direction) and plot results as a function of k_y
    #
    # Following the ideas in
    # A.A. Soluyanov and D. Vanderbilt, PRB 83, 235401 (2011)
    # R. Yu, X.L. Qi, A. Bernevig, Z. Fang and X. Dai, PRB 84, 075119 (2011)
    # the connectivity of these curves determines the Z2 index

    wan_cent = berry_phase(model,[1,2],2,contin=false,berry_evals=true)
    wan_cent ./= 2π

    ky = LinRange(0,1,nk)

    # draw shifted Wannier center positions
    fig2 = plot(framestyle=:box,legend=false)
    for shift in -2:2
        plot!(fig2,ky,wan_cent[:,1] .+ shift,marker=:circle,c="black")
        plot!(fig2,ky,wan_cent[:,2] .+ shift,marker=:circle,c="black")
    end
    ylims!(fig2,(-1.0,1.0))
    ylabel!(fig2,"Wannier center along x")
    xlabel!(fig2,"ky")
    labels = ["0","\\pi", "2\\pi"]
    xticks!(fig2,[0.0,0.5,1.0],labels)
    xlims!(fig2,(0.0,1.0))
    vline!(fig2,[0.5],linewidth=0.5,alpha=0.1, c="black")
    title!(fig2,"1D Wannier centers: "*topo_index*" phase")

    if topo_index == "odd"
        fig1_odd = fig1
        fig2_odd = fig2
    elseif topo_index == "even"
        fig1_even = fig1
        fig2_even = fig2
    end
end
