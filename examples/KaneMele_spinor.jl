using PythTB
using PyPlot

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
    #get the model
    model = kane_mele(topo_index)
    #list of symmetry point in k-space
    path = [[0,0],[2/3,1/3],[0.5,0.5],[1/3,2/3],[0,0]]
    #labels
    labels = [L"$\Gamma$",L"$K$",L"$M$",L"$K^\prime$",L"$\Gamma$"]
    #genereate kpoints
    nk = 201
    k_vec,k_dist,k_node = k_path(model,path,nk)

    #solve the model
    eigvals = solve_eig(model,k_vec)

    #plot
    fig, (ax1, ax2) = subplots(1,2,figsize=(8,4))

    ax1.plot(k_dist,eigvals',"k-")
    ax1.set_title("Kane-Mele: "*topo_index*" phase")
    ax1.set_xticks(k_node)
    ax1.set_xticklabels(labels)
    ax1.set_xlim(k_node[1],k_node[end])

    for n in 1:length(k_node)
        ax1.axvline(x=k_node[n],linewidth=0.5, color="k")
    end
    ax1.set_xlabel("k-space")
    ax1.set_ylabel("Band Energy")

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
    for shift in -2:2
        ax2.plot(ky,wan_cent[:,1] .+ shift,"k.")
        ax2.plot(ky,wan_cent[:,2] .+ shift,"k.")
    end
    ax2.set_ylim(-1.0,1.0)
    ax2.set_ylabel("Wannier center along x")
    ax2.set_xlabel(L"$k_y$")
    ax2.set_xticks([0.0,0.5,1.0])
    ax2.set_xlim(0.0,1.0)
    ax2.set_xticklabels([L"$0$",L"$\pi$", L"$2\pi$"])
    ax2.axvline(x=0.5,linewidth=0.5, color="k")
    ax2.set_title("1D Wannier centers: "*topo_index*" phase")

    if topo_index == "even"
        fig.savefig("examples/kane_mele_even.pdf")
    else
        fig.savefig("examples/kane_mele_odd.pdf")
    end
end
