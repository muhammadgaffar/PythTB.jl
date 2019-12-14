using PythTB
using PyPlot

# Calculates Berry phases for the Haldane model and compares it
# to the hybrid Wannier centers for a ribbon of the Haldane model.

# First, compute bulk Wannier centers along direction 1
# Then, cut a ribbon that extends along direcion 0, and compute
# both the edge states and the finite hybrid Wannier centers
# along direction 1.

#hexagonal lattice
lat = [[1,0],[0.5,sqrt(3)/2]]
#put orbital
orb = [[1/3,1/3],[2/3,2/3]]

#generate model
dim_k,dim_r = (2,2)
haldane = tb_model(dim_k,dim_r,lat,orb)

#model parameters
Δ = -0.2
t = -1.0
t2 = 0.05 - 0.15im
t2c = conj(t2)
efermi = 0.25

#set onsite
set_onsite!(haldane,[-Δ,Δ])
#set nearest neighbour hopping
set_hop!(haldane,t,1,2,[0,0])
set_hop!(haldane,t,2,1,[1,0])
set_hop!(haldane,t,2,1,[0,1])
#set second neighbour complex hopping
set_hop!(haldane,t2,1,1,[1,0])
set_hop!(haldane,t2,2,2,[1,-1])
set_hop!(haldane,t2,2,2,[0,1])
set_hop!(haldane,t2c,2,2,[1,0])
set_hop!(haldane,t2c,1,1,[1,-1])
set_hop!(haldane,t2c,1,1,[0,1])

#nk for kx and ky
nkx = 100
nky = 10

#compute bulk wannier center
setup_wf_array!(haldane,[nkx,nky])
solve_on_grid!(haldane,[0,0])
phi_1 = berry_phase(haldane,[1],2,contin=true)

#create haldane ribbon, cut along y
ribbon = cut_piece(haldane,nky,2,glue_edgs=false)
#generate kpoint along periodic direction (kx)
k_vec,k_dist,k_node = k_path(ribbon,"full",nkx)
k_label = [L"$0$",L"\pi",L"2\pi"]

#solve ribbon model
rib_eval,rib_evec = solve_eig(ribbon,k_vec,eig_vec=true)
#shift bands
rib_eval .-= efermi

#find k-points at which number of states below the Fermi level changes
jump_k = []
for i in 1:(size(rib_eval,2)-1)
    nocc_i = sum(rib_eval[:,i] .< 0)
    nocc_ip = sum(rib_eval[:,i+1] .< 0)
    if nocc_i != nocc_ip
        push!(jump_k,i)
    end
end

# plot expectation value of position operator for states in the ribbon
# and hybrid Wannier function centers
fig, (ax1, ax2) = subplots(2,1,figsize=(8,10))

# plot bandstructure of the ribbon
for n in 1:size(rib_eval,1)
  ax1.plot(k_dist,rib_eval[n,:],c="k", zorder=-50)
end

# color bands according to expectation value of y operator
# (red=top, blue=bottom)
for i in 1:size(rib_evec,2)
  global s
  # get expectation value of the position operator for states at i-th kpoint
  pos_exp = position_expectation(ribbon,rib_evec[:,i,:],2)
  # plot states according to the expectation value
  s = ax1.scatter(ones(size(rib_eval,1)) .* k_vec[i], rib_eval[:,i], c=pos_exp,
                s=30,marker="o", cmap="coolwarm", edgecolors="none", vmin=0.0,
                vmax=nky, zorder=-100)
end

#color scale
fig.colorbar(s,nothing,ax1,ticks=[0,nky])

# plot Fermi energy
ax1.axhline(0,c="g",zorder=-200)

# vertical lines show crossings of surface bands with Fermi energy
for ax in (ax1,ax2)
  for i in jump_k
    ax.axvline(x=(k_vec[i]+k_vec[i+1])/2.0, linewidth=0.7,
              color="k",zorder=-150)
  end
end


# tweaks
ax1.set_ylabel("Ribbon band energy")
ax1.set_ylim(-2.3,2.3)

# bottom plot shows Wannier center flow
# bulk Wannier centers in green lines
# finite-ribbon Wannier centers in black dots
# compare with Fig 3 in Phys. Rev. Lett. 102, 107603 (2009)

# plot bulk hybrid Wannier center positions and their periodic images
for j in -1:nky
    ax2.plot(k_vec,j .+ phi_1 ./ (2π),"k-",zorder=-50)
end

phi_1

# plot finite centers of ribbon along direction 1
for i in 1:size(rib_evec,2)
  # get occupied states only (those below Fermi level)
  occ_evec=rib_evec[rib_eval[:,i] .< 0.0,i,:]
  # get centers of hybrid wannier functions
  hwfc = position_hwf(ribbon,occ_evec,2)
  # plot centers
  s=ax2.scatter(ones(size(hwfc,1)).*k_vec[i], hwfc, c=hwfc, s=30,
                marker="o", cmap="coolwarm", edgecolors="none",
                vmin=0, vmax=nky, zorder=-100)
end

# color scale
fig.colorbar(s,nothing,ax2,ticks=[0,nky])

# tweaks
ax2.set_xlabel("k vector along direction 1")
ax2.set_ylabel("Wannier center along direction 2")
ax2.set_ylim(-0.5,nky+0.5)


# label both axes
for ax in (ax1,ax2)
  ax.set_xlim(k_node[1],k_node[end])
  ax.set_xticks(k_node)
  ax.set_xticklabels(k_label)
end
fig.savefig("examples/haldane_hwf.pdf")
