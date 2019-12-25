using PythTB

# Calculates Berry phases for the Haldane model and compares it
# to the hybrid Wannier centers for a ribbon of the Haldane model.

# First, compute bulk Wannier centers along direction 1
# Then, cut a ribbon that extends along direcion 0, and compute
# both the edge states and the finite hybrid Wannier centers
# along direction 1.

#hexagonal lattice
lat = [1 0; 0.5 sqrt(3)/2]
#put orbital
orb = [1/3 1/3; 2/3 2/3]

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
k_label = ["0","\\pi","2\\pi"]

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
fig = plot(framestyle=:box,legend=false)

# plot bandstructure of the ribbon
plot!(fig,k_dist,rib_eval',c="black")

# color bands according to expectation value of y operator
# (red=top, blue=bottom)
for i in 1:size(rib_evec,2)
  # get expectation value of the position operator for states at i-th kpoint
  pos_exp = position_expectation(ribbon,rib_evec[:,i,:],2)
  # plot states according to the expectation value
  scatter!(fig,ones(size(rib_eval,1)) .* k_vec[i], rib_eval[:,i],
          ms = 2, zcolor = pos_exp, m =(:coolwarm,0.3), alpha = 0.4, lab = "grad")
end


# plot Fermi energy
hline!(fig,[0],c="black",lw=0.2,ls=:dash)

# vertical lines show crossings of surface bands with Fermi energy
for i in jump_k
  vline!(fig,[(k_vec[i]+k_vec[i+1])/2.0],color="black",alpha=0.5)
end

#tweaks
xlims!(fig,(0,1))
ylabel!(fig,"Ribbon band energy")
ylims!(fig,(-2.3,2.3))
xlims!(fig,(0,1))
xticks!(fig,k_node,k_label)

# bottom plot shows Wannier center flow
# bulk Wannier centers in green lines
# finite-ribbon Wannier centers in black dots
# compare with Fig 3 in Phys. Rev. Lett. 102, 107603 (2009)

# plot bulk hybrid Wannier center positions and their periodic images
fig2 = plot(framestyle=:box,legend=false)
for j in -1:nky
    plot!(fig2,k_vec,j .+ phi_1 ./ (2π),c="black")
end

# plot finite centers of ribbon along direction 1
for i in 1:size(rib_evec,2)
  # get occupied states only (those below Fermi level)
  occ_evec=rib_evec[rib_eval[:,i] .< 0.0,i,:]
  # get centers of hybrid wannier functions
  hwfc = position_hwf(ribbon,occ_evec,2)
  # plot centers
  scatter!(fig2,ones(size(hwfc,1)).*k_vec[i],
          hwfc, zcolor=hwfc, ms=2,m=(:coolwarm,0.6))
end

for i in jump_k
  vline!(fig2,[(k_vec[i]+k_vec[i+1])/2.0],color="black",alpha=0.5)
end

# tweaks
xlabel!(fig2,"k vector along direction 1")
ylabel!(fig2,"Wannier center along direction 2")
ylims!(fig2,(-0.5,nky+0.5))
xlims!(fig2,(0,1))
xticks!(fig2,k_node,k_label)
