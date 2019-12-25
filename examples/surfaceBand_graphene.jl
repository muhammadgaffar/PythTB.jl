using PythTB
using Plots

# surface graphene using supercell

#define lattice vectors for honeycomb lattice
lat = [1 0; 0.5 sqrt(3)/2]

#put orbital
orb = [1/3 1/3; 2/3 2/3]

#graphene model
dim_k = 2
dim_r = 2
grp = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1
set_hop!(grp,t,1,2,[0,0])
set_hop!(grp,t,2,1,[0,1])
set_hop!(grp,t,2,1,[1,0])

#supercell graphene
sc_gr = make_supercell(grp,[[2,1],[-1,2]],to_home=true)

#cut model, finite along y direction
slab_gr = cut_piece(sc_gr,6,2,glue_edgs=false)

#visualize the finite model
fig = visualize_2d(slab_gr,1,2)
xlims!(fig,(-5,7))
ylims!(fig,(-1,12))

title!(fig,"Graphene, arbitrary surface")
xlabel!(fig,"x coordinate")
ylabel!(fig,"y coordinate")

#compute band structure
nk = 300
k_vec,k_dist,k_node = k_path(slab_gr,"full",nk)
eigvals = solve_eig(slab_gr,k_vec)

#plot
fig = plot(framestyle=:box,legend=false)
plot!(fig,k_dist,eigvals',c="black")
# zoom in close to the zero energy
xlims!(fig,(k_dist[1],k_dist[end]))
ylims!(fig,(-1.0,1.0))
# put title on top
title!(fig,"Graphene arbitrary surface band structure")
xlabel!(fig,"k parallel to edge")
ylabel!(fig,"Band energy")
xticks!(fig,k_node,["0","\\pi","2\\pi"])
