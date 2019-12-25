using PythTB
using PyPlot

#define lattice vectors
lat = [2 0; 0 1]

#put orbital position
orb = [0 0; 0.5 1]

#trestle model
dim_k = 1 # only 1 direction that periodic
dim_r = 2
trestle = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t1 = 0.8 + 0.6im
t2 = 2.0
set_hop!(trestle,t2,1,1,[1,0])
set_hop!(trestle,t2,2,2,[1,0])
set_hop!(trestle,t1,1,2,[0,0])
set_hop!(trestle,t1,2,1,[1,0])

#summary of the models
show(trestle)

#visualize the lattice
fig = visualize_2d(trestle,1,2)
title!(fig,"Trestle Lattice")
xlabel!(fig,"X")
ylabel!(fig,"Y")

#calculate energy
nk = 201
k_vec,k_dist,k_node = k_path(trestle,"fullc",nk)
eigvals = solve_eig(trestle,k_vec)

#plot band
fig = plot(framestyle=:box,legend=false)
plot!(fig,k_dist,eigvals')

labels = ["-\\pi","0","\\pi"]
title!(fig,"Trestle band structure")
ylabel!(fig,"Band Energy (eV)")
xlabel!(fig,"Path in K-space")
xticks!(fig,k_node,labels)
xlims!(fig,(k_node[1],k_node[end]))
vline!(fig,k_node,lw=0.4,alpha=0.5,color="black")
