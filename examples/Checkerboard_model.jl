using PythTB
using Plots

#define two-dimensional lattice vectors
lat = [1 0; 0 1]

#put orbital positions for checkerboard
orb = [0.2 0.2; 0.7 0.7]

#TB model for checkerboard
dim_k = 2
dim_r = 2
CBoard = tb_model(2,2,lat,orb)

#set onsite energy
Δ = 1.1
set_onsite!(CBoard,[-Δ,Δ])

#set hopping
t = -0.6
set_hop!(CBoard,t,2,1,[0,0])
set_hop!(CBoard,t,2,1,[1,0])
set_hop!(CBoard,t,2,1,[0,1])
set_hop!(CBoard,t,2,1,[1,1])

#summary of the models
show(CBoard)

#visualize chekcerboard lattice (with hoppings)
fig = visualize_2d(CBoard,1,2)
title!(fig,"CheckBoard Lattice and Hoppings")
xlabel!(fig,"X")
ylabel!(fig,"Y")
fig

#choose k-points path in two-dimensional lattice
path = [0.0 0.0; 0.0 0.5; 0.5 0.5; 0.0 0.0]
nk = 301
k_vec,k_dist,k_node = k_path(CBoard,path,nk)
eigvals = solve_eig(CBoard,k_vec)

#plot band
fig2 = plot(framestyle=:box)
plot!(fig2,k_dist,eigvals')

labels = ["\\Gamma","X","M","\\Gamma"]
title!(fig2,"Checkerboard band structure")
ylabel!(fig2,"Band Energy (eV)")
xlabel!(fig2,"Path in K-space")
xticks!(k_node,labels)
xlims!(fig2,k_node[1],k_node[end])
vline!(k_node, linewidth = 0.5, linealpha = 0.4, c = "black")
