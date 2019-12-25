using PythTB

#Haldane model is graphene model with complex second neighbour hopping

#define lattice vectors for honeycomb lattice
lat = [1 0; 0.5 sqrt(3)/2]
#put orbital
orb = [1/3 1/3; 2/3 2/3]

#generate model
dim_k,dim_r = (2,2)
haldane = tb_model(dim_k,dim_r,lat,orb)

#model parameters
Δ = 0.2
t = -1.0
t2 = 0.15 * exp(1im*π/2)
t2c = conj(t2)

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

#show model
show(haldane)

#visualize the model
using Plots
fig = visualize_2d(haldane,1,2)

#generate k point along symmetry point
path = [0 0; 2/3 1/3; 0.5 0.5; 1/3 2/3; 0 0]
labels = ["\\Gamma","K","M","K'","\\Gamma"]

#kpath
nk = 201
k_vec,k_dist,k_node = k_path(haldane,path,nk)

#solve eigenvalue
eigvals = solve_eig(haldane,k_vec)

#plot bandstructure
fig = plot(framestyle=:box,legend=false)
plot!(fig,k_dist,eigvals')
xlims!(fig,(k_node[1],k_node[end]))
xticks!(k_node,labels)
vline!(k_node,linewidth=0.4,linealpha=0.5,c="black")
title!(fig,"Haldane model band structure")
xlabel!(fig,"Path in k-space")
ylabel!(fig,"Band energy")

#Density of state
nk = 500
wmesh = 175
w,DOS = calc_DOS(haldane,nk,wmesh)

#DOS plot
fig = plot()
plot!(fig,w,DOS)
ylims!(fig,(0,maximum(DOS)))
title!(fig,"Haldane model density of states")
xlabel!(fig,"Band energy")
ylabel!(fig,"Density of state")
