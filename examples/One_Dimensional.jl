using PythTB
using Plots

#define one-dimensional lattice vectors
lat = [[1.0]]

#put atom/orbital position
orb = [[0.0]]

#make new model
dim_k = 1 #dimension of k-space (this shows the periodicity)
dim_r = 1 #dimension of real-space
model1D = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1.0 #hopping amplitude
set_hop!(model1D,t,1,1,[1]) #orbital 1 to another orbital 1

#show summary of the model
show(model1D)

#generate k-points
nk = 201
k_vec,k_dist,k_node = k_path(model1D,"full",nk)

#evaluate energy over k
eigvals = solve_eig(model1D,k_vec)

#plot and save the figure
fig = plot(framestyle=:box,legend=false)
plot!(fig,k_dist,eigvals')

labels = ["0","\\pi","2\\pi"]
title!(fig,"Band Structure of One Dimensional Chain")
ylabel!(fig,"Band Energy (eV)")
xlabel!(fig,"Path in K-space")
xlims!(fig,k_node[1],k_node[end])
xticks!(fig,k_node,labels)
vline!(fig,k_node,alpha=0.5,lw=0.5,c="black")
