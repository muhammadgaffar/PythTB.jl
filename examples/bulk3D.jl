using PythTB
using PyPlot

# example of 3d tight binding
# BCC crystal

#define BCC lattice
lat = [ 1 0 0; 0 1 0; 0 0 1]

#put atom
orb = [0 0 0; 0.5 0.5 0.5]

#get the model
bcc = tb_model(3,3,lat,orb)

#set onsite
Δ = 0.5
set_onsite!(bcc,[-Δ,Δ])
#set hop
t = -1.0
set_hop!(bcc,t,1,2,[0,0,0])
set_hop!(bcc,t,1,2,[0,-1,0])
set_hop!(bcc,t,1,2,[1,0,0])
set_hop!(bcc,t,1,2,[1,-1,0])
set_hop!(bcc,t,1,2,[0,0,-1])
set_hop!(bcc,t,1,2,[0,-1,-1])
set_hop!(bcc,t,1,2,[1,0,-1])
set_hop!(bcc,t,1,2,[1,-1,-1])

#generate kpoint
path = [0 0 0; -1/2 1/2 1/2; 0 1/2 0; 0 0 0;
        1/4 1/4 1/4; 0 1/2 0; 1/4 1/4 1/4; -1/2 1/2 1/2]
labels = ["\\Gamma","H","N","\\Gamma","P","N","P","H"]
k_vec,k_dist,k_node = k_path(bcc,path,101)

#solve band
eigvals = solve_eig(bcc,k_vec)

#plot
fig = plot(framestyle=:box,legend=false)
plot!(fig,k_dist,eigvals',c="black")
title!(fig,"BandStructure for body-centered cubic")
ylabel!(fig,"Band energy")
xlabel!(fig,"Path in k-space")
xlims!(fig,k_node[1],k_node[end])
xticks!(k_node,labels)
for n in 1:length(k_node)
    vline!(fig,[k_node[n]],linealpha=0.4,c="black")
end
display(fig)

#get density of state
w,dos = calc_DOS(bcc,80,150,window=[-10,10])

#plot
fig2 = plot(legend=false)
plot!(fig2,w,dos)
title!(fig2,"DOS of BCC")
xlabel!(fig2,"Energy")
ylabel!(fig2,"DOS")
ylims!(fig2,(0,maximum(dos)))
