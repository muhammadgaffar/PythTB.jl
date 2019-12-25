using PythTB

#define lattice vectors for honeycomb lattice
lat = [1 0; 0.5 sqrt(3)/2]

#put orbital
orb = [1/3 1/3; 2/3 2/3]

#graphene model
dim_k = 2
dim_r = 2
gr = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1.0
set_hop!(gr,t,1,2,[0,0])
set_hop!(gr,t,2,1,[1,0])
set_hop!(gr,t,2,1,[0,1])

#model summary
show(gr)

#path in kspace (symmetry point)
path = [0 0;2/3 1/3; 0.5 0.5; 0 0]
labels = ["\\Gamma","K","M","\\Gamma"]

#solve band
nk = 201
k_vec,k_dist,k_node = k_path(gr,path,nk)
eigvals = solve_eig(gr,k_vec)

#plot band
fig = plot(framestyle=:box,legend=false)
plot!(fig,k_dist,eigvals')

title!(fig,"Graphene band structure")
ylabel!(fig,"Band Energy (eV)")
xlabel!(fig,"Path in K-space")
xticks!(fig,k_node,labels)
xlims!(fig,(k_node[1],k_node[end]))
vline!(k_node,linealpha=0.2, c="black")
