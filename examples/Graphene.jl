using PythTB

#define lattice vectors for honeycomb lattice
lat = [[1,0],[0.5,sqrt(3)/2]]

#put orbital
orb = [[1/3,1/3],[2/3,2/3]]

#graphene model
dim_k = 2
dim_r = 2
gr = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1.0
set_hop!(gr,t,1,2,[1,1])
set_hop!(gr,t,2,1,[2,1])
set_hop!(gr,t,2,1,[1,2])

#model summary
show(gr)

#path in kspace (symmetry point)
path = [[0,0],[2/3,1/3],[0.5,0.5],[0,0]]
labels = [L"\Gamma","K","M",L"\Gamma"]

#solve band
nk = 201
k_vec,k_dist,k_node = k_path(gr,path,nk)
eigvals = solve_eig(gr,k_vec)

#plot band
fig, ax = subplots(figsize=(6,4))
ax.plot(k_dist,eigvals')

ax.set_title("Graphene band structure")
ax.set_ylabel("Band Energy (eV)")
ax.set_xlabel("Path in K-space")
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
ax.set_xlim(k_node[1],k_node[end])
for k in k_node
    ax.axvline(x = k, linewidth = 0.5, color = "k")
end
fig.savefig("examples/Graphene_Band.pdf")
