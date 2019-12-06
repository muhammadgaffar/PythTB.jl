using PythTB
using PyPlot

#define two-dimensional lattice vectors
lat = [[1.0,0.0],[0.0,1.0]]

#put orbital positions for checkerboard
orb = [[0.0,0.0],[0.5,0.5]]

#TB model for checkerboard
dim_k = 2
dim_r = 2
CBoard = tb_model(2,2,lat,orb)

#set onsite energy
Δ = 1.1
set_onsite!(CBoard,[-Δ,Δ])

#set hopping
t = -0.6
set_hop!(CBoard,t,2,1,[1,1])
set_hop!(CBoard,t,2,1,[2,1])
set_hop!(CBoard,t,2,1,[1,2])
set_hop!(CBoard,t,2,1,[2,2])

#summary of the models
show(CBoard)

#visualize chekcerboard lattice (with hoppings)
fig, ax = visualize(CBoard,1,2)
ax.set_title("CheckBoard Lattice and Hoppings")
ax.set_xlabel("X")
ax.set_ylabel("Y")
fig.savefig("examples/Checkerboard Lattice.pdf")

#choose k-points path in two-dimensional lattice
path = [[0.0,0.0],[0.0,0.5],[0.5,0.5],[0.0,0.0]]
nk = 301
k_vec,k_dist,k_node = k_path(CBoard,path,nk)
eigvals = solve_eig(CBoard,k_vec)

#plot band
fig, ax = subplots(figsize=(6,4))
ax.plot(k_dist,eigvals')

labels = [L"\Gamma","X","M",L"\Gamma"]
ax.set_title("Checkerboard band structure")
ax.set_ylabel("Band Energy (eV)")
ax.set_xlabel("Path in K-space")
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
ax.set_xlim(k_node[1],k_node[end])
for k in k_node
    ax.axvline(x = k, linewidth = 0.5, color = "k")
end
fig.savefig("examples/Checkerboard_Band.pdf")
