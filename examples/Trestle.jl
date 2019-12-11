using PythTB
using PyPlot

#define lattice vectors
lat = [[2.0,0.0],[0.0,1.0]]

#put orbital position
orb = [[0.0,0.0],[0.5,1.0]]

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
fig, ax = visualize(trestle,1,2)
ax.set_title("Trestle Lattice")
ax.set_xlabel("X")
ax.set_ylabel("Y")
fig.savefig("examples/Trestle Lattice.pdf")

#calculate energy
nk = 201
k_vec,k_dist,k_node = k_path(trestle,"fullc",nk)
eigvals = solve_eig(trestle,k_vec)

#plot band
fig, ax = subplots(figsize=(6,4))
ax.plot(k_dist,eigvals')

labels = [L"-\pi","0",L"\pi"]
ax.set_title("Trestle band structure")
ax.set_ylabel("Band Energy (eV)")
ax.set_xlabel("Path in K-space")
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
ax.set_xlim(k_node[1],k_node[end])
for k in k_node
    ax.axvline(x = k, linewidth = 0.5, color = "k")
end
fig.savefig("examples/Trestle_Band.pdf")
