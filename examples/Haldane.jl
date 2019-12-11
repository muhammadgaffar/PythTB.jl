using PythTB
using PyPlot

#Haldane model is graphene model with complex second neighbour hopping

#hexagonal lattice
lat = [[1,0],[0.5,sqrt(3)/2]]
#put orbital
orb = [[1/3,1/3],[2/3,2/3]]

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
fig, ax = visualize(haldane,1,2)
fig.savefig("examples/haldane_lattice.pdf")

#generate k point along symmetry point
path = [[0,0],[2/3,1/3],[0.5,0.5],[1/3,2/3],[0,0]]
labels = [L"\Gamma","K","M",L"K^\prime",L"\Gamma"]

#kpath
nk = 201
k_vec,k_dist,k_node = k_path(haldane,path,nk)

#solve eigenvalue
eigvals = solve_eig(haldane,k_vec)

#bandstructure
fig,ax = subplots()
ax.set_xlim(k_node[1],k_node[end])
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
for n in 1:length(k_node)
    ax.axvline(x=k_node[n],linewidth=0.5,color="k")
end
ax.set_title("Haldane model band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

#plot band
ax.plot(k_dist,eigvals')
fig.savefig("examples/haldane_band.pdf")

#Density of state
w,DOS = calc_DOS(haldane,nk)

#DOS plot
fig,ax = subplots()
ax.plot(w,DOS)
ax.set_ylim(0,maximum(DOS))
ax.set_title("Haldane model density of states")
ax.set_xlabel("Band energy")
ax.set_ylabel("Density of state")
fig.savefig("examples/haldane_dos.pdf")
