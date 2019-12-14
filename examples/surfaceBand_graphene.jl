using PythTB
using PyPlot

# surface graphene using supercell

#define lattice vectors for honeycomb lattice
lat = [[1,0],[0.5,sqrt(3)/2]]

#put orbital
orb = [[1/3,1/3],[2/3,2/3]]

#graphene model
dim_k = 2
dim_r = 2
gr = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1
set_hop!(gr,t,1,2,[0,0])
set_hop!(gr,t,2,1,[0,1])
set_hop!(gr,t,2,1,[1,0])

#supercell graphene
sc_gr = make_supercell(gr,[[2,1],[-1,2]],to_home=true)

#cut model, finite along y direction
slab_gr = cut_piece(sc_gr,6,2,glue_edgs=false)

#visualize the finite model
fig,ax = visualize(slab_gr,1,2)
ax.set_title("Graphene, arbitrary surface")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.savefig("examples/supercell_gr.pdf")

#compute band structure
nk = 300
k_vec,k_dist,k_node = k_path(slab_gr,"full",nk)
eigvals = solve_eig(slab_gr,k_vec)

#plot
fig,ax = subplots(figsize=(6,4))
ax.plot(k_dist,eigvals',"k-")
# zoom in close to the zero energy
ax.set_xlim(k_dist[1],k_dist[end])
ax.set_ylim(-1.0,1.0)
# put title on top
ax.set_title("Graphene arbitrary surface band structure")
ax.set_xlabel("k parallel to edge")
ax.set_ylabel("Band energy")
ax.xaxis.set_ticks(k_node)
ax.set_xticklabels((L"$0$",L"$\pi$",L"$2\pi$"))
# make an PDF figure of a plot
fig.savefig("examples/supercell_band.pdf")
