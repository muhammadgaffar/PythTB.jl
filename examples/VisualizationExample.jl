using PythTB

#lattice visualization example
#honeycomb lattice

#define lattice vectors for honeycomb lattice
lat = [[1,0],[0.5,sqrt(3)/2]]

#put orbital
orb = [[1/3,1/3],[2/3,2/3]]

#graphene model
dim_k = 2
dim_r = 2
model = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1.0
set_hop!(model,t,1,2,[0,0])
set_hop!(model,t,2,1,[1,0])
set_hop!(model,t,2,1,[0,1])

#visualize infinite model
fig,ax = visualize(model,1,2)
ax.set_title("Graphene, bulk")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.savefig("examples/visualize_bulk.pdf")

#finite model along direction x
num_atom = 8
finite_x = cut_piece(model,num_atom,1,glue_edgs=false)

#visualize finite model in x
fig,ax = visualize(finite_x,1,2)
ax.set_title("Graphene, ribbon")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.savefig("examples/visualize_ribbon.pdf")

# finite model along direction x and y
finite_xy = cut_piece(finite_x,num_atom,2,glue_edgs=false)

#visualize finite model in x-y
fig,ax= visualize(finite_xy,1,2)
ax.set_title("Graphene, finite")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.savefig("examples/visualize_finite.pdf")
