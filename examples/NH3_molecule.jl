using PythTB
using PyPlot

#Define lattice vectors in 3d space
lat = [[1,0,0],[0,1,0],[0,0,1]]

#put atoms
sq32 = sqrt(3) / 2
orb = [[ 2/3 * sq32,  0,   0],
       [-1/3 * sq32,  1/2, 0],
       [-1/3 * sq32, -1/2, 0],
       [ 0         ,  0,   1]]

#NH3 model
dim_k = 0 #because of single molecule
dim_r = 3
NH3 = tb_model(dim_k,dim_r,lat,orb)

#set onsite
Δ = 0.5
set_onsite!(NH3,[-Δ,-Δ,-Δ,Δ])

#set hopping
t = 1.0
set_hop!(NH3,t,1,2)
set_hop!(NH3,t,1,3)
set_hop!(NH3,t,1,4)
set_hop!(NH3,t,2,3)
set_hop!(NH3,t,2,4)
set_hop!(NH3,t,3,4)

#summary of the models
show(NH3)

#solve orbital energy
eigvals = solve_eig(NH3)

#plot orbital energy
fig, ax = subplots(figsize=(6,4))
ax.plot(eigvals,"o",markersize = 15)

ax.set_title("Orbital Energy for NH3 Molecule")
ax.set_xlabel("Orbital")
ax.set_ylabel("Orbital Energy")
ax.set_xticks([0,1,2,3])
fig.savefig("examples/NH3_energy.pdf")
