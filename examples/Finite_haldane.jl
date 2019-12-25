using PythTB
using Plots

#finite haldane model with periodic and non-periodic boundary condition

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

#cut model to make it as finite model along direction x and y
#with no periodic boundary
num_x = 20
num_y = 20
finite_haldane = cut_piece(haldane,num_x,1,glue_edgs=false)
finite_haldane = cut_piece(finite_haldane,num_y,2,glue_edgs=false)

#solve energy
eigvals_finite = solve_eig(finite_haldane)

#plot DOS
w,DOS = calc_DOS(finite_haldane,finite=true)

fig = plot(w,DOS)
ylims!(fig,(0,maximum(DOS)))

title!(fig,"Finite Haldane model without PBC")
xlabel!(fig,"Energy")
ylabel!(fig,"Density of state")
display(fig)

#visualize the finite model
fig2 = visualize_2d(finite_haldane,1,2)
title!(fig2,"Finite model of haldane (20 x 20)")
ylabel!(fig,"Y")
xlabel!(fig,"X")

#cut model to make it as finite model along direction x and y
#with periodic boundary
num_x = 20
num_y = 20
finite_haldane = cut_piece(haldane,num_x,1,glue_edgs=true)
finite_haldane = cut_piece(finite_haldane,num_y,2,glue_edgs=true)

#solve energy
eigvals_finite = solve_eig(finite_haldane)

#plot DOS
w,DOS = calc_DOS(finite_haldane,finite=true)

fig = plot(w,DOS)
ylims!(fig,(0,maximum(DOS)))

title!(fig,"Finite Haldane model with PBC")
xlabel!(fig,"Energy")
ylabel!(fig,"Density of state")
display(fig)
