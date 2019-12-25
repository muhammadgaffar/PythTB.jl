using PythTB
using Plots

#Calculate berry phases along k_x as function of k_y
#same as calculate 1d wannier center along x

#hexagonal lattice
lat = [1 0; 0.5 sqrt(3)/2]
#put orbital
orb = [1/3 1/3; 2/3 2/3]

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

#setup wf array that will be used for berry curvature
nk = 31
setup_wf_array!(haldane,[31,31])
#solve the wf array on grid, BZ start at -1/2 -1/2 point
solve_on_grid!(haldane,[-0.5,-0.5])

#calculate berry phases along kx direction
phi_a_1 = berry_phase(haldane,[1],1,contin=true)
phi_b_1 = berry_phase(haldane,[2],1,contin=true)
phi_c_1 = berry_phase(haldane,[1,2],1,contin=true)

#berry flux for first band
flux_a_1 = berry_flux(haldane,[0])

#plot berry phases
ky = LinRange(0,1,length(phi_a_1))

fig = plot()
plot!(fig,ky,phi_a_1,marker=:circle)
plot!(fig,ky,phi_b_1,marker=:circle)
plot!(fig,ky,phi_c_1,marker=:circle)
title!(fig,"Berry phase for lower (red), top (green), both bands (blue)")
xlabel!(fig,"ky")
ylabel!(fig,"Berry phase along kx")
xlims!(fig,0,1)
ylims!(fig,-7,7)
yticks!([-2π,-π,0,π,2π],["-2\\pi","-\\pi","0","\\pi","2\\pi"])

#printout berry flux
println("Berry flux = ", flux_a_1)
