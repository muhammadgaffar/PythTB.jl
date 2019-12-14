using PythTB

#Calculate berry phases along k_x as function of k_y
#same as calculate 1d wannier center along x

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
using PyPlot
fig, ax = subplots()
ky = LinRange(0,1,length(phi_a_1))

ax.plot(ky,phi_a_1,"ro")
ax.plot(ky,phi_b_1,"go")
ax.plot(ky,phi_c_1,"bo")
ax.set_title("Berry phase for lower (red), top (green), both bands (blue)")
ax.set_xlabel(L"$k_y$")
ax.set_ylabel(L"Berry phase along $k_x$")
ax.set_xlim(0.,1.)
ax.set_ylim(-7.,7.)
ax.yaxis.set_ticks([-2π,-π,0,π,2π])
ax.set_yticklabels((L"$-2\pi$",L"$-\pi$",L"$0$",L"$\pi$", L"$2\pi$"))
fig.savefig("examples/haldane_bp_phase.pdf")

#printout berry flux
println("Berry flux = ", flux_a_1)
