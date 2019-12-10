using PythTB

# This example considers a simple three-site one-dimensional tight-binding model
# parametrized by some parameter λ. As λ is changed from 0 to 1,
# the deepest onsite term is moved from the first to second, then to the third,
# and then back to the first tight-binding orbital. Therefore, we expect that
# Wannier center of the lowest band will shift by one lattice vector
# as λ changes from 0 to 1.

#one dimensional lattice
lat = [[1.0]]

#put three orbital
orb = [[0.0],[1/3],[2/3]]

#generate model
dim_k,dim_r = (1,1)
model1d = tb_model(dim_k,dim_r,lat,orb)

#set onsite
Δ = 2.0

#set hop
t = -1.0
set_hop!(model1d,t,1,2,[1])
set_hop!(model1d,t,2,3,[1])
set_hop!(model1d,t,3,1,[2])

#initialize figure for onsite and band
using PyPlot
fig_onsite, ax_onsite = subplots(figsize=(6,4))
fig_band, ax_band = subplots(figsize=(6,4))

#perform the effect on change of onsite term
path_steps = 31
#parametrized by λ
λs = LinRange(0,1,path_steps)
#number of kpoints
nk = 31
#setup wf_array, save for both kpoint and λ
setup_wf_array!(model1d,[nk,path_steps])

for i in 1:path_steps
    λ = λs[i]
    onsite_0 = -Δ*cos(2π* (λ - 0))
    onsite_1 = -Δ*cos(2π* (λ - 1/3))
    onsite_2 = -Δ*cos(2π* (λ - 2/3))

    #update onsite
    set_onsite!(model1d,[onsite_0,onsite_1,onsite_2],mode="reset")

    #kmesh
    k_vec,k_dist,k_node = k_path(model1d,"fullc",nk)

    #wavefunction
    eigval, eigvec = solve_eig(model1d,k_vec,eig_vec=true)
    for j in 1:nk
        model1d.array[j,i] = eigvec[:,j,:]
    end

    #plot onsite
    ax_onsite.scatter(λ,onsite_0,c="r")
    ax_onsite.scatter(λ,onsite_1,c="g")
    ax_onsite.scatter(λ,onsite_2,c="b")
    #plot band
    for k in 1:size(eigval,1)
        ax_band.plot(k_dist,eigval[k,:],"k-",linewidth=2.0)
    end
end

# impose periodic boundary condition along k-space direction only
# (so that |ψ(k)⟩ at k=0 and k=1 have the same phase)
impose_pbc!(model1d,1,1)

#compute berry phase
phase = berry_phase(model1d,[1],1)

#wannier center phase
wann_center = phase ./ 2π

#plot wannier
fig_wann, ax_wann = subplots(figsize=(6,4))
ax_wann.plot(λs,wann_center,"ko-")

#compute berry flux (integrated curvature)
final = berry_flux(model1d,[1])
println("Berry flux in k-λ space = ", final)

#onsite figure
ax_onsite.set_title("Onsite energy for all three orbitals")
ax_onsite.set_xlabel("λ")
ax_onsite.set_ylabel("Onsite terms")
ax_onsite.set_xlim(0,1)
fig_onsite.savefig("examples/3site_onsite.pdf")
#band figure
ax_band.set_title("Band structure")
ax_band.set_xlabel("Path in k-vector")
ax_band.set_ylabel("Band energies")
ax_band.set_xlim(0,1)
fig_band.savefig("examples/3site_band.pdf")
#wannier center figure
ax_wann.set_title("Center of Wannier function")
ax_wann.set_xlabel("Lambda parameter")
ax_wann.set_ylabel("Center (reduced coordinate)")
ax_wann.set_xlim(0,1)
fig_wann.savefig("examples/3site_wann.pdf")
