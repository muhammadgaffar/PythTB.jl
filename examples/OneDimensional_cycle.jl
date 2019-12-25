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
set_hop!(model1d,t,1,2,[0])
set_hop!(model1d,t,2,3,[0])
set_hop!(model1d,t,3,1,[1])

#initialize figure for onsite and band
using Plots
fig_onsite = plot(legend=false)
fig_band = plot(framstyle=:box,legend=false)

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
    scatter!(fig_onsite,[λ],[onsite_0],c="red")
    scatter!(fig_onsite,[λ],[onsite_1],c="green")
    scatter!(fig_onsite,[λ],[onsite_2],c="blue")
    #plot band
    for k in size(eigval,1)
        plot!(fig_band,k_dist,eigval[k,:],c="black",lw=2)
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
fig_wann = plot(legend=false)
plot!(fig_wann,λs,wann_center,marker=:circle)

#compute berry flux (integrated curvature)
final = berry_flux(model1d,[1])
println("Berry flux in k-λ space = ", final)

#onsite figure
title!(fig_onsite,"Onsite energy for all three orbitals")
xlabel!(fig_onsite,"\\lambda")
ylabel!(fig_onsite,"Onsite terms")
xlims!(fig_onsite,(0,1))

#band figure
title!(fig_band,"Band structure")
xlabel!(fig_band,"Path in k-vector")
ylabel!(fig_band,"Band energies")
xlims!(fig_band,(0,1))

#wannier center figure
title!(fig_wann,"Center of Wannier function")
xlabel!(fig_wann,"Lambda parameter")
ylabel!(fig_wann,"Center (reduced coordinate)")
xlims!(fig_wann,(0,1))
