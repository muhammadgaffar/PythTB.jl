using PythTB
using Plots

# This example below will read the tight-binding model from the Wannier90
# calculation, create a simplified model in which some small hopping terms
# are ignored, and finally plot the interpolated band structure.

# read output from wannier90 from folder
path = "examples/Wannier90_plus_pythTB/wannier90_example"
prefix = "silicon"
#fermi level
fermi_ev = 6.2285135
#generate TB model
silicon = w90_model(path,prefix,min_hopping_norm=0.01,zero_energy=fermi_ev)
#show the model
show(silicon)

# path on silicon
path = [0.5 0.5 0.5; 0 0 0; 0.5 -0.5 0.0; 0.375 -0.375 0.0; 0 0 0]
labels = ["L","\\Gamma", "X", "K", "\\Gamma"]
k_vec,k_dist,k_node = k_path(silicon,path,190)

#solve model
eigvals = solve_eig(silicon,k_vec)

#get w90 band file
w90_kpt, w90_evals = bands_consistency(silicon)

#plot and compare to w90 band
fig = plot(framestyle=:box,legend=false)
plot!(fig,k_dist,eigvals',c="black")
plot!(fig,k_dist,w90_evals' .- fermi_ev,c="red")
vline!(fig,k_node,lw=0.5,c="black",alpha=0.5)
xlabel!(fig,"Path in k-space")
ylabel!(fig,"Band energy (eV)")
title!(fig,"Quick example of Wannier90 + PythTB")
xlims!(fig,(k_node[1],k_node[end]))
xticks!(k_node,labels)
# now instead, we use kpoints from w90_kpt than manual generated path
# because the difference lie in here

#solve model
eigvals = solve_eig(silicon,w90_kpt)

#plot and compare to w90 band
fig2 = plot(framestyle=:box,legend=false)
plot!(fig2,k_dist,eigvals',c="black")
plot!(fig2,k_dist,w90_evals' .- fermi_ev,c="red")
xlabel!(fig2,"Path in k-space")
ylabel!(fig2,"Band energy (eV)")
title!(fig2,"Wannier90 + PythTB with kpoints consistency")
xlims!(fig2,(k_node[1],k_node[end]))
xticks!(fig2,k_node,labels)
