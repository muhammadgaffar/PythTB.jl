using PythTB
using PyPlot

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
path = [[0.5,0.5,0.5],[0,0,0],[0.5,-0.5,0.0],[0.375,-0.375,0.0],[0,0,0]]
labels = [L"$L$", L"$\Gamma$", L"$X$", L"$K$", L"$\Gamma$"]
k_vec,k_dist,k_node = k_path(silicon,path,190)

#solve model
eigvals = solve_eig(silicon,k_vec)

#get w90 band file
w90_kpt, w90_evals = bands_consistency(silicon)

#plot and compare to w90 band
fig,ax = subplots()
ax.plot(k_dist,eigvals',"k-")
ax.plot(k_dist,w90_evals' .- fermi_ev,"r")

for n in 1:length(k_node)
    ax.axvline(x=k_node[n],linewidth=0.5,color="k")
end
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy (eV)")
ax.set_title("Quick example of Wannier90 + PythTB")
ax.set_xlim(k_node[1],k_node[end])
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
fig.savefig("examples/Wannier90_plus_pythTB/silicon_quick.pdf")

# now instead, we use kpoints from w90_kpt than manual generated path
# because the difference lie in here

#solve model
eigvals = solve_eig(silicon,w90_kpt)

#plot and compare to w90 band
fig,ax = subplots()
ax.plot(k_dist,eigvals',"k-")
ax.plot(k_dist,w90_evals' .- fermi_ev,"r")

for n in 1:length(k_node)
    ax.axvline(x=k_node[n],linewidth=0.5,color="k")
end
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy (eV)")
ax.set_title("Wannier90 + PythTB with kpoints consistency")
ax.set_xlim(k_node[1],k_node[end])
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
fig.savefig("examples/Wannier90_plus_pythTB/silicon_consistency.pdf")
