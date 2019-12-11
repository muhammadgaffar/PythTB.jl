using PythTB
using PyPlot

#define one-dimensional lattice vectors
lat = [[1.0]]

#put atom/orbital position
orb = [[0.0]]

#make new model
dim_k = 1 #dimension of k-space (this shows the periodicity)
dim_r = 1 #dimension of real-space
model1D = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1.0 #hopping amplitude
set_hop!(model1D,t,1,1,[1]) #orbital 1 to another orbital 1

#show summary of the model
show(model1D)

#generate k-points
nk = 201
k_vec,k_dist,k_node = k_path(model1D,"full",nk)

#evaluate energy over k
eigvals = solve_eig(model1D,k_vec)

#plot and save the figure
fig, ax = subplots(figsize=(6,4))
ax.plot(k_dist,eigvals')

labels = [L"$0$", L"$\pi$", L"$2\pi$"]
ax.set_title("Band Structure of One Dimensional Chain")
ax.set_ylabel("Band Energy (eV)")
ax.set_xlabel("Path in K-space")
ax.set_xlim(k_node[1],k_node[end])
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
for k in k_node
    ax.axvline(x = k, linewidth = 0.5, color = "k")
end
fig.savefig("examples/Band1D.pdf")
