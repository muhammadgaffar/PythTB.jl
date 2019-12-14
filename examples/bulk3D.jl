using PythTB
using PyPlot

# example of 3d tight binding
# BCC crystal

#define BCC lattice
lat = [[1,0,0],[0,1,0],[0,0,1]]

#put atom
orb = [[0,0,0],[0.5,0.5,0.5]]

#get the model
bcc = tb_model(3,3,lat,orb)

#set onsite
Δ = 0.5
set_onsite!(bcc,[-Δ,Δ])
#set hop
t = -1.0
set_hop!(bcc,t,1,2,[0,0,0])
set_hop!(bcc,t,1,2,[0,-1,0])
set_hop!(bcc,t,1,2,[1,0,0])
set_hop!(bcc,t,1,2,[1,-1,0])
set_hop!(bcc,t,1,2,[0,0,-1])
set_hop!(bcc,t,1,2,[0,-1,-1])
set_hop!(bcc,t,1,2,[1,0,-1])
set_hop!(bcc,t,1,2,[1,-1,-1])

#generate kpoint
path = [[0,0,0],[-1/2,1/2,1/2],[0,1/2,0],[0,0,0],[1/4,1/4,1/4],[0,1/2,0],[1/4,1/4,1/4],[-1/2,1/2,1/2]]
labels = [L"$\Gamma$","H","N",L"\Gamma","P","N","P","H"]
k_vec,k_dist,k_node = k_path(bcc,path,101)

#solve band
eigvals = solve_eig(bcc,k_vec)

#plot
fig,ax = subplots(figsize=(6,4))
ax.plot(k_dist,eigvals',"k-")

ax.set_title("BandStructure for body-centered cubic")
ax.set_ylabel("Band energy")
ax.set_xlabel("Path in k-space")

ax.set_xlim(k_node[1],k_node[end])
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
for n in 1:length(k_node)
    ax.axvline(x=k_node[n],linewidth=0.5,color="k")
end
fig.savefig("examples/bcc_band.pdf")

#get density of state
dos = calc_DOS(bcc,80,150,window=[-10,10])

#plot
fig,ax = subplots(figsize=(6,4))
ax.plot(dos[1],dos[2])
ax.set_title("DOS of BCC")
ax.set_xlabel("Energy")
ax.set_ylabel("DOS(ω)")
ax.set_ylim(0,maximum(dos[2]))
fig.savefig("examples/bcc_dos.pdf")
