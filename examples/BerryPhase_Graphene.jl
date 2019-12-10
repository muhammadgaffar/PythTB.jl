using PythTB

#Compute berry phase around dirac cone in graphene

#hexagonal lattice vectors
lat = [[1.0,0.0],[0.5,sqrt(3.0)/2.0]]

#put orbital
orb = [[1/3,1/3],[2/3,2/3]]

#generate the model
dim_k = 2
dim_r = 2
gr = tb_model(dim_k,dim_r,lat,orb)

#set onsite
Δ = -0.1
set_onsite!(gr,[-Δ,Δ])

#set hop
t = -1.0
set_hop!(gr,t,1,2,[1,1])
set_hop!(gr,t,2,1,[2,1])
set_hop!(gr,t,2,1,[1,2])

#summary of the models
show(gr)

#path aroud dirac cone
#circular path
circ_step = 51
circ_center = [1/3,2/3]
circ_radius = 0.05

#setup wave function array on the path
setup_wf_array!(gr,[circ_step])
#fill wf array
for i in 1:circ_step
    #kpoint of the path
    ang = 2π * (i-1) / (circ_step - 1)
    kpt = circ_radius .* [cos(ang),sin(ang)]
    #shift to center
    kpt += circ_center
    #calculate eigenvector (wavefunction)
    eigval, eigvec = solve_eig(gr,kpt,eig_vec=true)
    gr.array[i] = eigvec
end
#make sure the last point and first is same
gr.array[51] = gr.array[1]

#Compute berry phase along circular path
println("Berry phase along circle with radius = ", circ_radius)
println("Centered at k point = ", round.(circ_center,digits=4))
println("for band 1 equals =", berry_phase(gr,[1],1))
println("for band 2 equals =", berry_phase(gr,[2],1))
println("for both bands equals =", berry_phase(gr,[1,2],1))

#squre path
sq_step = 51
sq_center = [1/3,2/3]
sq_length = 0.1

#setup wave function array on the path
setup_wf_array!(gr,[sq_step,sq_step])
#save kpoint of the path for plot
kpts = zeros(sq_step,sq_step,2)
#fill wf array
for j in 1:sq_step, i in 1:sq_step
    kpt = sq_length .* [-0.5 + (i-1) / (sq_step - 1),
                        -0.5 + (j-1) / (sq_step - 1)]
    kpt += sq_center
    kpts[i,j,:] = kpt
    #calculate eigenvector (wavefunction)
    eigval, eigvec = solve_eig(gr,kpt,eig_vec=true)
    gr.array[i,j] = eigvec
end

#Compute berry flux along square path
println("Berry phase along circle with radius = ", sq_length)
println("Centered at k point = ", round.(sq_center,digits=4))
println("for band 1 equals =", berry_flux(gr,[1]))
println("for band 2 equals =", berry_flux(gr,[2]))
println("for both bands equals =", berry_flux(gr,[1,2]))

#plot berry phase
using PyPlot
plaq = berry_flux(gr,[1],individual_phases=true)

fig, ax = subplots(figsize=(6,4))
ax.imshow(plaq',origin="lower",
        extent=(kpts[1,1,1],kpts[50,1,1,],
                kpts[1,1,2],kpts[1,50,1]))
ax.set_title("Berry curvature near Dirac Cone")
ax.set_xlabel(L"$k_x$")
ax.set_ylabel(L"$k_y$")
fig.savefig("examples/cone_phases.pdf")
