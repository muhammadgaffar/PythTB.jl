using PythTB

#one dimensional TB with finite chain

#we set function to get model where parameterized by λ
function set_model(t,Δ,λ)
    # one dimensional lattice
    lat = [[1.0]]
    # put three atom
    orb = [[0.0],[1/3],[2/3]]
    #the model
    model = tb_model(1,1,lat,orb)
    #set hopping
    set_hop!(model,t,1,2,[1])
    set_hop!(model,t,2,3,[1])
    set_hop!(model,t,3,1,[2])
    #set onsite
    onsite_0 = -Δ*cos(2π * (λ - 0))
    onsite_1 = -Δ*cos(2π * (λ - 1/3))
    onsite_2 = -Δ*cos(2π * (λ - 2/3))
    set_onsite!(model,[onsite_0,onsite_1,onsite_2])

    return model
end

#set parameters
Δ = 2.0
t = -1.3

#effect of λ
path_steps = 31
λs = LinRange(0,1,path_steps)

#init model
model1d = set_model(t,Δ,0)

#kmesh in BZ
nk = 51
k_vec,k_dist,k_node = k_path(model1d,"fullc",nk)

#init wf
setup_wf_array!(model1d,[path_steps,nk])

#fill wf
for iλ in 1:path_steps
    λ = λs[iλ]
    model = set_model(t,Δ,λ)
    #solve wavefunction
    eigval,eigvec = solve_eig(model,k_vec,eig_vec=true)
    for ik in 1:nk
        model1d.array[iλ,ik] = eigvec[:,ik,:]
    end
end

#compute berry flux
println("Chern numbers for rising fillings")
println("Band 1 = ",berry_flux(model1d,[1]) ./ 2π)
println("Band 1,2 = ",berry_flux(model1d,[1,2]) ./ 2π)
println("Band 1,2,3 = ",berry_flux(model1d,[1,2,3]) ./ 2π)
println("")
println("Chern numbers for individual band")
println("Band 1 = ",berry_flux(model1d,[1]) ./ 2π)
println("Band 2 = ",berry_flux(model1d,[2]) ./ 2π)
println("Band 3 = ",berry_flux(model1d,[3]) ./ 2π)
println("")

#loop parameter λ again now for finite chain
path_steps = 241
λs = LinRange(0,1,path_steps)

#length of chain
num_cells = 10
num_orb = 3*num_cells

#initialize array for eigenvalues and x expectation for 1d finite chain
ch_eval = zeros(num_orb,path_steps)
ch_xexp = zeros(num_orb,path_steps)

#fill array
for iλ in 1:path_steps
    λ = λs[iλ]

    #construct dummy model
    model = set_model(t,Δ,λ)
    #cut model for finite chain
    ch_model = cut_piece(model,num_cells,1)
    eigval, eigvec = solve_eig(ch_model,eig_vec=true)

    #save
    ch_eval[:,iλ] = eigval
    ch_xexp[:,iλ] = position_expectation(ch_model,eigvec,1)
end

#plot
using PyPlot
fig, ax = subplots()

for n in 1:num_orb
    # diminish the size of the ones on the borderline
    xcut = 2   # discard points below this
    xfull = 4  # use sybols of full size above this
    size = (ch_xexp[n,:] .- xcut) ./ (xfull - xcut)
    for i in 1:path_steps
        size[i] = minimum([size[i],1])
        size[i] = maximum([size[i],0.1])
    end
    ax.scatter(λs,ch_eval[n,:],edgecolors="none",s=size*6,c="k")
end

ax.set_title("Eigenenergies for finite chain of 3-site-model")
ax.set_xlabel(L"Parameter $\lambda$")
ax.set_ylabel("Energy")
ax.set_xlim(0,1)
fig.savefig("examples/3site_endstates.pdf")
