using Plots

import Base: show

export show, visualize_2d

function show(io::IO, model::model)
    println(io,"Tight-binding model summary")
    println(io,"===========================================================")
    println(io,"k-space dimension = ", model.dim_k)
    println(io,"r-space dimension = ", model.dim_r)
    println(io,"periodic direction = ", model.per)
    println(io,"number of orbitals = ", model.norb)
    println(io,"number of electronic states = ", model.nstates)
    println(io,"------------------------------------------")
    println(io,"lattice vectors:")
    for i in 1:size(model.latvec,1)
        o = model.latvec[i,:]
        println(io,"# â[$i] ==> ", round.(o,digits=4))
    end
    println(io,"------------------------------------------")
    println(io,"position of orbitals:")
    for i in 1:size(model.orb,1)
        o = model.orb[i,:]
        println(io,"# r[$i] ==> ", round.(o,digits=4))
    end
    println(io,"------------------------------------------")
    println(io,"site energies:")
    if model.nspin == 2
        for i in 1:size(model.site_energies,1)
            o = model.site_energies[:,:,i]
            println(io,"# ε[$i] : ")
            show(IOContext(io, :limit=>false), MIME"text/plain"(), o)
            println(io,"")
        end
    else
        for i in 1:size(model.site_energies,1)
            o = model.site_energies[i]
            println(io,"# ε[$i] ==> ", round.(o,digits=4))
        end
    end
    println(io,"------------------------------------------")
    println(io,"hoppings:")
    if model.hoppings == []
        println(io,"no hopping")
    elseif model.nspin == 2
        for i in 1:size(model.hoppings,1)
            o = model.hoppings[i,:]
            if length(o) == 4
                println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)+$(o[4])⟩ : ")
                show(IOContext(io, :limit=>false), MIME"text/plain"(), o[1])
                println(io,"")
                println(io,"")
            else
                println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)⟩ : ")
                show(IOContext(io, :limit=>false), MIME"text/plain"(), o[1])
                println(io,"")
                println(io,"")
            end
        end
    else
        for i in 1:size(model.hoppings,1)
            o = model.hoppings[i,:]
            if length(o) == 4
                println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)+$(o[4])⟩ = ", round.(complex(o[1]),digits=4))
            else
                println(io,"⟨$(o[2] + 1)|H|$(o[3] + 1)⟩ = ", round.(complex(o[1]),digits=4) )
            end
        end
    end

end

function show(io::IO,pauli::pauli)
    println(io,"Pauli matrices (x,y,z)")
end

function visualize_2d(model::model,dir_first,dir_second;
                    eig_dr=nothing,draw_hoppings=true)

    #check eig_dr format
    if eig_dr != nothing
        if length(eig_dr) != model.norb
            error("eig_dr must be array of size norb")
        end
    end

    #check if dir_second had to be specified
    if dir_second == nothing && model.dim_r > 1
        error("Need to specify index of second coordinate for projection")
    end

    function proj(v)
        coord_x = v[dir_first,:]
        if dir_second == nothing
            coord_y = 0.0
        else
            coord_y = v[dir_second,:]
        end
        return [coord_x coord_y]'
    end

    function cart(redc)
        return redc' * model.latvec
    end

    #initialize plot
    p = plot(legend=false)
    linestyle = :dash

    #draw latvec
    v = proj(model.latvec)
    for i in model.per
        x = [0,v[i,1]]
        y = [0,v[i,2]]
        plot!(x,y,color="red",ls=linestyle)

        if length(model.per) == 2
            if i == 1
                origin = 2
            elseif i == 2
                origin = 1
            end
            x .+= v[origin,1]
            y .+= v[origin,2]
            plot!(x,y,color="red",ls=linestyle)
        end
    end

    #draw orbital
    for i in 1:model.norb
        orb = vec(proj(cart(model.orb[i,:])'))
        x = [orb[1]]
        y = [orb[2]]
        plot!(x,y,marker=:circle,color="blue")
    end

    #draw hoppings
    if draw_hoppings == true
        for i in 1:size(model.hoppings,1)
            h = model.hoppings[i,:]
            for s in 1:2
                pos_i = copy(model.orb[h[2]+1,:])
                pos_j = copy(model.orb[h[3]+1,:])
                if model.dim_k != 0
                    if s == 1
                        pos_j[model.per] = pos_j[model.per] .+ h[4][model.per]
                    elseif s == 2
                        pos_i[model.per] = pos_i[model.per] .- h[4][model.per]
                    end
                end
                pos_i = vec(proj(cart(pos_i)'))
                pos_j = vec(proj(cart(pos_j)'))
                hop_lines = hcat(pos_i,pos_j)
                plot!(hop_lines[1,:],hop_lines[2,:],color="green")
            end
        end
    end

    #draw eigenstate
    if eig_dr != nothing
        for i in 1:model.norb
            pos = cart(model.orb[i,:])
            pos = vec(proj(pos'))
            nrm = real(eig_dr[i]' * eig_dr[i])
            nrm_rad = 2nrm*model.norb
            phase = angle(eig_dr[i])
            plot!([pos[1]],[pos[2]],color="grey",marker=:circle,
                    markersize=nrm_rad,markeralpha=0.5)
        end
    end

    return p
end
