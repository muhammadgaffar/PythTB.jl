import Base: show

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
