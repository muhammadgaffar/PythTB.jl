using PyPlot
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
    for (i,o) in enumerate(model.latvec)
        println(io,"# â[$i] ==> ", round.(o,digits=4))
    end
    println(io,"------------------------------------------")
    println(io,"position of orbitals:")
    for (i,o) in enumerate(model.orb)
        println(io,"# r[$i] ==> ", round.(o,digits=4))
    end
    println(io,"------------------------------------------")
    println(io,"site energies:")
    for (i,o) in enumerate(model.site_energies)
        println(io,"# ε[$i] ==> ", round(0,digits=4))
    end
    println(io,"------------------------------------------")
    println(io,"hoppings:")
    if model.hoppings == 0 || model.hoppings == nothing
        println(io,"no hopping")
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

function visualize(model::model,dir_first,dir_second=nothing;eig_dr=nothing,draw_hoppings=true,ph_color="black")
    return model.model.visualize(dir_first - 1,dir_second - 1,eig_dr,draw_hoppings,ph_color)
end
