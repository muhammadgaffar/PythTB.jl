using PyPlot
import Base: show

function show(io::IO, model::model)
    println("Tight-binding model")
    println("-------------------------")
    println("k-space dimension = ", model.dim_k)
    println("r-space dimension = ", model.dim_r)
    println("periodic direction = ", model.per)
    println("number of orbitals = ", model.norb)
    println("number of electronic states = ", model.nstates)
    println("lattice vectors:")
    for (i,o) in enumerate(model.latvec)
        println("# â[$i] ==> ", round.(o,digits=4))
    end
    println("position of orbitals:")
    for (i,o) in enumerate(model.orb)
        println("# r[$i] ==> ", round.(o,digits=4))
    end
    println("site energies:")
    for (i,o) in enumerate(model.site_energies)
        println("# ε[$i] ==> ", round(0,digits=4))
    end
    println("hoppings:")
    for i in 1:size(model.hoppings,1)
        o = model.hoppings[i,:]
        if length(o) == 4
            println("⟨$(o[2])|H|$(o[3])+$(o[4])⟩ = ", round.(complex(o[1]),digits=4))
        else
            println("⟨$(o[2])|H|$(o[3])⟩ = ", round.(complex(o[1]),digits=4) )
        end
    end
end

function visualize(model::model,dir_first,dir_second=nothing;eig_dr=nothing,draw_hoppings=true,ph_color="black")
    return model.model.visualize(dir_first,dir_second,eig_dr,draw_hoppings,ph_color)
end
