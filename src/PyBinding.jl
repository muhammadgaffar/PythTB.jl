function to_pybinding(tb)
    if tb.dim_r  == 1
        throw("current dimension for pybinding must 2 or 3")
    end

    latvec = tb.latvec * 0.1 # pybinding use nm, convert to nm.
    if tb.dim_r == 2
        model = pb.Lattice( a1 = latvec[1,:],
                            a2 = latvec[2,:] )
    elseif tb.dim_r == 3
        model = pb.Lattice( a1 = latvec[1,:],
                            a2 = latvec[2,:],
                            a3 = latvec[3,:] )
    end

    for i in 1:tb.norb
        pos = tb.orb[i,:]' * latvec
        onsite = tb.site_energies[i,:,:]
        model.add_one_sublattice(string(i),pos,onsite)
    end

    for i in 1:size(tb.hoppings,1)
        h = tb.hoppings[i,:]
        i = string(h[2] .+ 1)
        j = string(h[3] .+ 1)
        lat = h[4]
        hop = h[1]
        model.add_one_hopping(lat,i,j, hop)
    end

    return model
end
