using Plots

export visualize_2d

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
