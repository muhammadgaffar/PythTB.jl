using Shell
using DelimitedFiles

struct command
    x
    tools
    py
end

struct config
    pb
    conf
    calc
    kite
end

function KITE(PATH::String)
    pushfirst!(PyVector(pyimport("sys")."path"),PATH)

    x = PATH * "/build/KITEx "
    tools = PATH * "/tools/build/KITE-tools "

    cmd = "sed '788,794d' " * PATH * "/kite.py > " * PATH * "/_kite.py"
    Shell.run(cmd)

    pushfirst!(PyVector(pyimport("sys")."path"),PATH)
    py = pyimport("_kite")

    return command(x,tools,py)
end

function configuration(tb::model,kite::command; parallel=[4,4],lat_size=[1024,1024],
                    PBC=[true,true],is_complex=false,precision=1,w_range=nothing)

    if tb.dim_r == 1
        model = pb.Lattice(a1 = [tb.latvec[1],0], a2 = [0,1])
        parallel = [4,1]
        lat_size = [4*1024,64]
    elseif tb.dim_r == 2
        model = pb.Lattice(a1 = tb.latvec[1,:], a2 = tb.latvec[2,:])
    elseif tb.dim_r == 3
        model = pb.Lattice( a1 = tb.latvec[1,:],
                            a2 = tb.latvec[2,:],
                            a3 = tb.latvec[3,:] )
    end
    for i in 1:tb.norb
        if tb.dim_r == 1
            model.add_one_sublattice(string(i),[tb.orb[i],0],tb.site_energies[i])
        else
            if tb.nspin==1
                model.add_one_sublattice(string(i),tb.orb[i,:],tb.site_energies[i])
            elseif tb.nspin==2
                model.add_one_sublattice(string(i),tb.orb[i,:],tb.site_energies[i,:,:])
            end
        end
    end
    for i in 1:size(tb.hoppings,1)
        h = tb.hoppings[i,:]
        if tb.dim_r == 1
            model.add_one_hopping([h[4],0],string(h[2] .+ 1),string(h[3] .+ 1),h[1])
        else
            model.add_one_hopping(h[4],string(h[2] .+ 1),string(h[3] .+ 1),h[1])
        end
    end

    complex = is_complex
    prec = precision

    conf = kite.py.Configuration(divisions=parallel,length=lat_size,boundaries=PBC,
                            is_complex=complex,precision=prec,spectrum_range=w_range)
    calc = kite.py.Calculation(conf)

    return config(model,conf,calc,kite)
end

function calc_dos(config::config; nw, num_random=1, num_moments=512)
    nrand = num_random
    nmomt = num_moments
    prefix = "tmp.h5"

    config.calc.dos(num_points=nw,num_random=nrand,num_moments=nmomt)
    config.kite.py.config_system(config.pb,config.conf,config.calc,filename=prefix)

    cmd = config.kite.x * prefix
    Shell.run(cmd)
    cmd = config.kite.tools * prefix
    Shell.run(cmd)

    x = readdlm("dos.dat")

    rm = "rm -rf dos.dat"
    Shell.run(rm)
    rm = "rm -rf tmp.h5"
    Shell.run(rm)

    return x[:,1], x[:,2]
end

function calc_opticalConductivity(config::config; nhw,num_random=1,num_moments=512,T,direction)
    nrand = num_random
    nmomt = num_moments
    dir   = direction
    prefix = "tmp.h5"

    config.calc.conductivity_optical(num_points=nhw,num_random=nrand,
                                num_moments=nmomt,direction=dir,temperature=T)
    config.kite.py.config_system(config.pb,config.conf,config.calc,filename=prefix)

    cmd = config.kite.x * prefix
    Shell.run(cmd)
    cmd = config.kite.tools * prefix
    Shell.run(cmd)

    x = readdlm("optcond.dat")

    rm = "rm -rf optcond.dat"
    Shell.run(rm)
    rm = "rm -rf tmp.h5"
    Shell.run(rm)

    return x[:,1], x[:,2] .+ 1im .* x[:,3]
end

function calc_dcConductivity(config::config; nw,num_random=1,num_moments=512,T,direction)
    nrand = num_random
    nmomt = num_moments
    dir   = direction
    prefix = "tmp.h5"

    config.calc.conductivity_dc(num_points=nw, num_moments=nmomt, num_random=nrand,
                                direction=dir, temperature=T)
    config.kite.py.config_system(config.pb,config.conf,config.calc,filename=prefix)

    println("done")

    cmd = config.kite.x * prefix
    Shell.run(cmd)
    cmd = config.kite.tools * prefix
    Shell.run(cmd)

    x = readdlm("condDC.dat")

    rm = "rm -rf condDC.dat"
    Shell.run(rm)
    rm = "rm -rf tmp.h5"
    Shell.run(rm)

    return x[:,1], x[:,2], x[:,3]
end
