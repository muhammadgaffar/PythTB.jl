using Shell
using DelimitedFiles
using Suppressor

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

function cleanup()
    run(`rm -rf tmp.h5`)
    run(`rm -rf optcond.dat`)
    run(`rm -rf condDC.dat`)
    nothing
end

function configuration(tb::model,kite::command; parallel,lat_size,
                    PBC,complex=false,precision=1,w_range=nothing)

    model = to_pybinding(tb)

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

    println("Calculating Density of State...")

    @suppress config.calc.dos(num_points=nw,num_random=nrand,num_moments=nmomt)
    @suppress config.kite.py.config_system(config.pb,config.conf,config.calc,filename=prefix)

    cmd = config.kite.x * prefix
    @suppress Shell.run(cmd)
    cmd = config.kite.tools * prefix
    @suppress Shell.run(cmd)

    println("Done...!")

    x = readdlm("dos.dat")

    cleanup()

    return x[:,1], x[:,2]
end

function calc_opticalConductivity(config::config; nhw,num_random=1,num_moments=512,T,dir)
    nrand = num_random
    nmomt = num_moments
    prefix = "tmp.h5"

    println("Calculating Optical Conductivity...")

    @suppress config.calc.conductivity_optical(num_points=nhw,num_random=nrand,
                                num_moments=nmomt,direction=dir,temperature=T)
    @suppress config.kite.py.config_system(config.pb,config.conf,config.calc,filename=prefix)

    cmd = config.kite.x * prefix
    @suppress Shell.run(cmd)
    cmd = config.kite.tools * prefix
    @suppress Shell.run(cmd)

    println("Done...!")

    x = readdlm("optcond.dat")

    cleanup()

    return x[:,1], x[:,2] .+ 1im .* x[:,3]
end

function calc_dcConductivity(config::config; nw,num_random=1,num_moments=512,T,dir)
    nrand = num_random
    nmomt = num_moments
    prefix = "tmp.h5"

    println("Calculating DC Conductivity...")

    @suppress config.calc.conductivity_dc(num_points=nw, num_moments=nmomt, num_random=nrand,
                                direction=dir, temperature=T)
    @suppress config.kite.py.config_system(config.pb,config.conf,config.calc,filename=prefix)

    cmd = config.kite.x * prefix
    @suppress Shell.run(cmd)
    cmd = config.kite.tools * prefix
    @suppress Shell.run(cmd)

    x = readdlm("condDC.dat")

    println("Done...!")

    cleanup()

    return x[:,1], x[:,2], x[:,3]
end
