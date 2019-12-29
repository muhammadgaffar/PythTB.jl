using PythTB

# This is example how to use KITE for DOS and transport calculation in graphene

#my kite path installation
PATH = "/home/mgaffar/envs/kite"
kite = KITE(PATH)

#honeycomb lattice
lat = [1 0; 0.5 sqrt(3)/2]
#put atom position for graphene
orb = [1/3 1/3; 2/3 2/3]

#tb model
tb = tb_model(2,2,lat,orb)

#set hop
t = -1.0
set_hop!(tb,t,1,2,[0,0])
set_hop!(tb,t,2,1,[1,0])
set_hop!(tb,t,2,1,[0,1])

#interfacing and KITE configuration for graphene
conf = configuration(tb,kite)

#calculate the transport properties
w,dos = calc_dos(conf; nw=1000)
hw,sxx = calc_opticalConductivity(conf,nhw=1000,T=300,direction="xx")
ww,cxx = calc_dcConductivity(conf, nw=1000,T=300,direction="xx")

using Plots

#plot DOS
dosplot = plot(w,dos)
dosplot

#plot dc conductivity
dcplot = plot(ww,cxx)
dcplot

#plot opt conductivity
oplot = plot(hw,real(sxx))
plot!(oplot,hw,imag(sxx))
ylims!(oplot,(-1,3))
xlims!(oplot,(0,0.2))
