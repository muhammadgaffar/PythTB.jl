using PythTB

lat = [1.0 0.0;
       0.5 sqrt(3)/2]
orb = [1/3 1/3;
       2/3 2/3]
model = tb_model(2,2,lat,orb,nspin=2)

# set parameters
Δ =  0.7
t = -1.0
soc = -0.24
rashba = 0.05
soc_2 = [-0.06, -0.24]

#set onsite
set_onsite!(model, [Δ,-Δ])

# soc effect
r3h = sqrt(3)/2
sigma_a = 0.5 * σ.x - r3h * σ.y
sigma_b = 0.5 * σ.x + r3h * σ.y
sigma_c = -1.0 * σ.x

# spin independent nearest neighbour hoppings
for latvec in ([0,0],[-1,0],[0,-1])
      set_hop!(model, t, 1,2, latvec)
end
# spin dependent next-nearest neighbour hoppings
for latvec in ([1,0],[-1,1],[0,-1])
      set_hop!(model, 1im*soc*σ.z, 1,1, latvec)
      set_hop!(model, 1im*soc*σ.z, 2,2, -latvec)
end
# spin-flip nearest neighbour hoppings
set_hop!(model, 1im*rashba*sigma_a, 1,2, [0,0], mode="add")
set_hop!(model, 1im*rashba*sigma_b, 1,2, [-1,0], mode="add")
set_hop!(model, 1im*rashba*sigma_c, 1,2, [0,-1], mode="add")

# set up ribbon model
width = 20
nk = 101
rib_model = cut_piece(model, width, 2, glue_edgs=false)
kvec,kdist,knode = k_path(rib_model,"full",nk)
rib_eval = solve_eig(rib_model, kvec)

# show plot
fig = plot(framestyle=:box,legend=false)
plot!(fig,kdist,rib_eval',c="black")
ylabel!(fig,"Edge band structure")
xlabel!(fig,"k / 2\\pi")
ylims!(fig,(-2.5,2.5))
xlims!(fig,(0,1))
