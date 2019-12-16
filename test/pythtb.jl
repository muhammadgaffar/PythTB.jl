using PythTB
using Test

function test_model_nspin1()
    lat = [1 0; 0.5 0.5]
    orb = [0 0; 0.5 2/3]
    model = tb_model(2,2,lat,orb)
    set_onsite!(model,[-1.2,0.9])
    set_hop!(model,1.0,1,2,[0,0])
    set_hop!(model,1.0,1,1,[1,0])
    set_hop!(model,1.0,1,1,[0,1])
    return model
end

function test_model_nspin2()
    lat = [1 0; 0.5 0.5]
    orb = [0 0; 0.5 2/3]
    model = tb_model(2,2,lat,orb,nspin=2)
    set_onsite!(model,[-1.2,0.9])
    set_hop!(model,1.0,1,2,[0,0])
    set_hop!(model,1.0,1,1,[1,0])
    set_hop!(model,1.0,1,1,[0,1])
    return model
end

#test functionality on generating model with spin = 1
@testset "nspin1" begin
    model = test_model_nspin1()
    @test model.nspin == 1
    @test length(model.site_energies[1]) == 1
    @test length(model.hoppings[1,1]) == 1
    @test size(hamiltonian(model,[rand(),rand()])) == (2,2)
end

#test functionality on generating model with spin = 2
@testset "nspin2" begin
    model = test_model_nspin2()
    @test model.nspin == 2
    @test size(model.site_energies) == (2,2,2)
    @test size(model.hoppings[1]) == (2,2)
    @test size(hamiltonian(model,[rand(),rand()])) == (2,2,2,2)
    @test size(val_to_block(model,2*σ.x)) == (2,2)
end
