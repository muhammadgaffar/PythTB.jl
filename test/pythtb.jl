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

function test_model_kmesh()
    lat = [0.8 0.2 0; 0.5 3/2 1; 0 0 2]
    orb = [0 0 0; 1 1 1]
    model = tb_model(3,3,lat,orb)
    set_hop!(model,1,1,1,[1,0,0])
    set_hop!(model,1,1,1,[0,1,1])
    set_hop!(model,0.8,1,2,[0,0,0])
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

@testset "kmesh and pauli" begin
    model = test_model_kmesh()
    nk = 5
    @test size(k_uniform_mesh(model,[nk,nk,nk])) == (nk^3,3)
    path = [1 2 3; 2/3 3 1/2; 0 0 0]
    nk = 24
    k_vec,k_dist,k_node = k_path(model,path,nk)
    @test size(k_vec) == (nk,model.dim_k)
    @test length(k_dist) == nk
    @test length(k_node) == size(path,2)
    @test typeof(σ) == PythTB.pauli
end
