# PythTB.jl

PythTB.jl is Julia Wrapper for PythTB 1.7 Package in Pytho via PyCall.jl

## Installation
Make sure your current python environment is Anaconda, then install
PythTB python package first

```
pip install pythtb --upgrade
```
then open your julia REPL, add PythTB.jl package
```julia
] add https://github.com/muhammadgaffar/PythTB.jl
```

## Features
PythTB.jl can solve tight binding model for
- 0D cluster
- 1D chains and ladder
- 2D layers
- 3D crystal
- cluster, ribbon, slabs (finite model)

PythTB.jl include capabilities for
- computing electron eigenvalues and eigenvectors at selected k-points or on a mesh of k-points
- generating band-structure plots
- generating density-of-states plots
- calculate Berry Phase, and Berry curvature

PythTB.jl also provides an interface to the Wannier90 code, which can be used to take the output of a first-principles density-functional calculation and construct from it a tight-binding model, in the basis of Wannier functions, that accurately reproduces the first-principles bandstructure.

## Usage
This is example usage in honeycomb lattice

```julia
using PythTB

#define lattice vectors for honeycomb lattice
lat = [ 1.0  0.0;            # lattice vector one
        0.5  sqrt(3)/2]]     # lattice vector two

#put atom in reduced coordinate of lattice vector
atom = [ 1/3 1/3;  # atom one
         2/3 2/3]  # atom two

# honeycomb TightBinding
dim_k = 2 # dimension of k-space
dim_r = 2 # dimension of real-space
gr = tb_model(dim_k,dim_r,lat,orb)

#set onsite energy
Δ = 0.4
set_onsite!(gr,[-Δ,Δ])

#set hopping
t = -1.0 #hopping amplitude
# hopping from atom 1 to atom 2 in first BZ
set_hop!(gr,t,1,2,[0,0])
#hopping from atom 2 to atom 1 in neighbour BZ
set_hop!(gr,t,2,1,[1,0])
#hopping from atom 2 to atom 1 in neighbour BZ
set_hop!(gr,t,2,1,[0,1])

# see the model summary
show(gr)

# generate path in kspace (symmetry point) of honeycomb lattice
path = [[0,0],[2/3,1/3],[0.5,0.5],[0,0]]

#solve band
nk = 201
k_vec,k_dist,k_node = k_path(gr,path,nk)
eigvals = solve_eig(gr,k_vec)

#plot band
using Plots
plot(k_dist,eigvals')
```

See more advanced examples in 'examples' folder.

## Acknowledgement
All credits are goes to PythTB developer:
- Sinisa Coh (University of California at Riverside)
- David Vanderbilt (Rutgers University)

## TO DO
- documentation
- remake visualization of the TB model in Plots.jl (visualize.jl)
- transport calculation (calc.jl), i.e dc conductivity, optical conductivity

## Goal
PythTB.jl aims for fast, easy, and rich calculation of physical observable
based on Tight Binding approach.
