# PythTB.jl

[![Build Status](https://travis-ci.org/muhammadgaffar/PythTB.jl.svg?branch=master)](https://travis-ci.org/muhammadgaffar/PythTB.jl)
[![Coverage Status](https://coveralls.io/repos/github/muhammadgaffar/PythTB.jl/badge.svg?branch=master)](https://coveralls.io/github/muhammadgaffar/PythTB.jl?branch=master)

PythTB.jl is Julia Wrapper for PythTB 1.7 Package in Python.

## Installation
```julia
] add https://github.com/muhammadgaffar/PythTB.jl
```

## Usage
This is example usage of Graphene TightBinding
```julia
using PythTB

#define lattice vectors for honeycomb lattice
lat = [[1,0],[0.5,sqrt(3)/2]]

#put orbital
orb = [[1/3,1/3],[2/3,2/3]]

#graphene model
dim_k = 2
dim_r = 2
gr = tb_model(dim_k,dim_r,lat,orb)

#set hopping
t = -1.0
set_hop!(gr,t,1,2,[0,0])
set_hop!(gr,t,2,1,[1,0])
set_hop!(gr,t,2,1,[0,1])

#model summary
show(gr)

#path in kspace (symmetry point)
path = [[0,0],[2/3,1/3],[0.5,0.5],[0,0]]
labels = [L"\Gamma","K","M",L"\Gamma"]

#solve band
nk = 201
k_vec,k_dist,k_node = k_path(gr,path,nk)
eigvals = solve_eig(gr,k_vec)

#plot band
using Plots
plot(k_dist,eigvals')
```
