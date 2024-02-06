using Gurobi
using JuMP
using LinearAlgebra
using DataFrames, CSV
using DelimitedFiles
using Statistics
using Random
Random.seed!(123);


cd("C:\\Users\\alexn")
pwd()
s = Matrix(CSV.read("s.csv", DataFrame))
zᵖ = Matrix(CSV.read("z_p.csv", DataFrame) )
zᶜ = Matrix(CSV.read("z_c.csv", DataFrame) )
zᵐ = Matrix(CSV.read("z_m.csv", DataFrame) )

T = 10
H = 24
J = 12
Ω = 2
a1=[2, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2]
b1= a1 .+ 1
a2= 2 .+ a1
@show sum(a1)
@show sum(b1)
@show sum(a2)
a=hcat(a1,a1,a1,a1,a1,a1,a1,a1,a1,a1)
b=hcat(b1,b1,b1,b1,b1,b1,b1,b1,b1,b1) 
c=hcat(a2,a2,a2,a2,a2,a2,a2,a2,a2,a2)
D=fill(0.0,12,10,Ω)
#D[:, :, 1]=a
#D[:, :, 2]=b
#D[:, :, 3]=c

D[:, :, 1]=b
D[:, :, 2]=c

p = [1/Ω for ω in 1:Ω]
