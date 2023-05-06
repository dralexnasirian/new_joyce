using GLPK, Gurobi
using JuMP
using CSV, DataFrames
using Distributions, Random
using Printf
using LinearAlgebra
Random.seed!(123)

H = 2
J = 2
T = 2
Ω = 2
p = [1/Ω for ω in 1:Ω]

# s ##################
s = Diagonal( ones(H) )
s = Matrix(s)
# zᵖ ##################
l = range(100, stop=100 + (H-1)*10 , step=10)
zᵖ = collect(l)
# zᶜ ##################
l = range(200, stop=200 + (H-1)*10 , step=10)
zᶜ = collect(l)
# zᵐ ##################
l = fill(10, H)
zᵐ = Diagonal(l)

D1 = hcat(fill(1, H), fill(2, H))

D2 = hcat(fill(2, H), fill(3, H))
D = cat(D1, D2, dims=3)

@show s
@show zᵖ
@show zᶜ
@show zᵐ
@show D;

m = Model(Gurobi.Optimizer)
@variable(m , ψ[1:H], Bin)
@variable(m , χ[1:H , 1:J], Bin)
@variable(m , α[1:H , 1:J , 1:T , 1:Ω], Bin);
@variable(m , 0 ≤ γ[1:J , 1:T , 1:Ω], Int )

@objective(m , Min , 
    sum(T * zᵖ[h] * ψ[h] for h in 1:H ) + 
    sum(zᵐ[h,j] *χ[h,j] for h in 1:H for j in 1:J ) +
    sum( p[ω] * zᶜ[j] * γ[j,t,ω]  for j in 1:J for t in 1:T for ω in 1:Ω )
)


c1 = @constraint(m, [h in 1:H , j in 1:J], χ[h,j] ≤ ψ[h]  )
cc = @constraint(m, [h in 1:H , j in 1:J], s[h,j] + ψ[h] ≤ χ[h,j] + 1  )
c2 = @constraint(m, [h in 1:H , j in 1:J , t in 1:T , ω in 1:Ω], α[h,j,t,ω] ≤ χ[h,j] )
c3 = @constraint(m, [h in 1:H , t in 1:T , ω in 1:Ω], sum( α[h,j,t,ω] for j in 1:J) ≤ 1 )
c4 = @constraint(m, [j in 1:J , t in 1:T , ω in 1:Ω], sum( α[h,j,t,ω] for h in 1:H) + γ[j,t,ω] ≥ D[j,t,ω]);

optimize!(m)


no_var = num_variables(m)
no_con = num_constraints(m; count_variable_in_set_constraints = false);
time =  solve_time(m);
obj = objective_value(m);
ψᵏ = value.(ψ)
χᵏ = value.(χ)
γᵏ = value.(γ)
no_permanent =  sum(ψᵏ)
no_skills = sum(χᵏ)
ave_casual = sum(γᵏ) / (J*T*Ω)

DataFrame(
    H=H,
    J=J,
    T=T,
    Ω=Ω,
    no_var=no_var,
    no_con=no_con,
    time=time,
    obj=obj,
    no_permanent=no_permanent,
    no_multi=no_skills-no_permanent,
    ave_casual = ave_casual
    )

