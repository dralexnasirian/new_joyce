using Gurobi
using JuMP
using LinearAlgebra
using DataFrames, CSV
using DelimitedFiles
using Statistics
using Random
# Random.seed!(123);

#############################################
# functions for parameters
#############################################
function vertical_concatenate(s, n)
        result = copy(s)  # Create a copy of s to avoid modifying the original matrix
        for i in 2:n
            result = vcat(result, s)
        end
        return result
end

    function f_s(H, J)
        l = length(H) / length(J)
        l = ceil(l)
        s = ones(length(J))
        s = Diagonal(s)
        s = Matrix(s)
        s = vertical_concatenate(s::Matrix{Float64}, l)
        s = s[1:length(H), 1:length(J)]
        return s
    end

set_χ = []    
set_D = []
      
set_mean_D = []
set_mean_zp = []
set_mean_zc = []
set_mean_zm = []

set_std_D = []
set_std_zp = []
set_std_zc = []
set_std_zm = []

set_H= []
set_J= []
set_T = []
set_O= []
set_U= []
set_noPer= []
set_noSkil = []
set_noCas = []
set_costPer= []
set_costCas= []
set_costTrain= []
set_costTot = []
result = []

num_iterations = 1_000
for i in 1:num_iterations
    global H = collect(1:rand(3:10))
    global J = collect(1:rand(3:length(H)))
    global s = f_s(H, J) 
    global T = collect(1:rand(200:300))
    global Ω = collect(1:rand(1:1))
    global p = [1/length(Ω) for ω in 1:length(Ω)]
    global U = rand(1:2)
    global D = rand(0:U,length(J), length(T), length(Ω))
    global zᵖ = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    global zᶜ = (1+ rand() ) * zᵖ
    a = rand(201:700)
    b = a - 200
    global zᵐ = rand(b:a, length(H), length(J))
    l = length(H) / length(J)
    L = floor(l)
    remaining = length(H) - ( L * length(J) )
    for l in 1:Int(L)                                                                                                                                                                                                                                                  
        for j in 1:length(J)                                                                                                                                                                                                                                      
            zᵐ[Int((l-1) * length(J) + j), j] = 0                                                                                                                                                                                                                  
        end                                                                                                                                                                                                                                                       
    end
    for i in 1:Int(remaining)
        zᵐ[Int(L*length(J) + i), i] = 0
    end
    ##############################
    deter = Model(Gurobi.Optimizer)
    @variable(deter , ψ[1:length(H)], Bin)
    @variable(deter , χ[1:length(H), 1:length(J)], Bin)
    @variable(deter , α[1:length(H) , 1:length(J) , 1:length(T), 1:length(Ω) ], Bin)
    @variable(deter , 0 ≤ γ[1:length(J), 1:length(T), 1:length(Ω)], Int )
    @objective(deter , Min ,
        sum( length(T) * zᵖ[h] * ψ[h] for h in H ) +
        sum(zᵐ[h,j] * χ[h,j] for h in H for j in J ) +
        sum( p[ω] * zᶜ[j] * γ[j,t,ω] for t in T for j in J for ω in Ω )
    )
    #############################
    # flexibility

    c1 = @constraint(deter, [h in H , j in J], χ[h,j] ≤ ψ[h]  )
    c2 = @constraint(deter, [h in H , j in J], s[h,j] + ψ[h] ≤ χ[h,j] + 1  )
    c3 = @constraint(deter, [h in H , j in J, t in T, ω in Ω ], α[h,j,t,ω] ≤ χ[h,j] )
    c4 = @constraint(deter, [h in H , t in T, ω in Ω], sum( α[h,j,t,ω] for j in J) ≤ 1 )
    c5 = @constraint(deter, [j in J , t in T, ω in Ω], sum( α[h,j,t,ω] for h in H) + γ[j,t,ω] ≥ D[j,t,ω])
    optimize!(deter)
    χᵏ = value.(χ)
    χᵏ = Int.(round.(χᵏ))
    dt_χᵏ = DataFrame(χᵏ, :auto)
    print("\n\n ***********\n dt_χᵏ \n *********** \n ", dt_χᵏ, "\n\n" )
    ψᵏ = value.(ψ)
    γᵏ = value.(γ)
    no_per_deter = sum(ψᵏ)
    no_skill_deter = sum(χᵏ)
    no_train_deter = no_skill_deter - no_per_deter
    no_aveCas_deter = sum(p[ω] * γᵏ[j,t,ω]  for j in J for t in T for ω in Ω ) / (length(J) * length(T) )
    no_aveCas_deter = round(no_aveCas_deter, digits = 2)
    cost_permanent_deter = sum(length(T) * zᵖ[h] * ψᵏ[h] for h in H )
    cost_multi_deter = sum(zᵐ[h,j] *χᵏ[h,j] for h in H for j in J )
    cost_cacual_deter = sum( p[ω] * zᶜ[j] * γᵏ[j,t,ω]  for j in J for t in T for ω in Ω )
    cost_total_deter = objective_value(deter)
    time_deter = solve_time(deter)
    noVar_deter = num_variables(deter)
    noCon_deter = num_constraints(deter; count_variable_in_set_constraints = false)

    push!(set_χ, dt_χᵏ)
    push!(set_D, D)

    push!(set_mean_D, mean(D) )
    push!(set_mean_zp, mean(zᵖ) )
    push!(set_mean_zc, mean(zᶜ) )
    push!(set_mean_zm, mean(zᵐ) )
    push!(set_std_D, std(D) )
    push!(set_std_zp, std(zᵖ) )
    push!(set_std_zc, std(zᶜ) )
    push!(set_std_zm, std(zᵐ) )

    push!(set_H, length(H) )
    push!(set_J, length(J))
    push!(set_T, length(T))
    push!(set_O, length(Ω))
    push!(set_U, U)
    push!(set_noPer, sum(ψᵏ))
    push!(set_noSkil, sum(χᵏ) - sum(ψᵏ))
    push!(set_noCas, no_aveCas_deter)
    push!(set_costPer, Int(sum(length(T) * zᵖ[h] * ψᵏ[h] for h in H)) )
    push!(set_costCas, round(sum(p[ω] * γᵏ[j,t,ω]  for j in J for t in T for ω in Ω ) , digits=2))
    push!(set_costTrain, Int(sum(zᵐ[h,j] * χᵏ[h,j] for h in H for j in J)))
    push!(set_costTot, Int(round(objective_value(deter),digits=0) ) )

    df = DataFrame(
        mean_D = set_mean_D,
        mean_zp = set_mean_zp,
        mean_zc = set_mean_zc,
        mean_zm = set_mean_zm,
        std_D = set_std_D,
        std_zp = set_std_zp,
        std_zc = set_std_zc,
        std_zm = set_std_zm,
        H = set_H,
        J = set_J,
        T = set_T,
        O = set_O,
        U = set_U,
        noPer = set_noPer,
        noSkill = set_noSkil,
        noCas_perSkill = set_noCas,
        costPer = set_costPer,
        costCas = set_costCas,
        costTrain = set_costTrain,
        costTot = set_costTot
    )
    # cd("C:\\Users\\alexn\\OneDrive - RMIT University\\0. Araz + Joyce\\codes\\21. General managerial insights")
    CSV.write("df_.csv", df)
end

#set_χ[41]


#set_D[41]
#set_D[19]
