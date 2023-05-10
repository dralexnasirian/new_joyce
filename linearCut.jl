using GLPK, Gurobi
using JuMP
using CSV, DataFrames
using Distributions, Random
using Printf
using LinearAlgebra
Random.seed!(123)

R = 200
H = 4
J = 4
T = 2
Ω = 1

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
l = Diagonal(l)
zᵐ = fill(10,H,J)
zᵐ = zᵐ - l

D = [2 1 3 3; 3 2 2 1 ]

D = reshape(D, J, T, Ω)

p = [1/Ω for ω in 1:Ω]

function index_generator(χᵏ)
    Hᵏ = [ [] for j in 1:J ] # for 1
    Hᵖ = [ [] for j in 1:J ] # for 0
    Hᵀ = [h for h in 1:H]

        for j in 1:J
        
            for h in 1:H
                if χᵏ[h,j] == 1
                    push!(Hᵏ[j] , h)
                end
            end
        end

        for j in 1:J
            Hᵖ[j] = setdiff( Hᵀ , Hᵏ[j] )
        end
   Jᵀ = [j for j in 1:J]
    Jᵏ = [ [] for h in 1:H]
    Jᵛ = [ [] for h in 1:H]
    for h in 1:H
        for j in 1:J
            if χᵏ[h,j] == 1
                push!(Jᵏ[h] , j)
            end
        end
    end
    for h in 1:H    
        Jᵛ[h] = setdiff(Jᵀ , Jᵏ[h])
    end
    return (Hᵏ , Hᵖ , Jᵏ , Jᵛ)
end

#######################################################################################
# f_zʳ
#######################################################################################
function f_zʳ(γᵏ) 
    s = [ [] for t in 1:T, ω in 1:Ω]
    for t in 1:T
        for ω in 1:Ω
            append!(s[t,ω] , [ones(Int(γᵏ[j,t,ω]) )  * zᶜ[j] for j in 1:J] )
            s[t,ω] = reduce(vcat , s[t,ω])
            s[t,ω] = sort(s[t,ω] , rev=true)
        end
    end
    for t in 1:T
        for ω in 1:Ω
            if length(s[t,ω]) < H
                append!(s[t,ω] , [0.0 for i in 1:(H - length(s[t,ω]) ) ] )
            end
        end
    end
    return s
end

#######################################################################################
# f_dʳ
#######################################################################################
function f_dʳ(zʳ)    
    d=[fill(0.0 , H) for t in 1:T, ω in 1:Ω]
    for t in 1:T
        for ω in 1:Ω
            for k in 1:H
                d[t,ω][k]= sum( zʳ[t,ω][l] for l in 1:k)
            end
        end
    end
    return d
end


#######################################################################################
# main problem
#######################################################################################
m = Model(Gurobi.Optimizer)
@variable(m, ψ[1:H], Bin)
@variable(m, χ[1:H , 1:J], Bin)
@variable(m, 0 ≤ β[1:T , 1:Ω])
# auxilary variable #########################



@objective(m, Min, 
    sum( T*zᵖ[h]*ψ[h] for h in 1:H ) +
    sum(χ[h,j] * zᵐ[h,j] for h in 1:H for j in 1:J) +
    sum( p[ω] * β[t,ω] for t in 1:T for ω in 1:Ω )
)
c1 = @constraint(m, [h in 1:H , j in 1:J], χ[h,j] ≤ ψ[h])
c2 = @constraint(m, [h in 1:H , j in 1:J], s[h,j] + ψ[h] ≤ χ[h,j] + 1);

#######################################################################################
# sub problem
#######################################################################################
function solve_sub(χᵏ)
    θ_t = fill(0.0,T,Ω)
    γ_t = fill(0.0,J,T,Ω)
    α_t = fill(0.0,H,J,T,Ω)
   sp = [Model(Gurobi.Optimizer) for t in 1:T, ω in 1:Ω] 
    for v in 1:T
        for u in 1:Ω
            @variable(sp[v,u], 0 ≤ γ[1:J , [v] , [u] ], Int)
            @variable(sp[v,u], α[1:H , 1:J , [v] , [u] ], Bin)
            @objective(sp[v,u], Min, sum(zᶜ[j] * γ[j,v,u] for j in 1:J ) )
            c1 = @constraint(sp[v,u] , [ h in 1:H , j in 1:J], α[h,j,v,u] ≤ χᵏ[h,j] )
            c2 = @constraint(sp[v,u] , [ h in 1:H ], sum( α[h,j,v,u] for j in 1:J ) <= 1 )
            c3 = @constraint(sp[v,u] , [ j in 1:J ], sum( α[h,j,v,u] for h in 1:H) + γ[j,v,u] ≥ D[j,v,u] )
            optimize!(sp[v,u])
            θ_t[v,u] = round.(objective_value(sp[v,u]), digits=0)
            γ_t[:,v,u] = round.(value.(γ) , digits=2)
            α_t[:,:,v,u] = value.(α) 
        end
    end
    return (θ_t , γ_t , α_t )
end


rSet = []
LBSet = []
UBSet = []
Θ_minSet = []
noPerSet = []
noMultiSet = []
costPerSet = []
costMultiSet = []
costβSet = []
costCasSet = []

function initiate()
   for r in 1:R
        optimize!(m) 
        ψʳ=value.(ψ)
        χʳ=value.(χ)
        βʳ=value.(β)
        
        
        ss = solve_sub(χʳ)
        θʳ = ss[1]  
        γʳ = ss[2]
        αʳ = ss[3]
        
        ig = index_generator(χʳ) 
        Hʳ = ig[1]
        Hᵖ = ig[2]
        Jʳ = ig[3]
        Jᵖ = ig[4]
        
        zʳ = f_zʳ(γʳ)

        dʳ = f_dʳ(zʳ)
        
        noPer = sum(ψʳ)
        noMulti = sum(χʳ) - sum(ψʳ)  
        costPer =  sum( T*zᵖ[h]*ψʳ[h] for h in 1:H ) 
        costMulti = sum(χʳ[h,j] * zᵐ[h,j] for h in 1:H for j in 1:J) 
        costβ = sum( p[ω] * βʳ[t,ω] for t in 1:T for ω in 1:Ω )
        costCas = sum(p[ω] * θʳ[t,ω] for t in 1:T for ω in 1:Ω)
        LB = objective_value(m)
        UB = costPer + costMulti + costCas
        
        push!(rSet, r) 
        push!(noPerSet, noPer) 
        push!(noMultiSet, noMulti) 
        push!(costPerSet,  costPer)
        push!(costMultiSet, costMulti)
        push!(costβSet, costβ)
        push!(costCasSet, costCas)
        push!(LBSet, LB) 
        push!(UBSet, UB) 
        Θ_min = minimum(UBSet)
        push!(Θ_minSet, Θ_min)
        
        
        
        dt =     DataFrame( r = rSet, 
                            noPer = noPerSet, 
                            noMulti = noMultiSet,
                            costPer = costPerSet,
                            costMulti = costMultiSet,
                            costβ = costβSet,
                            costCas = costCasSet,
                            LB = LBSet, 
                            UB = UBSet,
                            Θ_min = Θ_minSet
                            )
        
        @show dt
        print("\n \n ")
        @show DataFrame(χʳ, :auto)
        
        if LB ≥ Θ_min
            print("             OPTIMALITY                ")
            break
        end
        
        cʳ = @variable(m, [1:T, 1:Ω] )
        @constraint(m, [t in 1:T, ω in 1:Ω], cʳ[t,ω] ≥ 0 )
        
        @constraint(m, [t in 1:T, ω in 1:Ω], 
            β[t,ω] ≥ θʳ[t,ω] - cʳ[t,ω] - 
            θʳ[t,ω] * sum(1-χ[h,j] for j in 1:J for h in Hʳ[j]) )
        
        bʳ = @variable(m, [1:H], Bin) # bʳ is not defined for every t and ω
        
        @constraint(m, [h in 1:H],
                     J * bʳ[h] ≥ sum( χ[h,j] for j in Jᵖ[h] )   )
        
        @constraint(m, [h in 1:H],
                    bʳ[h] ≤ sum(χ[h,j] for j in Jᵖ[h]) )
        
        aʳ =  @variable(m, [k in 1:H], Bin)
       
        @constraint(m, sum( k * aʳ[k] for k in 1:H) == sum(bʳ[h] for h in 1:H) )
            
        @constraint(m, sum(aʳ[k] for k in 1:H) ≤ 1 ) # I think this should not be == and it should be ≤
        
        @constraint(m,  [t in 1:T, ω in 1:Ω], 
             cʳ[t,ω] ≤ sum( aʳ[k] * dʳ[t,ω][k] for k in 1:H  )   )
    
        
        @constraint(m, [t in 1:T, ω in 1:Ω], 
             cʳ[t,ω] ≤ sum(zʳ[t,ω][i] for i in 1:H)   )
    end
end


initiate()
