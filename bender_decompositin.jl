using GLPK, Gurobi
using JuMP
using CSV, DataFrames
using Printf
using LinearAlgebra
##########################################
# parameters
##########################################


H = 5
J = 5
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
#################
D1 = hcat(fill(1, H), fill(2, H))

D2 = hcat(fill(2, H), fill(3, H))
D = cat(D1, D2, dims=3)

@show s
@show zᵖ
@show zᶜ
@show zᵐ
@show D;
#######################
K = 10
L = 6
##################

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
    return s
end

#######################################################################################
# max_zʳ
#######################################################################################
function max_zʳ(zʳ)
    l = fill(0,T,Ω)
    for t in 1:T
        for ω in 1:Ω
            l[t,ω] = count(i->(i>0), zʳ[t,ω])
        end
    end
    return maximum(l)
end

#######################################################################################
# main problem
#######################################################################################
mp = Model(Gurobi.Optimizer)
@variable(mp, ψ[1:H], Bin)
@variable(mp, χ[1:H , 1:J], Bin)
@variable(mp, 0 ≤ β[1:T , 1:Ω])
# auxilary variable #########################



@objective(mp, Min, 
    sum( T*zᵖ[h]*ψ[h] for h in 1:H ) +
    sum(χ[h,j] * zᵐ[h,j] for h in 1:H for j in 1:J) +
    sum( p[ω] * β[t,ω] for t in 1:T for ω in 1:Ω )
)
c1 = @constraint(mp, [h in 1:H , j in 1:J], χ[h,j] ≤ ψ[h])
c2 = @constraint(mp, [h in 1:H , j in 1:J], s[h,j] + ψ[h] ≤ χ[h,j] + 1);



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
    return(θ_t , γ_t , α_t )
end

#######################################################

#######################################################################################
# initiate
#######################################################################################
function initiate()
    β_set = []
    β_t1_ω1_set = []
    β_t2_ω1_set = []
    β_set = []
    LB_set = []
    UB_set = []
    Θₘᵢₙ_set = []
    no_permanent_set = []
    no_skill_set = []
    no_multiskilling_set = []
    for k in 1:K

        optimize!(mp)

       
        @show k
     
        ψᵏ = value.(ψ)
        @show ψᵏ
        no_permanent = sum(ψᵏ)
        push!(no_permanent_set , no_permanent)
        χᵏ = value.(χ)
        no_skill = sum(χᵏ)
        push!(no_skill_set , no_skill)
        no_multiskilling = no_skill - no_permanent 
        push!(no_multiskilling_set , no_multiskilling)
        dt_χᵏ = DataFrame(χᵏ , :auto)
        @show dt_χᵏ
        βᵏ =value.(β)
        βᵏ = sum( p[ω] * βᵏ[t,ω] for t in 1:T for ω in 1:Ω )
        @show βᵏ
        push!(β_set, βᵏ)
        LB = objective_value(mp)
        @show LB
        
        ss = solve_sub(χᵏ)
  
        θᵏ = ss[1]
        γᵏ = ss[2]
        for ω in 1:Ω
            @show DataFrame(γᵏ[:,:,ω] , :auto)
        end
        αᵏ = ss[3]
        casual_cost = sum(p[ω] * zᶜ[j] * γᵏ[j,t,ω] for j in 1:J for t in 1:T for ω in 1:Ω )

        cost_permanent = sum( T*zᵖ[h]*ψᵏ[h] for h in 1:H )
        @show cost_permanent 
        cost_training = sum(χᵏ[h,j] * zᵐ[h,j] for h in 1:H for j in 1:J) 
        @show cost_training
        UB = cost_permanent + cost_training + casual_cost
        @show βᵏ
        
        @show LB

        @show UB
        push!(UB_set , UB)
        push!(LB_set , LB)
        Θₘᵢₙ = minimum(UB_set)
        push!(Θₘᵢₙ_set , Θₘᵢₙ)
        @show Θₘᵢₙ
        @show casual_cost
        @show θᵏ
        # ******************************* 1 
        dt = DataFrame(LB = LB_set , UB = UB_set, c_UB = Θₘᵢₙ_set , 
            num_per = no_permanent_set, num_skill = no_skill_set, num_multi = no_multiskilling_set, beta = β_set )
        @show dt
        # ******************************* 2
        ig = index_generator(χᵏ)
        Hᵏ = ig[1]
        Hᵖ = ig[2]
        Jᵏ = ig[3]
        Jᵖ = ig[4]
        @show Hᵏ
        @show Hᵖ
        @show Jᵏ
        @show Jᵖ

        if LB ≥ Θₘᵢₙ
            println("#######       OPTIMALITY      ##############")
            break
        end
       
        y = @variable(mp, [1:T , 1:Ω , 1:H , k:k], Bin)
        @show y
        
        zʳ = f_zʳ(γᵏ) 
        @show zʳ
        
        L = [length(zʳ[t,ω]) for t in 1:T, ω in 1:Ω]
       @show L
        
        δ = @variable(mp, [t = 1:T , ω = 1:Ω , 1:length(zʳ[t,ω]) , k:k], Bin)
        @show δ
        
         con1 = @constraint(mp, [t in 1:T , ω in 1:Ω , h in 1:H],
                                                sum(χ[h,j] - χᵏ[h,j] for j in 1:J) ≥ 1 - (1-y[t,ω,h,k]) * 1000
        )
        
        con2 = @constraint(mp, [t in 1:T , ω in 1:Ω , h in 1:H], 
                 sum(χ[h,j] - χᵏ[h,j] for j in 1:J) ≤ y[t,ω,h,k] * 1000
        )
        
                con3 = @constraint(mp, [t=1:T , ω=1:Ω , l = 2:length(zʳ[t,ω]) ], δ[t,ω,l-1,k] ≥ δ[t,ω,l,k] )
      
        

        con4 = @constraint(mp, [t in 1:T , ω in 1:Ω ], 
            sum(y[t,ω,h,k] for h in 1:H) == sum( δ[t,ω,l,k] for l in 1:length(zʳ[t,ω]) ) )
      
        

         con5 = @constraint(mp, [t in 1:T , ω in 1:Ω ], 
                                                β[t,ω] ≥ θᵏ[t,ω] - 
                                                           sum(  δ[t,ω,l,k] * zʳ[t,ω][l] for l in 1:length(zʳ[t,ω]) ) -
                                                           θᵏ[t,ω] * sum( (1 - χ[h,j] ) for j in 1:J for h in Hᵏ[j] ) ) 
        
     #@show num_constraints(mp; count_variable_in_set_constraints = false)
     #@show num_variables(mp)
    end
end

initiate()
