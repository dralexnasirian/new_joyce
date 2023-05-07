#using GLPK, Gurobi
#using JuMP
#using CSV, DataFrames
#using Printf
#using LinearAlgebra

##########################################
# parameters
##########################################
R = 1
H = 2
J = 2
T = 3
Ω = 1
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
l = Diagonal(l)
zᵐ = fill(10,H,J)
zᵐ = zᵐ - l

D = [
 3  0  3  3  3  0  3  1  2  1
 2  0  3  0  2  3  2  0  0  3
 3  2  2  1  2  0  2  0  3  1
 1  3  2  2  0  0  0  3  3  1
 0  2  1  3  2  0  3  0  3  0
 1  2  1  2  0  3  2  0  1  2
 0  1  0  2  1  0  2  3  2  0
 1  3  2  3  3  0  3  0  2  0
 2  0  1  0  2  1  3  2  2  2
 1  2  2  1  1  1  0  0  0  2
]

@show s
@show zᵖ
@show zᶜ
@show zᵐ
@show D;
###############################
K = 1
###############################

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
# main problem
#######################################################################################
mp = Model(GLPK.Optimizer)
@variable(mp, ψ[1:H], Bin)
@variable(mp, χ[1:H , 1:J], Bin)
@variable(mp, 0 ≤ β[1:T , 1:Ω])



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
    solveTimeSub = fill(0.0,T,Ω)
    noVarSub = fill(0.0,T,Ω)
    noConSub = fill(0.0,T,Ω)
    sp = [Model(GLPK.Optimizer) for t in 1:T, ω in 1:Ω] 
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
            solveTimeSub[v,u] = solve_time( sp[v,u] )
            noVarSub[v,u] = num_variables( sp[v,u] )
            noConSub[v,u] = num_constraints( sp[v,u] ; count_variable_in_set_constraints = false)
        end
    end
    return(θ_t , γ_t , α_t, solveTimeSub, noVarSub, noConSub )
end

#######################################################

#######################################################################################
# initiate
#######################################################################################

rSet = []

LBSet=[]
noPerSet=[]
noMultiSet=[]
perCostSet=[]
multiCostSet=[]
βSet=[]
noVarMainSet=[]
noConMainSet=[]
timeMainSet=[]

casCostSet=[]
UBSet=[]
Θ_minSet=[]
noVarSubSet=[]
noConSubSet=[]
timeSubSet=[]

timeTotalSet=[]

function initiate()
    for r in 1:R
        
        optimize!(mp)
        ψʳ = value.(ψ)
        χʳ = value.(χ)
        βʳ = value.(β)
        LB = objective_value(mp)
        
        ss = solve_sub(χʳ)
        θʳ = ss[1]
        γʳ = ss[2]
        αʳ = ss[3]
        timeSub = ss[4]
        noVarSub = ss[5]
        noConSub = ss[6]
        
        noPer = sum(ψʳ)
        noMulti = sum(χʳ) - sum(ψʳ)
        perCost = sum( zᵖ[h] * ψʳ[h] for h in 1:H )
        multiCost = sum(zᵐ[h,j] * χʳ[h,j] for h in 1:H for j in 1:J)
        βʳ = sum( βʳ[t,ω] for t in 1:T for ω in 1:Ω )
        noVarMain = num_variables(mp)
        noConMain = num_constraints(mp; count_variable_in_set_constraints = false )
        timeMain = solve_time(mp)
        
        casCost = sum(θʳ[t,ω] for t in 1:T for ω in 1:Ω)
        UB = perCost + multiCost + casCost
        timeSub = sum(timeSub[t,ω] for t in 1:T for ω in 1:Ω)
        noVarSub = sum( noVarSub[t,ω] for t in 1:T for ω in 1:Ω )
        noConSub = sum(noConSub[t,ω] for t in 1:T for ω in 1:Ω)
        
        timeTotal = timeMain + timeSub
        
        ig = index_generator(χʳ)
        Hʳ = ig[1]
        Hᵖ = ig[2]
        Jʳ = ig[3]
        Jᵖ = ig[4]
        
        zʳ = f_zʳ(γʳ)
        
        push!(rSet, r)
        
        push!(LBSet, LB)
        push!(noPerSet, noPer)
        push!(noMultiSet, noMulti)
        push!(perCostSet, perCost)
        push!(multiCostSet, multiCost)
        push!(βSet, βʳ)
        push!(noVarMainSet, noVarMain)
        push!(noConMainSet, noConMain)
        push!(timeMainSet, timeMain)
        
        push!(casCostSet, casCost)
        push!(UBSet, UB)
        Θ_min = minimum(UBSet)
        push!(Θ_minSet, Θ_min)
        push!(timeSubSet, timeSub)
        push!(noVarSubSet, noVarSub)
        push!(noConSubSet, noConSub)
        
        push!(timeTotalSet, timeTotal)
        
        for j in 1:J
            print( "Hʳ[$j] = ", Hʳ[j])
            print(" \n \n \n ")
        end
        
        for t in 1:T
            for ω in 1:Ω
                print("zʳ[$t,$ω] = ", zʳ[t,ω] )
                print("\n\n\n")
            end
        end
        
        df = DataFrame(r = rSet, 
                       LB = LBSet,
                       noPer = noPerSet, 
                       noMulti=noMultiSet,
                       perCost = perCostSet,
                       multiCost = multiCostSet,
                       β = βSet,
                       noVarMain = noVarMainSet,
                       noConMain = noConMainSet,
                       timeMain = timeMainSet,
                       casCost = casCostSet,
                       UB = UBSet,
                       Θ_min = Θ_minSet,
                       noVarSub = noVarSubSet,
                       noConSub = noConSubSet,
                       timeSub = timeSubSet,
                       timeTotal = timeTotalSet
                       )
        
        print("$df");
        print("\n \n ")
        
        @variable(mp, y[1:T, 1:Ω, 1:H], Bin)
        @variable(mp, δ[t=1:T, ω=1:Ω, 1:length( zʳ[t,ω] ) ], Bin)
        
        l5 = @constraint(mp, [t in 1:T, ω in 1:Ω], β[t,ω] ≥ θʳ[t,ω] - 
                                            sum(δ[t,ω,l] * zʳ[t,ω][l] for l in 1:length(zʳ[t,ω])) -
                                            θʳ[t,ω] * sum(1-χ[h,j] for j in 1:J for h in Hʳ[j])
                        )
        
        l4 = @constraint(mp, [t in 1:T, ω in 1:Ω, l in 2:length(zʳ[t,ω]) ],
                        δ[t,ω,l-1] ≥ δ[t,ω,l]
                        )
        
        l3 = @constraint(mp, [t in 1:T, ω in 1:Ω],
                           sum(y[t,ω,h] for h in 1:H) == sum(δ[t,ω,l] for l in 1:length(zʳ[t,ω]) )
                        )
        l2 = @constraint(mp, [t in 1:T, ω in 1:Ω, h in 1:H],
                        sum(χ[h,j] - χʳ[h,j] for j in 1:J ) ≤ J * y[t,ω,h]
                        )
        l1 = @constraint(mp, [t in 1:T, ω in 1:Ω, h in 1:H],
                        sum( χ[h,j] - χʳ[h,j] for j in 1:J ) ≥ 1 - (1 - y[t,ω,h])
                        )
        
        print("\n **********    y    ********  \n")
        for t in 1:T
            for ω in 1:Ω
                for h in 1:H
                    print("y[$t,$ω,$h]")
                    print("\n")
                end
            end
        end
        print("\n **********    δ    ********  \n")
        for t in 1:T
            for ω in 1:Ω
                for l in 1:length(zʳ[t,ω])
                    print("δ[$t,$ω,$l]")
                    print("\n")
                end
            end
        end
        
        for t in 1:T
            for ω in 1:Ω
                for l in 2:length(zʳ[t,ω])
                    for h in 1:H
                        print(" \n  \n")
                        print("l5[$t,$ω] = ", l5[t,ω] )
                        print(" \n  \n")
                        print("l4[$t,$ω,$l] = ", l4[t,ω,l] )
                        print(" \n  \n")
                        print("l3[$t,$ω] = ", l3[t,ω]) 
                        print(" \n  \n")
                        print("l2[$t,$ω,$h] = ", l2[t,ω,h]) 
                        print(" \n  \n")
                        print("l1[$t,$ω,$h] = ", l1[t,ω,h]) 
                    end
                end
            end
        end
        
        
    end
end

initiate()



