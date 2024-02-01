
using Gurobi
using JuMP
using LinearAlgebra
using DataFrames, CSV
using DelimitedFiles
using Statistics
using Random
Random.seed!(123);

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

function solve_sub(χʳ, zᵖ, zᶜ, zᵐ,s, H, J, T, D, Ω)
    θ_t = fill(0.0, length(T), length(Ω) )
    γ_t = fill(0.0, length(J), length(T), length(Ω))
    α_t = fill(0.0, length(H), length(J), length(T), length(Ω))
    solveTimeSub = fill(0.0, length(T), length(Ω))
    noVarSub = fill(0.0, length(T), length(Ω))
    noConSub = fill(0.0, length(T), length(Ω))
    for ω in 1:length(Ω)
        for t in 1:length(T)
            sp = Model(Gurobi.Optimizer)
            @variable(sp, 0 ≤ γ[J], Int)
            @variable(sp, α[H, J], Bin)
            @objective(sp, Min, sum(zᶜ[j] * γ[j] for j in J  ) )
            c1 = @constraint(sp, [ h in H, j in J], α[h,j] ≤ χʳ[h,j] )
            c2 = @constraint(sp, [ h in H], sum( α[h,j] for j in J) <= 1 )
            c3 = @constraint(sp, [ j in J], sum( α[h,j] for h in H) + γ[j] ≥ D[j,t,ω] )
            optimize!(sp)
            θ_t[t,ω] = round.(objective_value(sp), digits=0)
            γ_t[:,t,ω] = round.(value.(γ) , digits=2)
            α_t[:,:,t,ω] = value.(α)
            solveTimeSub[t,ω] = solve_time( sp )
            noVarSub[t,ω] = num_variables( sp )
            noConSub[t,ω] = num_constraints( sp ; count_variable_in_set_constraints = false)
        end
    end
    return(θ_t , γ_t, α_t, sum(solveTimeSub), sum(noVarSub), sum(noConSub) )
end

function f_SP_expectaion(χʳ,D)
    θ_t = fill(0.0, length(T), length(Ω) )
    γ_t = fill(0.0, length(J), length(T), length(Ω) )
    α_t = fill(0.0, length(H), length(J), length(T), length(Ω) )
    solveTimeSub = fill(0.0, length(T), length(Ω) )
    noVarSub = fill(0.0, length(T), length(Ω) )
    noConSub = fill(0.0, length(T), length(Ω) )
    for ω in Ω
        for t in T
            for j in J
                set_normalized_rhs(c3[j], D[j,t,ω] )
                for h in H
                    set_normalized_rhs(c1[h,j], χʳ[h,j])
                end
            end
        optimize!(sp)
        θ_t[t,ω] = round.(objective_value(sp), digits=0)
        γ_t[:,t,ω] = round.(value.(γ) , digits=2)
        α_t[:,:,t,ω] = value.(α)
        solveTimeSub[t,ω] = solve_time( sp )
        noVarSub[t,ω] = num_variables( sp )
        noConSub[t,ω] = num_constraints( sp ; count_variable_in_set_constraints = false)
        end
    end
    return(θ_t , γ_t, α_t, sum(solveTimeSub), sum(noVarSub), sum(noConSub) )
end

function f_D_expectation()
    D = fill(0, length(J), length(T), length(Ω) )
    l = []
    for a in 0:1
        for b in 0:1
            for c in 0:1
                for d in 0:1
                    for e in 0:1
                        for f in 0:1
                            push!(l, [a, b, c, d, e, f])
                        end
                    end
                end
            end
        end
    end
    l = hcat(l...)
    for i in T
        D[:,i,:] = l
    end
    return D
end

function f_H¹_H⁰(H, J, χʳ)
    H¹ = [ [] for j in J]
    H⁰ = [ [] for j in J]
    for h in H
        for j in J
            if χʳ[h,j] == 1
                push!(H¹[j], h)
            else
                push!(H⁰[j], h)
            end
        end
    end
    return H¹, H⁰
end

function f_J¹_J⁰(H, J, χʳ)
    J¹ = [ [] for h in H]
    J⁰ = [ [] for h in H]
    for j in J
        for h in H
            if χʳ[h,j] == 1
                push!(J¹[h], j)
            else
                push!(J⁰[h], j)
            end
        end
    end
    return J¹, J⁰
end

function f_df_info(r, χʳ, βʳ, θʳ, γʳ, set_LB, set_UB, set_UB_min, UB)
    noPer = sum(χʳ)
    costPer = sum( length(T) * zᵖ[j] * χʳ[h,j] * s[h,j] for j in J for h in H )
    costTrain = sum(zᵐ[h,j] * χʳ[h,j] for h in H for j in J) 
    
    push!(set_r, r)
    push!(set_costPer, costPer)
    push!(set_costTrain, costTrain)
    push!(set_costTot, costPer + costTrain + UB)

    df_χʳ = DataFrame( Int.(round.(χʳ)), :auto)
    push!(set_χ, df_χʳ)

    # df_υʳ = DataFrame( Int.(υʳ), :auto)
    # df_δʳ = DataFrame( Int.(δʳ), :auto)

    # df_zʳ = DataFrame(zʳ, :auto)
    # df_θʳ = DataFrame(θʳ, :auto)
    # df_γʳ = DataFrame( reshape(Int.(γʳ), length(J), length(T)) , :auto)
   
    df = DataFrame(
        r = set_r,
        costPer = set_costPer,
        costTrain = set_costTrain,
        costTot = set_costTot,
        LB = set_LB,
        UB = set_UB,
        UB_min = set_UB_min
    )
    print("\n\n ***********\n   r   \n *********** \n $(r)\n\n")
    print("\n\n ***********\n df_χʳ \n *********** \n $(df_χʳ)\n\n")

    # print("\n\n ***********\n df_υʳ \n *********** \n $(df_υʳ)\n\n")

    # print("\n\n ***********\n df_δʳ \n *********** \n $(df_δʳ)\n\n")

    # print("\n\n ***********\n df_θʳ \n *********** \n $(df_θʳ)\n\n")

    # print("\n\n ***********\n df_γʳ \n *********** \n $(df_γʳ)\n\n")

    # print("\n\n ***********\n df_zʳ \n *********** \n $(df_zʳ)\n\n")

    print("\n\n ***********\n df_info \n *********** \n $(df)\n\n")
end

function f_Jʳ(γʳ, zᵖ, zᶜ, zᵐ,s, H, J, T, D, Ω)
    qʳ = [[] for t in 1:length(T), ω in 1:length(Ω) ]
    for ω in Ω
        for t in T
            for j in J
                if γʳ[j,t,ω] ≥ 1
                    append!(qʳ[t,ω], j)
                end
            end
        end
    end
    return qʳ
end

function f_zʳ(H, J, T, Ω, γʳ, zᶜ)
    l = [ fill(zᶜ[j], Int(γʳ[j,t,ω]) ) for j in J, t in T, ω in Ω]
    zʳ = [ [] for t in T, ω in Ω ]
    for t in T
        for ω in Ω
            zʳ[t,ω] = sort(append!( reduce(vcat, l[:,t,ω]) ),rev=true)
            if length(zʳ[t,ω]) < length(H)
                append!( zʳ[t,ω], fill(0, length(H) - length(zʳ[t,ω] ) ) )
            elseif length(zʳ[t,ω]) > length(H)
                zʳ[t,ω] = zʳ[t,ω][1:length(H)]
            end
        end
    end
    return zʳ
end

function f_sets()
    global set_r = []
    global set_χ = []
    global set_LB = []
    global set_UB = []
    global set_UB_min = []
    global set_r = []
    global set_noPer = []
    global set_noTrain = []
    global set_costPer = [] 
    global set_costTrain = [] 
    global set_costCas = []
    global set_costTot = []
    global set_noVarMP = []
    global set_noConMP = []
    global set_timeMP = []
    global set_noVarSP = []
    global set_noConSP = []
    global set_timeSP = []
end

function f_indicator(x)
    if round(x) ≤ 0
        return 0
    elseif round(x) > 0
        return 1
    end 
end

function f_U(γʳ,zᶜ,T,Ω,J)
    U = [ [] for t in T, ω in Ω]
    for t in T
        for ω in Ω
            for j in J
                push!(U[t,ω], f_indicator(γʳ[j,t,ω]) * zᶜ[j] )
            end
        end
    end
    return U 
end


function f_I(γʳ,zᶜ,T,Ω,J)
    I = [ [] for t in T, ω in Ω]
    l = f_indicator.(γʳ)
    for t in T
        for ω in Ω
            for j in J
                if l[j,t,ω] > 0
                    push!(I[t,ω], j)
                end
            end
        end
    end
    sort!.(I, rev=true)
    return I
end


H = collect(1:6)
J = collect(1:6)
s = f_s(H, J) 
T = collect(1:20)
Ω = collect(1:40_000)
p = [1/length(Ω) for ω in 1:length(Ω)]
U = 1
R = 50

D = rand(0:U,length(J), length(T), length(Ω))
#D = [1 0 1 1;
#     1 0 0 1;
#     0 1 0 0;
#     0 0 0 1]
#D = reshape(D, length(J) , length(T) , length(Ω) ) 

zᵖ = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
zᵖ = zᵖ ./ 2
zᶜ = zᵖ .* 1.5
zᵐ = [0 120 130 140 150 160 170 180 190 200;
    210 0   230 240 250 260 270 280 290 300;
    310 320 0   340 350 360 370 380 390 400;
    410 420 430 0   450 460 470 480 490 500;
    510 520 530 540 0 560 570 580 590 600;
    610 620 630 640 640 0   670 680 690 700;
    710 720 730 740 750 760 0   780 790 800;
    810 820 830 840 850 860 870 0   890 900;
    1010 1020 1030 1040 1050 1060 1070 1080 0 1100;
    1110 1120 1130 1140 1150 1160 1170 1180 1190 0 ]
s = f_s(H, J)   
zᵖ = zᵖ[1:length(H)]
zᶜ = zᶜ[1:length(J)]
zᵐ = zᵐ[1:length(H), 1:length(J)]



#############################################
# Deterministic
#############################################

elapsed_time = @time begin   
    deter = Model(Gurobi.Optimizer)
    @variable(deter , ψ[1:length(H)], Bin)
    @variable(deter , χ[1:length(H), 1:length(J)], Bin)
    @variable(deter , α[1:length(H) , 1:length(J) , 1:length(T), 1:length(Ω) ], Bin)
    @variable(deter , 0 ≤ γ[1:length(J), 1:length(T), 1:length(Ω)], Int )
    @objective(deter , Min ,
        sum( length(T) * zᵖ[h] * ψ[h] for h in H ) +
        sum(zᵐ[h,j] * χ[h,j] for h in H for j in J ) +
        sum( p[ω] * zᶜ[j] * γ[j,t,ω]  for j in J for t in T for ω in Ω )
    )
    #############################
    # flexibility
    # @constraint(deter, [h in H], sum( χ[h,j] for j in J) ≤ flexibility )

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
    no_aveCas_deter = sum(p[ω] * γᵏ[j,t,ω]  for j in J for t in T for ω in Ω ) / (length(T) * length(Ω))
    no_aveCas_deter = round(no_aveCas_deter, digits = 2)
    cost_permanent_deter = sum(length(T) * zᵖ[h] * ψᵏ[h] for h in H )
    cost_multi_deter = sum(zᵐ[h,j] *χᵏ[h,j] for h in H for j in J )
    cost_cacual_deter = sum( p[ω] * zᶜ[j] * γᵏ[j,t,ω]  for j in J for t in T for ω in Ω )
    cost_total_deter = objective_value(deter)
    time_deter = solve_time(deter)
    noVar_deter = num_variables(deter)
    noCon_deter = num_constraints(deter; count_variable_in_set_constraints = false)

    df_para = DataFrame(
        H = length(H), 
        J = length(J), 
        T = length(T), 
        Ω = length(Ω), 
        U = U
    )

    print("\n\n df_para: \n\n $(df_para)\n\n")

    df_result = DataFrame(
        NoPer = no_per_deter, 
        NoTrain = no_train_deter, 
        NoCas = no_aveCas_deter, 
        CostPer = Int(round(cost_permanent_deter)),
        CostTrain = Int(round(cost_multi_deter)), 
        CostCas = Int(round(cost_cacual_deter)), 
        CostTot = Int(round(cost_total_deter)), 
        noVar = noVar_deter, 
        noCon = noCon_deter,
        Time = round(time_deter, digits = 3)
    )

    print("\n\n Results: \n\n $(df_result)\n\n")

end
#print("elapsed time = ", elapsed_time )

#######################################################
# cut 5 2 1
#######################################################
V = 2^length(H)
Ω = collect(1:V)
p = [1/V for i in 1:V]
D = f_D_expectation()

#########################################
# MP
m = Model(Gurobi.Optimizer)
set_silent(m)
set_optimizer_attribute(m, "OutputFlag", 0)
@variable(m, χ[H, J], Bin)
@variable(m, β[T, Ω] ≥ 0 )
@objective(m, Min, 
    sum( length(T) * zᵖ[j] * χ[h,j] * s[h,j] for j in J for h in H ) +
    sum(zᵐ[h,j] * χ[h,j] for h in H for j in J) + 
    sum( p[ω] * β[t,ω] for t in T for ω in Ω )
)
@constraint(m, [h in H, j in J], sum(s[h,l] * χ[h,l] for l in J ) ≥ χ[h,j] )
# @constraint(m, [h in H], sum( χ[h,j] for j in J) ≤ flexibility  )

f_sets()

#########################################
# SP
sp = Model(Gurobi.Optimizer)
set_silent(sp)
set_optimizer_attribute(sp, "OutputFlag", 0)
@variable(sp, 0 ≤ γ[J], Int)
@variable(sp, α[H, J], Bin)
@objective(sp, Min, sum(zᶜ[j] * γ[j] for j in J  ) )
x = fill(0, length(H), length(J))
@constraint(sp, c1[ h in H, j in J], α[h,j] ≤ x[h,j] )
@constraint(sp, c2[ h in H], sum( α[h,j] for j in J) <= 1 )
x = fill(0,length(J))
@constraint(sp, c3[ j in J], sum( α[h,j] for h in H) + γ[j] ≥ x[j] )


function initiate_cut521()
    for r in 1:R
        optimize!(m)
        χʳ = value.(χ)
        βʳ = value.(β)
        LB = sum( p[ω] * βʳ[t,ω] for t in T for ω in Ω )
        push!(set_LB, LB)
        timeMP = solve_time(m)
        noVarMP = num_variables(m)
        noConMP = num_constraints(m; count_variable_in_set_constraints = false)
        push!(set_noVarMP ,noVarMP) 
        push!(set_noConMP ,noConMP) 
        push!(set_timeMP ,timeMP )

        H¹, H⁰ = f_H¹_H⁰(H, J, χʳ)
        J¹, J⁰ = f_J¹_J⁰(H, J, χʳ)

        # θʳ, γʳ, αʳ, timeSP, noVarSP, noConSP = solve_sub(χʳ, zᵖ, zᶜ, zᵐ,s, H, J, T, D, Ω)
        θʳ, γʳ, αʳ, timeSP, noVarSP, noConSP = f_SP_expectaion(χʳ, D)

        UB = sum( p[ω] * θʳ[t,ω] for t in T for ω in Ω )
        push!(set_UB, UB)
        UB_min = minimum(set_UB)
        push!(set_UB_min, UB_min)

        push!(set_noVarSP, noVarSP) 
        push!(set_noConSP, noConSP)
        push!(set_timeSP, timeSP)

        zʳ = f_zʳ(H, J, T, Ω, γʳ, zᶜ)
        Jʳ = f_Jʳ(γʳ, zᵖ, zᶜ, zᵐ,s, H, J, T, D, Ω)    

        f_df_info(r, χʳ, βʳ, θʳ, γʳ, set_LB, set_UB, set_UB_min, UB)
        
        if LB + 1 ≥ UB_min
            print("\n\n******************\nOPTIMALITY\n******************\n\n")
            df = DataFrame(
                noVarMP = sum(set_noVarMP),
                noConMP = sum(set_noConMP),
                timeMP = sum(set_timeMP),
                noVarSP = sum(set_noVarSP),
                noConSP = sum(set_noConSP),
                timeSP = sum(set_timeSP),
                r = r,
                noVar = Int(sum(set_noVarMP) + sum(set_noVarSP)),
                noCon = Int(sum(set_noConMP) + sum(set_noConSP)),
                time = sum(set_timeMP) + sum(set_timeSP)
            )
            print(df)
            break
        end

        ###################################
        U = f_U(γʳ,zᶜ,T,Ω,J)
        I = f_I(γʳ,zᶜ,T,Ω,J)

        f = @variable(m, [H, J], Bin)
        y = @variable(m, [T, Ω, H, J], Bin)

        c_36 = @constraint(m, [h in H, j in J⁰[h] ], f[h,j] == χ[h,j] )
        c_37 = @constraint(m, [t in T, ω in Ω, h in H, i in I[t,ω]], y[t,ω,h,i] ≤ f[h,i] )
        c_38 = @constraint(m, [t in T, ω in Ω, h in H], sum(y[t,ω,h,i] for i in I[t,ω]) ≤ 1 )
        c_39 = @constraint(m, [t in T, ω in Ω, h in H], sum(y[t,ω,h,i] for i in I[t,ω]) ≤ sum(f[h,i] for i in I[t,ω]) )
        c_40 = @constraint(m, [t in T, ω in Ω, h in H], length(J)*sum(y[t,ω,h,i] for i in I[t,ω]) ≥ sum(f[h,i] for i in I[t,ω]) )
        c_42 = @constraint(m, 
        sum( p[ω] * β[t,ω] for t in T for ω in Ω )
        ≥ 
        sum(p[ω]*γʳ[j,t,ω]*zᶜ[j] for t in T for ω in Ω for j in J) +
        sum( (1 - χ[h,j]) * sum( p[ω]*αʳ[h,j,t,ω]*zᶜ[j] for t in T for ω in Ω) for h in H for j in J¹[h] )  - 
        sum(y[t,ω,h,i] * U[t,ω][i] for t in T for ω in Ω for h in H for i in I[t,ω])
        )

        #print("\n c_42 = $c_42 \n")
        #########################################
        # cut I
        υ = @variable(m, [H], Bin )
        δ = @variable(m, [H], Bin )
         con_27_a = @constraint(m, [h in H], sum(χ[h,j] for j in J⁰[h]) ≥ υ[h] )
         con_27_b = @constraint(m, [h in H], sum(χ[h,j] for j in J⁰[h]) ≤ length(J⁰[h]) * υ[h] )
         con_27_c = @constraint(m, sum(δ[h] for h in H) == sum( υ[h] for h in H ) )
         con_27_d = @constraint(m, [l in 2:length(H)], δ[l-1] ≥ δ[l] )    
         con_27_e = @constraint(m, [t in T, ω in Ω], β[t,ω] ≥ θʳ[t,ω] -   
             sum( zʳ[t,ω][l] * δ[l] for l in 1:length(zʳ[t,ω])  ) - 
             θʳ[t,ω] * sum(1- χ[h,j] for j in J for h in H¹[j] )
         )
        #########################################
        # cut II
        λ = @variable(m, [T,Ω,J] )
        con_34_a = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω] ], λ[t,ω,j] ≤ γʳ[j,t,ω] * zᶜ[j] )
        con_34_b = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω] ], λ[t,ω,j] ≤ sum( χ[h,j] * zᶜ[j] for h in H⁰[j] ) )
        cʳ = [isempty(Jʳ[t, ω]) ? 0 : sum(λ[t, ω, j] for j in Jʳ[t, ω]) for t in T, ω in Ω]
         con_34_c = @constraint(m, [t in T, ω in Ω], β[t,ω] ≥ θʳ[t,ω] -   
            cʳ[t,ω] - 
            θʳ[t,ω] * sum(1- χ[h,j] for j in J for h in H¹[j] )
        )

        #########################################
         # cut III
        # ϕ = fill(0.0, length(T), length(Ω), length(J), length(H) )
        # # for updating ϕ
        # for t in T
        #     for ω in Ω
        #         for j in Jʳ[t,ω]
        #             for h in H⁰[j]
        #                 if ( zᶜ[j] - sum( αʳ[h,l,t,ω] * zᶜ[l] for l in J) ) > 0
        #                     ϕ[t,ω,j,h] = zᶜ[j] - sum( αʳ[h,l,t,ω] * zᶜ[l] for l in J)
        #                 end
        #             end
        #         end
        #     end
        # end
        # ϵ = @variable(m, [T, Ω, J])
        # con_29_a = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω]], ϵ[t,ω,j] ≤ sum(χ[h,j] * ϕ[t,ω,j,h] for h in H⁰[j] ) )   
        # con_29_b = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω]], ϵ[t,ω,j] ≤ γʳ[j,t,ω] * zᶜ[j] )   
        # con_29_c = @constraint(m, [t in T, ω in Ω], β[t,ω] ≥ θʳ[t,ω] -
        #     sum( ϵ[t,ω,j] for j in Jʳ[t,ω] ) - 
        #     θʳ[t,ω] * sum(1 - χ[h,j] for j in J for h in H¹[j] )
        # )
    end
end
print("\n\n------------\n\n cut 5 cut 2 cut 1 \n\n------------\n\n")
initiate_cut521()

#######################################################
# cut 5 3 1
#######################################################
V = 2^length(H)
Ω = collect(1:V)
p = [1/V for i in 1:V]
D = f_D_expectation()

#########################################
# MP
m = Model(Gurobi.Optimizer)
set_silent(m)
set_optimizer_attribute(m, "OutputFlag", 0)
@variable(m, χ[H, J], Bin)
@variable(m, β[T, Ω] ≥ 0 )
@objective(m, Min, 
    sum( length(T) * zᵖ[j] * χ[h,j] * s[h,j] for j in J for h in H ) +
    sum(zᵐ[h,j] * χ[h,j] for h in H for j in J) + 
    sum( p[ω] * β[t,ω] for t in T for ω in Ω )
)
@constraint(m, [h in H, j in J], sum(s[h,l] * χ[h,l] for l in J ) ≥ χ[h,j] )
# @constraint(m, [h in H], sum( χ[h,j] for j in J) ≤ flexibility  )

f_sets()

#########################################
# SP
sp = Model(Gurobi.Optimizer)
set_silent(sp)
set_optimizer_attribute(sp, "OutputFlag", 0)
@variable(sp, 0 ≤ γ[J], Int)
@variable(sp, α[H, J], Bin)
@objective(sp, Min, sum(zᶜ[j] * γ[j] for j in J  ) )
x = fill(0, length(H), length(J))
@constraint(sp, c1[ h in H, j in J], α[h,j] ≤ x[h,j] )
@constraint(sp, c2[ h in H], sum( α[h,j] for j in J) <= 1 )
x = fill(0,length(J))
@constraint(sp, c3[ j in J], sum( α[h,j] for h in H) + γ[j] ≥ x[j] )

print("\n\n------------\n\n cut 5 cut 3 cut 1 \n\n------------\n\n")
function initiate_cut531()
    for r in 1:R
        optimize!(m)
        χʳ = value.(χ)
        βʳ = value.(β)
        LB = sum( p[ω] * βʳ[t,ω] for t in T for ω in Ω )
        push!(set_LB, LB)
        timeMP = solve_time(m)
        noVarMP = num_variables(m)
        noConMP = num_constraints(m; count_variable_in_set_constraints = false)
        push!(set_noVarMP ,noVarMP) 
        push!(set_noConMP ,noConMP) 
        push!(set_timeMP ,timeMP )

        H¹, H⁰ = f_H¹_H⁰(H, J, χʳ)
        J¹, J⁰ = f_J¹_J⁰(H, J, χʳ)

        # θʳ, γʳ, αʳ, timeSP, noVarSP, noConSP = solve_sub(χʳ, zᵖ, zᶜ, zᵐ,s, H, J, T, D, Ω)
        θʳ, γʳ, αʳ, timeSP, noVarSP, noConSP = f_SP_expectaion(χʳ, D)

        UB = sum( p[ω] * θʳ[t,ω] for t in T for ω in Ω )
        push!(set_UB, UB)
        UB_min = minimum(set_UB)
        push!(set_UB_min, UB_min)

        push!(set_noVarSP, noVarSP) 
        push!(set_noConSP, noConSP)
        push!(set_timeSP, timeSP)

        zʳ = f_zʳ(H, J, T, Ω, γʳ, zᶜ)
        Jʳ = f_Jʳ(γʳ, zᵖ, zᶜ, zᵐ,s, H, J, T, D, Ω)    

        f_df_info(r, χʳ, βʳ, θʳ, γʳ, set_LB, set_UB, set_UB_min, UB)
        
        if LB + 1 ≥ UB_min
            print("\n\n******************\nOPTIMALITY\n******************\n\n")
            df = DataFrame(
                noVarMP = sum(set_noVarMP),
                noConMP = sum(set_noConMP),
                timeMP = sum(set_timeMP),
                noVarSP = sum(set_noVarSP),
                noConSP = sum(set_noConSP),
                timeSP = sum(set_timeSP),
                r = r,
                noVar = Int(sum(set_noVarMP) + sum(set_noVarSP)),
                noCon = Int(sum(set_noConMP) + sum(set_noConSP)),
                time = sum(set_timeMP) + sum(set_timeSP)
            )
            print(df)
            break
        end

        ###################################
        U = f_U(γʳ,zᶜ,T,Ω,J)
        I = f_I(γʳ,zᶜ,T,Ω,J)

        f = @variable(m, [H, J], Bin)
        y = @variable(m, [T, Ω, H, J], Bin)

        c_36 = @constraint(m, [h in H, j in J⁰[h] ], f[h,j] == χ[h,j] )
        c_37 = @constraint(m, [t in T, ω in Ω, h in H, i in I[t,ω]], y[t,ω,h,i] ≤ f[h,i] )
        c_38 = @constraint(m, [t in T, ω in Ω, h in H], sum(y[t,ω,h,i] for i in I[t,ω]) ≤ 1 )
        c_39 = @constraint(m, [t in T, ω in Ω, h in H], sum(y[t,ω,h,i] for i in I[t,ω]) ≤ sum(f[h,i] for i in I[t,ω]) )
        c_40 = @constraint(m, [t in T, ω in Ω, h in H], length(J)*sum(y[t,ω,h,i] for i in I[t,ω]) ≥ sum(f[h,i] for i in I[t,ω]) )
        c_42 = @constraint(m, 
        sum( p[ω] * β[t,ω] for t in T for ω in Ω )
        ≥ 
        sum(p[ω]*γʳ[j,t,ω]*zᶜ[j] for t in T for ω in Ω for j in J) +
        sum( (1 - χ[h,j]) * sum( p[ω]*αʳ[h,j,t,ω]*zᶜ[j] for t in T for ω in Ω) for h in H for j in J¹[h] )  - 
        sum(y[t,ω,h,i] * U[t,ω][i] for t in T for ω in Ω for h in H for i in I[t,ω])
        )

        #print("\n c_42 = $c_42 \n")
        #########################################
        # cut I
        υ = @variable(m, [H], Bin )
        δ = @variable(m, [H], Bin )
         con_27_a = @constraint(m, [h in H], sum(χ[h,j] for j in J⁰[h]) ≥ υ[h] )
         con_27_b = @constraint(m, [h in H], sum(χ[h,j] for j in J⁰[h]) ≤ length(J⁰[h]) * υ[h] )
         con_27_c = @constraint(m, sum(δ[h] for h in H) == sum( υ[h] for h in H ) )
         con_27_d = @constraint(m, [l in 2:length(H)], δ[l-1] ≥ δ[l] )    
         con_27_e = @constraint(m, [t in T, ω in Ω], β[t,ω] ≥ θʳ[t,ω] -   
             sum( zʳ[t,ω][l] * δ[l] for l in 1:length(zʳ[t,ω])  ) - 
             θʳ[t,ω] * sum(1- χ[h,j] for j in J for h in H¹[j] )
         )
        #########################################
        # cut II
       # λ = @variable(m, [T,Ω,J] )
       # con_34_a = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω] ], λ[t,ω,j] ≤ γʳ[j,t,ω] * zᶜ[j] )
       # con_34_b = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω] ], λ[t,ω,j] ≤ sum( χ[h,j] * zᶜ[j] for h in H⁰[j] ) )
       # cʳ = [isempty(Jʳ[t, ω]) ? 0 : sum(λ[t, ω, j] for j in Jʳ[t, ω]) for t in T, ω in Ω]

        # con_34_c = @constraint(m, [t in T, ω in Ω], β[t,ω] ≥ θʳ[t,ω] -   
        #    cʳ[t,ω] - 
        #    θʳ[t,ω] * sum(1- χ[h,j] for j in J for h in H¹[j] )
        # )

        #########################################
         # cut III
         ϕ = fill(0.0, length(T), length(Ω), length(J), length(H) )
        # # for updating ϕ
         for t in T
             for ω in Ω
                 for j in Jʳ[t,ω]
                     for h in H⁰[j]
                         if ( zᶜ[j] - sum( αʳ[h,l,t,ω] * zᶜ[l] for l in J) ) > 0
                             ϕ[t,ω,j,h] = zᶜ[j] - sum( αʳ[h,l,t,ω] * zᶜ[l] for l in J)
                         end
                     end
                 end
             end
         end
         ϵ = @variable(m, [T, Ω, J])
         con_29_a = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω]], ϵ[t,ω,j] ≤ sum(χ[h,j] * ϕ[t,ω,j,h] for h in H⁰[j] ) )   
         con_29_b = @constraint(m, [t in T, ω in Ω, j in Jʳ[t,ω]], ϵ[t,ω,j] ≤ γʳ[j,t,ω] * zᶜ[j] )   
         con_29_c = @constraint(m, [t in T, ω in Ω], β[t,ω] ≥ θʳ[t,ω] -
             sum( ϵ[t,ω,j] for j in Jʳ[t,ω] ) - 
             θʳ[t,ω] * sum(1 - χ[h,j] for j in J for h in H¹[j] )
         )
    end
end

initiate_cut531()




