
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
