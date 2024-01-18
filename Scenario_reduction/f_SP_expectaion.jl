function f_SP_expectaion()
    D = f_D_expectation()
    θ_t = fill(0.0, length(T), 16 )
    γ_t = fill(0.0, length(J), length(T), 16 )
    α_t = fill(0.0, length(H), length(J), length(T), 16 )
    solveTimeSub = fill(0.0, length(T), 16 )
    noVarSub = fill(0.0, length(T), 16 )
    noConSub = fill(0.0, length(T), 16 )
    
    p = [1/16 for i in 1:16]
    set_obj = []
    for ω in 1:16
        for t in T
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
