H = 36
m = Model(Gurobi.Optimizer)
@variable(m , χ[1:H, 1:J] , Bin)
@variable(m , α[1:H , 1:J , 1:T , 1:Ω], Bin);
@variable(m , 0 ≤ γ[1:J , 1:T , 1:Ω], Int )


@objective(m , Min , 
    sum(T * zᵖ[h] * s[h,j] * χ[h,j] for j in 1:J for h in 1:H ) + 
    sum(zᵐ[h,j] * χ[h,j] for h in 1:H for j in 1:J ) +
    sum( p[ω] * zᶜ[j] * γ[j,t,ω]  for j in 1:J for t in 1:T for ω in 1:Ω )
)
# 24 existing workers constraint
existing_workers = 24 
@constraint(m, [h in 1:existing_workers], sum(χ[h,j] * s[h,j] for j in 1:J) == 1 )  

c1 = @constraint(m, [h in 1:H , j in 1:J], sum(s[h,l] * χ[h,l] for l in 1:J) ≥ χ[h,j]  );
c2 = @constraint(m, [h in 1:H , j in 1:J , t in 1:T , ω in 1:Ω], α[h,j,t,ω] ≤ χ[h,j] )
c3 = @constraint(m, [h in 1:H , t in 1:T , ω in 1:Ω], sum( α[h,j,t,ω] for j in 1:J) ≤ 1 )
c4 = @constraint(m, [j in 1:J , t in 1:T , ω in 1:Ω], sum( α[h,j,t,ω] for h in 1:H) + γ[j,t,ω] == D[j,t,ω]);


            # extra constraints
            # c5 = @constraint(m, [h in 1:H], sum(χ[h,j] for j in 1:J) ≤ flexibility )

optimize!(m)

χᵏ = value.(χ)
αᵏ = value.(α)
sum(αᵏ[h,j,t,ω] * p[ω] for h in 1:H for j in 1:J for t in 1:T for ω in 1:Ω ) / (T * H)
γᵏ = value.(γ)

No_P = sum(s[h,j] * χᵏ[h,j] for j in 1:J for h in 1:H ) 
No_S = sum(χᵏ)
No_M = No_S - No_P
No_C = sum( p[ω] * γᵏ[j,t,ω]  for j in 1:J for t in 1:T for ω in 1:Ω ) / T
cost_per = sum(T * zᵖ[h] * s[h,j] * χᵏ[h,j] for j in 1:J for h in 1:H ) 
cost_multi = sum(zᵐ[h,j] * χᵏ[h,j] for h in 1:H for j in 1:J ) 
cost_cas = sum( p[ω] * zᶜ[j] * γᵏ[j,t,ω]  for j in 1:J for t in 1:T for ω in 1:Ω)

no_Cas = sum( p[ω] * γᵏ[j,t,ω]  for j in 1:J for t in 1:T for ω in 1:Ω )
utilization = sum(αᵏ[h,j,t,ω] * p[ω] for h in 1:H for j in 1:J for t in 1:T for ω in 1:Ω ) / (T * H)
dt = DataFrame( 
    T= T,  
    No_P = No_P, 
    No_S = No_S,
    No_M = No_M,
    No_C = No_C,
    cost_per = cost_per,
    cost_multi = cost_multi,
    cost_cas = cost_cas,
    cost_tot = objective_value(m),
    utilization = utilization
    )

    