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
