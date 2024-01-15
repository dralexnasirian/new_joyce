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

        # cut III

        ϕ = fill(0.0, length(T), length(Ω), length(J), length(H) )
        # for updating ϕ
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
