U = f_U(γʳ,zᶜ,T,Ω,J)
        I = f_I(γʳ,zᶜ,T,Ω,J)

        f = @variable(m, [H, J], Bin)
        y = @variable(m, [T, Ω, H, J], Bin)

        c_36 = @constraint(m, [h in H, j in J⁰[h] ], f[h,j] == χ[h,j] )

        c_37 = @constraint(m, [t in T, ω in Ω, h in H, i in I[t,ω]], y[t,ω,h,i] ≤ f[h,i] )
        #for t in T, ω in Ω, h in H, i in I[t,ω]
        #    print("\n c_37[$t,$ω,$h,$i] = $(c_37[t,ω,h,i]) \n")
        #end

        c_38 = @constraint(m, [t in T, ω in Ω, h in H], sum(y[t,ω,h,i] for i in I[t,ω]) ≤ 1 )
        #for t in T, ω in Ω, h in H
        #    print("\n c_38[$t,$ω,$h] = $(c_38[t,ω,h]) \n")
        #end

        c_39 = @constraint(m, [t in T, ω in Ω, h in H], sum(y[t,ω,h,i] for i in I[t,ω]) ≤ sum(f[h,i] for i in I[t,ω]) )
        #for t in T, ω in Ω, h in H
        #    print("\n c_39[$t,$ω,$h] = $(c_39[t,ω,h]) \n")
        #end

        c_40 = @constraint(m, [t in T, ω in Ω, h in H], length(J)*sum(y[t,ω,h,i] for i in I[t,ω]) ≥ sum(f[h,i] for i in I[t,ω]) )
        #for t in T, ω in Ω, h in H
        #    print("\n c_40[$t,$ω,$h] = $(c_40[t,ω,h]) \n")
        #end

        c_42 = @constraint(m, 
        sum( p[ω] * β[t,ω] for t in T for ω in Ω )
        ≥ 
        sum(p[ω]*γʳ[j,t,ω]*zᶜ[j] for t in T for ω in Ω for j in J) +
        sum( (1 - χ[h,j]) * sum( p[ω]*αʳ[h,j,t,ω]*zᶜ[j] for t in T for ω in Ω) for h in H for j in J¹[h] )  - 
        sum(y[t,ω,h,i] * U[t,ω][i] for t in T for ω in Ω for h in H for i in I[t,ω])
        )
