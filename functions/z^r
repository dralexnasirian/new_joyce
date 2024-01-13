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
