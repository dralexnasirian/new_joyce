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

