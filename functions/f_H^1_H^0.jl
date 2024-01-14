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
