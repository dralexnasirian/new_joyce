function f_D_expectation()
    D = fill(0, length(J), length(T), 16)
    l = []
    for a in 0:1
        for b in 0:1
            for c in 0:1
                for d in 0:1
                    push!(l, [a, b, c, d])
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
