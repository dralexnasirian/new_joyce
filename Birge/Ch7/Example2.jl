using Distributions
using JuMP, Gurobi

dist = Poisson(4.6)

pdf(dist, 3)

cdf(dist,3)


sum([pdf(dist, i) for i in 0:3])

c = [ 10 12 8 6 5 14;
 8 5 10 15 9 12;
 7 14 4 11 15 8;
 5 8 12 10 10 10]

f = fill(10,4)

r = fill(40,6)

ξ = [4 4 5 3 3 8;
     5 2 4 8 5 6;
     2 8 3 4 7 5;
     3 5 6 4 6 5]

D = [8, 4, 6, 3, 5, 8]

dist = sum(ξ[:,1])

sum([pdf(Poisson(dist),i) for i in 0:7])


dist = 0
cdf(Poisson(dist),8)

2^10
