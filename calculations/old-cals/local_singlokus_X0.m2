restart
R = ZZ/1009[a,b,x_0..x_6]

M = matrix{{x_0,x_1,x_2},{x_4,x_0,x_3},{x_5,x_6,x_0}}
I = minors(2,M)

betti res I

N = {1,a,a^2*b, a^2*b^2,a*b^2,b}
N = {a,a^2*b, a^2*b^2,a*b^2,b}

h1 = sum N + sum toList(x_0..x_6)
h2 = sum toList apply(0..4, i-> random(0,R) * (-1)^i * N#i) + sum toList apply(0..5, i-> random(0,R) * (-1)^i * x_i)

J = I + ideal(h1,h2)


loadPackage "LocalRings"
J0 = (minimalPresentation J)
setMaxIdeal ideal gens ring J0
G = localMingens gens J0

loadPackage "VersalDeformations" 

CT^1(0, G)
