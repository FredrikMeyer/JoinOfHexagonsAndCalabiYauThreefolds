restart

R = QQ[a,b,c]
S' = QQ[x_0..x_9]
S'' = QQ[x_10..x_19]



f = map(R,S',basis(3,R))
f' = map(R,S'',basis(3,R))

S = QQ[x_0..x_19]
I = sub(ker f,S) + sub(ker f',S)

transpose mingens I
loadPackage "VersalDeformations"

CT^1(gens I)
h1 = sum gens ring J
h2 = sum toList apply(0..19, i-> (-1)^i * x_i)
J = I + ideal(h1,h2)
CT^1(0, gens oo)


ideal mingens ideal mingens(I + ideal (x_4) +ideal( x_14))

betti res oo

