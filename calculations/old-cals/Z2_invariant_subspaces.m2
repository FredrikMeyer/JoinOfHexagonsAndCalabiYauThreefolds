restart
Gr = Grassmannian(1,5, CoefficientRing => QQ)
P = ring Gr


I = ideal apply(subsets({0,1,2,3,4,5}, 3), S -> p_(toSequence S) - p_(toSequence sort toList ((set toList {0,1,2,3,4,5}) - set S)))


Inv =  I + Gr

loadPackage "MinimalPrimes"
installMinprimes()
Im = ideal mingens minimalPresentation Inv
transpose gens oo

betti res oo
X = Proj(ring oo/oo)

transpose mingens Im

I = id_(ZZ^3)
Z = 0*I
M = (Z | I) || (I | Z)

(last eigenvectors exteriorPower(2,M))_0

--- 1-egenrommet
apply(select(0..14, i-> (first eigenvectors exteriorPower(2,M))#i == 1), j -> (last eigenvectors exteriorPower(2,M))_j)

apply(select(0..14, i-> (first eigenvectors exteriorPower(2,M))#i == -1), j -> (last eigenvectors exteriorPower(2,M))_j)

eigenvectors M
