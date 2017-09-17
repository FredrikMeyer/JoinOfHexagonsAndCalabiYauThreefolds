restart
R = QQ[x_0..x_5]
M = matrix{{0, x_0, 0, x_1, 0,x_2},
          {0, 0 ,x_2, 0, x_3, x_3},
          {0, 0, 0, x_4, 0, x_5},
               {0,0,0,0, 0 ,  x_4},
                   {0,0,0,0,0,x_0},
    {0,0,0,0,0,0}}
M = M   - transpose M
I = pfaffians(6,M)
I = ideal sum apply(gens R, f -> f^3) --- fermat cubic
I = ideal sum toList apply(0..5,  i -> x_(i % 6)*x_((i+1) % 6)*x_((i+2) % 6))
I = ideal(x_0^3 + x_1^3+x_2^3 + x_3*x_4*x_5)
I = ideal(x_0*x_1*x_2-x_3*x_4*x_5)
X = Fano(1,I)

XX = Proj((ring X)/X)
apply(0..4, i-> HH^0(OO_XX))
    
transpose mingens X
degree X
hilbertPolynomial X

loadPackage "MinimalPrimes"
installMinprimes()

L = decompose X
#oo
apply(L, degree)
transpose mingens minimalPresentation last L
loadPackage "defMethods"

makeCone minimalPresentation last L
--hilbertPolynomial (minimalPresentation last L, Projective => false)
J = minimalPresentation last L
use ring J
minimalPresentation(J + ideal(p_6-1))
mingens ideal singularLocus oo

M = mutableMatrix makeCone J
rowAdd(M,4,1,0)
rowAdd(M,4,1,3)
rowAdd(M,4,1,2)
M = matrix M
N = transpose((transpose M)_{0,1,2,3})
loadPackage "NormalToricVarieties"
P =  convexHull N
-- interiorLatticePoints (3*oo)
normalFan P
V = normalToricVariety oo
isProjective V

interiorLatticePoints (3*P)
P2 = affineImage(3*P, transpose matrix{{-1,-1,-1,-1}})
select(cones(4,fan V), c -> isSmooth c  == false) --- seks singulære punkter
first oo
fan V
toricIdeal rays oo
W = makeSmooth V



use ring X
transpose mingens minimalPresentation(X + p_8)
decompose ideal oo

