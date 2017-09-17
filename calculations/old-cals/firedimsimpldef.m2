restart
loadPackage "SimplicialComplexes"
loadPackage "VersalDeformations"
kk = ZZ/32003
kk = QQ
w = {1, 1, 1, 4, 7, 4, 1, 1, 1, 4, 7, 4, 2}
T = kk[x_1..x_6,z_1..z_6,y_0,  Weights => {w}]
Ik = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)
f = map(T,T,{z_1,z_2,z_3,z_4,z_5,z_6,x_1,x_2,x_3,x_4,x_5,x_6,y_0})
J
J = Ik + f Ik
h = sum flatten entries vars T
J+h

loadPackage "MinimalPrimes"
installMinprimes()

J + h + random(1,T) + random(1,T) + random(1,T)
degree oo
CT^1(0, gens (J+ h))

singlist = {}

for i from 0 to 11 do {
    sz = sub(J, (gens T)_i => 1);
    singz = radical ideal mingens ideal singularLocus  minimalPresentation sz;
    sings = (decompose singz);
    invz = sz.cache.minimalPresentationMap;
    singlist = singlist| apply(sings, I -> homogenize(preimage(invz, singz),(gens T)_i));
    print i;
    }
loadPackage "SimplicialComplexes"
S = simplicialComplex monomialIdeal intersect singlist

CT^1(0, gens (h+J))



transpose gens gb (J+h)
leadTerm (J+h)
(transpose gens gb (J+h))
loadPackage "gfanInterface"
inL = {x_3*z_1*z_2*z_3,x_1*x_2*x_3*z_1, x_2*x_3*z_1*z_2*z_3}
L = flatten entries (gens gb (J+h))_{4,5,14}
weightVector(inL, L)
simplicialComplex monomialIdeal leadTerm J

h = sum flatten entries vars T
hilbertPolynomial(J + h, Projective => false)


CT^2(0, leadTerm J)
(f1,r1,g1,c1) = versalDeformation(leadTerm J, CT^1(0, leadTerm J), CT^2(0, leadTerm J), HighestOrder => 3, SmartLift => false);


decompose ideal sum g1

-----
J
loadPackage "defMethods"
loadPackage "Polyhedra"

P  = convexHull transpose ((transpose makeCone J)_{0..3})

(C,L,M) = minkSummandCone polar P
apply(values L, vertices)
