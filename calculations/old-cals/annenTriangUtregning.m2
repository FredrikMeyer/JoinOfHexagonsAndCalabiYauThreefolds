restart
loadPackage "SimplicialComplexes"
loadPackage "VersalDeformations"
kk = ZZ/32003
kk = QQ
w =  {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1}
w2 =  {7, 7, 8, 9, 8, 7, 1, 5, 14, 19, 14, 5, 5, 7}
w3 = {4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 1, 1}
T = kk[x_1..x_6,z_1..z_6,y_0,y_1,  Weights => w2]
Ik = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)
f = map(T,T,{z_1,z_2,z_3,z_4,z_5,z_6,x_1,x_2,x_3,x_4,x_5,x_6,y_1,-y_0})
J
J = Ik + f Ik
--- J er idealet til dP6 * dP6 (join)

sigma = map(T,T,{x_2,x_1,x_6,x_5,x_4,x_3,z_1,z_2,z_3,z_4,z_5,z_6,y_0,y_1})
tau   = map(T,T,{x_5,x_4,x_3,x_2,x_1,x_6,z_1,z_2,z_3,z_4,z_5,z_6,y_0,y_1})
alfa  = map(T,T,{z_1,z_2,z_3,z_4,z_5,z_6,x_1,x_2,x_3,x_4,x_5,x_6,y_1,y_0})

sigma J == J
tau J == J
alfa J == J
--- Abelsk gruppevirkning av orden 4x4x2=32.


h1 = x_1+x_2+x_3+x_4+x_5+x_6+z_1+z_2+z_3+z_4+z_5+z_6
h2 = x_1-x_2+x_3-x_4+x_5-x_6-z_1+z_2-z_3+z_4-z_5+z_6

sub(( sigma h2), {x_1 => 0, x_2 => 0, x_3 => 0, x_4 => 0, x_5 => 0, x_6 => 0})

sigmatau = map(T,T, sigma tau vars T)
iden = map(T,T, vars T)
z2z2_1 = {sigma, tau, sigmatau, iden}
z2z2_2 = {alfa * sigma * alfa , alfa * tau * alfa , alfa * sigmatau * alfa, iden}
z2     = {iden ,alfa}
--- product of the two grups
Z2Z2_2 =  apply(toList( (set z2z2_1) ** (set z2z2_2)), S -> S#0 * S#1)
G = apply(toList( (set Z2Z2_2) ** (set z2)), S -> S#0 * S#1)

--- make random  invarian deg 1
h1' = random(1,T)
h2' = random(1,T)
h1G = (sum apply(G, g -> g(h1')))*45*7
h2G = (sum apply(G, g -> g(h2')))*63*5
----
c1 = {1,1,2,1,1,2,1}  -- trennger alle coeffs å være forskjellige?
c2 = {2,2,1,2,2,1,1}
h1G = sum toList apply(0..5, i-> c1#i * (gens T)#i) + f sum toList apply(0..5, i-> c1#i * (gens T)#i) + y_0 + y_1
h2G = sum toList apply(0..5, i-> c2#i * (gens T)#i) + f sum toList apply(0..5, i-> c2#i * (gens T)#i) + y_0 + y_1

IX = J + h1G + h2G
sigma IX == IX
tau IX == IX
alfa IX == IX

fix1  = ideal(vars T - sigma vars T)
fix2 = ideal(vars T - tau vars T)
fix3 = ideal(vars T- alfa vars T)
transpose gens ideal mingens (fix1 + fix2 + fix3)
transpose gens ideal mingens (fix1 + fix2)
transpose gens ideal mingens (fix1 + fix3)
mingens(ideal oo + IX)
radical ideal oo

time mingens ideal singularLocus minimalPresentation(IX + ideal(x_1-1)) -- brukte 3.5 timer
SING = oo;
time dim ideal oo --- dim 0, 8 sek
radical ideal SING
--
--  Conclusion: exists a degeneration of X_i to something with many singularities, but with
-- many discrete symmetries as well


---- finne sing til IX (i én del PEzzo)
R = QQ[x_1..x_6,y_0]-- ,r]/(r^2+35)
Ik = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)
--h1 = x_1+x_2+x_3+x_4+x_5+x_6
--h2 = x_1-x_2+x_3-x_4+x_5-x_6
h1 = sub(h1G, R)
h1 = x_1+x_2+x_3+x_4+x_5+x_6+y_0
h2 = sub(h2G, R)

dim (Ik + h1 + h2)

IX1 = Ik + h1 + h2
L = decompose IX1
IX11 =  ideal mingens(last L + (y_0-1))
SS  = QQ[a,b,c,d,e,ff]
fss = map(SS,R,{a,b,c,d,e,ff,1})
mingens fss IX11
factor 84
factor 72

mI = ideal transpose gens ideal mingens minimalPresentation last L

eliminate(eliminate(mI, x_6), x_5)
mII = minimalPresentation ideal mingens(mI + (x_5-1))
use ring oo
eliminate(mII, x_4)

I11 = L#0
transpose mingens sub(L#1, y_0 => 1)

sigma1 = map(R,R,{x_2,x_1,x_6,x_5,x_4,x_3,y_0})
tau1   = map(R,R,{x_5,x_4,x_3,x_2,x_1,x_6,y_0})



L2y = ideal mingens (L#2 + ideal(y_0-1))
minimalPresentation L2y

sigma1 L#0 == L#5
tau1 L#1 == L#5
tau1 L#3 == L#2
sigma1 L#3 == L#3

JJ = minimalPresentation(J + ideal(y_0-y_1))

use ring JJ
--dim radical ideal mingens ideal singularLocus ideal mingens minimalPresentation(JJ + ideal(y_0-1))

singlist = {}
for i from 1 to 6 do {
    sz = sub(JJ, z_i => 1);
    sx = sub(JJ, x_i => 1);
    singz = radical ideal mingens ideal singularLocus  minimalPresentation sz;
    singx = radical ideal mingens ideal singularLocus  minimalPresentation sx;
    singsx = (decompose singx);
    singsz = (decompose singz);
    sings = singsx | singsz;
    invz = sz.cache.minimalPresentationMap;
    invx = sx.cache.minimalPresentationMap;
    singlist = singlist | apply(singsx, I -> homogenize(preimage(invx, singx),x_i)) | apply(singsz, I -> homogenize(preimage(invz, singz),z_i));
    print i;
    }

loadPackage "SimplicialComplexes"
simplicialComplex monomialIdeal intersect singlist
sum apply(decompose radical intersect singlist, degree)
R = ring intersect singlist
set vars R ** set {x_3,x_4}
S1 = set {x_1,x_2,x_3,x_4,x_5,x_6}
S2 = set {z_1,z_2,z_3,z_4,z_5,z_6}
prune HH ooo
S = flatten entries vars R
h2 = sum apply(S ** S, s -> s#0 * s#1)
h2 = sum apply(toList (S1 ** S2), s -> s#0 * s#1)
h1 = (x_1+x_2+x_3+x_4+x_5+x_6+z_1+z_2+z_3+z_4+z_5+z_6)
singsCY = ideal mingens ((intersect singlist) + ideal(x_1+x_2+x_3+x_4+x_5+x_6+z_1+z_2+z_3+z_4+z_5+z_6))
apply(decompose singsCY, I -> dim (I+h2))

surf =   (JJ + h2 + h1)
mingens ideal singularLocus minimalPresentation(surf + (x_1-1))

JJ1 = JJ + ideal(x_1-1)
sing1 = first decompose radical ideal mingens ideal singularLocus minimalPresentation(JJ1)

h = JJ1.cache.minimalPresentationMap (x_1+x_2+x_3+x_4+x_5+x_6+z_1+z_2+z_3+z_4+z_5+z_6)

jsnitt = JJ1+sub(h,ring JJ1)
JJmin = minimalPresentation(jsnitt)

transpose mingens JJmin

mingens(JJmin + jsnitt.cache.minimalPresentationMap sub(sing1, ring jsnitt))
use ring JJmin
JJtrans = ideal transpose mingens sub(JJmin, z_6 => z_6-1)

transpose mingens JJtrans
minimalPresentation tangentCone JJtrans --- kan tyde på at JJ har singulariteter som er C(P1xP1) ? 

loadPackage "LocalRings"
setMaxIdeal ideal gens ring JJtrans
localMingens  mingens JJtrans
gg = (ideal oo)_*
loadPackage "MinimalPrimes"
installMinprimes()

apply(decompose ideal gg, I -> mingens(I + ideal gens ring I))
(decompose ideal gg)#1
tangentCone oo

(JJtrans_*)#0

minimalPresentation ideal apply(flatten entries transpose gens sub(JJtrans, {z_3 => z_3/(z_6-1), z_2 => z_2/(z_6-1)}), f -> numerator f)
minimalPresentation oo
transpose mingens oo
ggg = ideal oo
last ggg_*
ggg
decompose radical ideal mingens ideal singularLocus ggg
tangentCone ggg
------
h = x_1+x_2+x_3+x_4+x_5+x_6+z_1+z_2+z_3+z_4+z_5+z_6
CT^2(0,(JJ + h));



CT^1(0, gens (JJ + h1))
--CT^2(0, gens (JJ + h1))

CT^2(0, gens JJ)
(f1,r1,g1,c1) = versalDeformation(gens JJ, CT^1(0, gens JJ), CT^2(0, gens JJ), HighestOrder => 2, SmartLift => false);
sum g1
decompose ideal sum g1


---
use ring J
minimalPresentation (J + ideal(y_0-y_1))


