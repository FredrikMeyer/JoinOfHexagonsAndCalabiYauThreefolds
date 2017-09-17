restart
----
R = QQ[x_1..x_6,z_1..z_6,y]

M = matrix{{y,x_1,x_2},
    {x_4,y,x_3},
    {x_5,x_6,y}}

M' = matrix{{y,z_1,z_2},
    {z_4,y,z_3},
    {z_5,z_6,y}}

I = minors(2,M) + minors(2,M')

BB = betti res (I)


hilbertPolynomial(I, Projective => false)
I + ideal(x_1,x_2,x_3,x_4,x_5,x_6,y)
--I + ideal(x_1-1,x_2-1,x_3-1,x_4-1,x_5-1,x_6-1,y-1)
mingens oo
decompose ideal mingens oo
apply(oo, dim)

loadPackage( "defMethods" , Reload => true)

N = makeCone I
N= transpose ((transpose N)_{0..3})
loadPackage "Polyhedra"
loadPackage "NormalToricVarieties"
C = posHull N
FanY = normalFan convexHull N

TV = normalToricVariety(FanY)

faceFan convexHull N
TVp = normalToricVariety(oo)
pic oo

C0 = first cones(3, fan TV)
dualCone C0
rays fan TV
rays fan TVp
TVsm = makeSmooth TV
fVector first cones(4,fan TV)

P = convexHull N
latticePoints polar P
rays fan TV

faceswithinterior = select(faces(2,polar P), f -> #interiorLatticePoints f == 1)


listofpos = apply(faceswithinterior,
    f -> positions(rays TV,
	v -> any(entries transpose vertices  f,
	    w -> w == v)))

TVbl = TV
apply(listofpos, pos -> TVbl = blowup(pos,TVbl))

pic TVbl

singcones = select(cones(4,fan TVbl), C -> not isSmooth C)
isFano TV



pfirst faces(2,polar P)

Mv = vertices oo
(Mv * transpose  matrix{{1,1,1,1,1,1}}) * (1/6)

(rays TV)_{0,6,12,18,24,


TVbl = blowup({0,6,12,18,24,30},TV)
fan TVbl
apply(cones(3,fan TVbl), isSmooth)

rays TVbl
rays TV

TVbl2 = blowup({1,7,13,19,25,31},TVbl)
TVbl3 = blowup({2,8,14,20,26,32},TVbl2)
fan oo
fan TVbl2
fan TVbl
isProjective TVbl2
hilbertBasis dualCone C
hilbertBasis  C
rays dualCone C
rays C
isSmooth C
isSmooth dualCone C


-- This is 4-dim toric with 1-dim sings.
-- 48 P^1s.
-- Intersecting as follows:
-- Hexagon (all lines between points in the two hexagons) hexagon

betti res I
betti res (I*I)
S = QQ[x_1..x_6,z_1..z_6]

M = matrix{{z_1+x_4,x_1,x_2},
    {x_4,z_1+z_3,x_3},
    {x_5,x_6,z_1+z_5}}

M' = matrix{{x_1+x_2+x_3+x_4,z_1,z_2},
    {z_4,x_1+x_3,z_3},
    {z_5,z_6,x_1+x_5+x_6}}


IX = minors(2,M) + minors(2,M')

torus = ideal apply(flatten apply(apply(apply(flatten entries gens IX, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)
primaryDecomposition torus
select(oo, I -> dim I == 1)


IT = mingens (I + ideal(x_1-x_3, x_3-x_5) + ideal(z_1-z_3, z_3-z_5))
ideal gens R

primaryDecomposition saturate(ideal IT, ideal gens R)
primdecomp = primaryDecomposition ((ideal IT)+ ideal(x_2-1))


select(primdecomp, I -> radical I == I)


IX= I + (x_1+x_2+x_3+x_4+x_5+x_6+z_1+z_2+z_3+z_4+z_5+z_6+y) + ideal(x_1-1)

IXm = minimalPresentation IX
torus = ideal apply(flatten apply(apply(apply(flatten entries gens IXm, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)
------
------

restart
kk = ZZ/1009
R = kk[x_1..x_6,z_1..z_6,y_0,y_1,y_2,y_3,y_4,y_5]
M = matrix{{y_0,x_1,x_2},
    {x_4,y_1,x_3},
    {x_5,x_6,y_2}}

M' = matrix{{y_3,z_1,z_2},
    {z_4,y_4,z_3},
    {z_5,z_6,y_5}}

I = minors(2,M) + minors(2,M') --- <- RIGID!!!! Har T^1 og T^2 lik null
IX = I + random(1,R) + random(1,R) + random(1,R) + random(1,R) + random(1,R) + random(1,R)
IXP = eliminate(toList(z_1..z_6) | {y_3,y_4,y_5},IX) + ideal(z_1..z_6,y_3,y_4,y_5)
--torus = ideal apply(flatten apply(apply(apply(flatten entries gens oo, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)
--radical torus
XP = Proj(R/IXP)
IXPM =  minimalPresentation IXP;
XPM = Proj(ring IXPM/IXPM)
HH^3(OO_XPM) -- 100
betti res IXPM
--PP = Proj ring IXPM
--HH^0(OO_PP(3))
--HH^0(OO_XPM(1)), 100, så antakelig er w_X = O_X(3)
IM = ideal mingens sub(I,ring IXPM) 
IM + IXPM == IXPM
S = kk[a_1,a_2,a_3,b_1,b_2,b_3]
phi = map(S,ring IXPM,{a_1*b_2,a_1*b_3,a_2*b_3,a_2*b_1,a_3*b_1,a_3*b_2,a_1*b_1,a_2*b_2, a_3*b_3})
degree phi IXPM
 mingens phi IXPM --- X er O(6,6)-snitt
euler  XPM
f = (gens minimalPresentation (IXPM + ideal(x_1-1)))_0_0;
ideal mingens IXPM
loadPackage "VersalDeformations"
prune (IXPM/IXPM^2);
CN = sheaf oo;
A = ring IXPM/IXPM;
NS = sheaf Hom(prune (IXPM/IXPM^2),A)


-- I <--  har sing.lok av dim 4. Vi kutter med 6 polynomer, så svaret er glatt!
7*14
-----
restart
w =  {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4}
R = ZZ/1009[x_0..x_11]--, Weights => w]
rms1 = apply(0..11, i -> random(R^3,R^3))
rms2 = apply(0..11, i -> random(R^3,R^3))

rms1 = apply(0..11, i-> matrix {{0,random(0,R),random(0,R)},{random(0,R),0,random(0,R)},{random(0,R),random(0,R),0}})
rms2 = apply(0..11, i-> matrix {{0,random(0,R),random(0,R)},{random(0,R),0,random(0,R)},{random(0,R),random(0,R),0}})

A = sum toList apply(0..11, i-> x_i * random(0,R) * sub(matrix{{1,0,0},{0,1,0},{0,0,1}},R))
B = sum toList apply(0..11, i-> x_i * random(0,R) * sub(matrix{{1,0,0},{0,1,0},{0,0,1}},R))


---- utspenner generisk P^11 i P^17= P(E o E) * P(E o E)
A = sum toList apply(0..11, i -> x_i * rms1#i)
B = sum toList apply(0..11, i -> x_i * rms2#i)

IX = minors(2,A) + minors(2,B)


A1 = matrix{{0,1,0},{0,0,0},{0,0,0}}
A2 = matrix{{0,0,1},{0,0,0},{0,0,0}}
A3 = matrix{{0,0,0},{0,0,1},{0,0,0}}
A4 = matrix{{0,0,0},{1,0,0},{0,0,0}}
A5 = matrix{{0,0,0},{0,0,0},{1,0,0}}
A6 = matrix{{0,0,0},{0,0,0},{0,1,0}}
D = matrix{{1,0,0},{0,1,0},{0,0,1}}
L1 = apply({A1,A2,A3,A4,A5,A6,0,0,D,D,D,0}, a -> sub(a,R))
L2 = apply({0,D,0,D,D,D,A1,A2,A3,A4,A5,A6}, a -> sub(a,R))

A = sum toList apply(0..11, i -> x_i * L1#i)
B = sum toList apply(0..11, i -> x_i * L2#i)

IX = (minors(2,A) + minors(2,B))
loadPackage "MinimalPrimes"
installMinprimes()

transpose mingens IX
mingens ideal singularLocus minimalPresentation (IX + (x_11-1))
decompose ideal oo



M1 = matrix{{x_1+x_2,x_1,x_2},{x_4,x_1+x_2,x_3},{x_5,x_6,x_1+x_2}}
M2 = matrix{{x_0+x_8,x_7,x_8},{x_10,x_7+x_8,x_9},{x_11,x_0,x_7+x_8}}

minors(2,M1) + minors(2,M2)
hilbertPolynomial oo

--------
restart
R = QQ[x_1..x_18]
M = matrix{{x_1,x_2,x_3,0,0,0},
    {x_4,x_5,x_6,0,0,0},
    {x_7,x_8,x_9,0,0,0},
    {0,0,0,x_10,x_11,x_12},
    {0,0,0,x_13,x_14,x_15},
    {0,0,0,x_16,x_17,x_18}}


ideal mingens minors(3,M)
apply(1..5, i-> degree minors(i,M))

saturate(minors(3,M), ideal(x_1..x_9)*ideal(x_10..x_18))
transpose mingens oo

loadPackage "MinimalPrimes"
installMinprimes()
decompose minors(3,M)

decompose minors(5,M)
apply(oo, dim)

---
--SR
restart
R = QQ[x_1..x_6,z_1..z_6]

M1 = matrix{{0,x_1,x_2},{x_4,0,x_3},{x_5,x_6,0}}
I1 = minors(2,M1)
M2 = matrix{{0,z_1,z_2},{z_4,0,z_3},{z_5,z_6,0}}
I2 = minors(2,M2)

betti res (I1+I2)
135^2

restart
R = QQ[x_1..x_9]
IS = minors(2,genericMatrix(R,3,3))
betti res(sub(I1,R))
betti res IS

betti res IS^2
betti ((res IS^2) ** (res IS^2))
BBB = oo
genericMatrix(R,3,3)
C = res IS
C.dd_2
M = genericMatrix(R,3,3)
gens IS
loadPackage "VersalDeformations"

betti res IS

P = sub(transpose matrix{{0,1,0},{0,0,1},{1,0,0}},R)
S = sub(transpose matrix{{1,0,0},{0,0,1},{0,1,0}},R)
G =  {P^3,P,P^2,S,S*P, S*P^2}

apply(G, g -> g*M*(inverse g))


betti res IS
transpose C.dd_2
ker C.dd_3
ker transpose C.dd_2
betti(C ** C)
BB
C =  res (IS)
betti(C**C)
BB
CtC = C ** C
betti C
betti (C**C)
16*9
viewHelp betti
betti(C ** C)
C

restart
R = QQ[x_1..x_9,z_1..z_9]
M1 = matrix{{x_1,x_2,x_3},{x_4,x_5,x_6},{x_7,x_8,x_9}}
M2 = matrix{{z_1,z_2,z_3},{z_4,z_5,z_6},{z_7,z_8,z_9}}
I = minors(2,M1) + minors(2,M2)

betti res I
betti res I^2

HH^9(sheaf prune (I/I^2))

IX = I + random(1,R)

betti res prune (IX^2/I^2)

restart
loadPackage "NormalToricVarieties"
M = transpose matrix{{1,0},{0,1},{-1,1},{-1,-1},{1,-1}
M = transpose matrix{{1,0},{0,1},{-1,1},{-1,-1},{1,-1},{0,0}}
M = transpose matrix{{1,0},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1},{0,0}}	
P = convexHull M
F = normalFan P
V = normalToricVariety F

loadPackage "defMethods"


toricIdeal(M || matrix{{1,1,1,1,1,1}})
toricIdeal(M || matrix{{1,1,1,1,1,1,1,1}})
betti res oo

