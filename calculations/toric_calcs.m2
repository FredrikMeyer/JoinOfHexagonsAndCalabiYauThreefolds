restart
loadPackage "NormalToricVarieties"
PP2 = projectiveSpace 2;

M  = transpose matrix rays (PP2 ** PP2)
P1 =vertices polar convexHull M


Z = 0*random(ZZ^4,ZZ^9)
N = (2*P1 | Z) || (Z | 2*P1) || matrix{{1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1}}
P = convexHull N
polar P
latticePoints P
L = ooo
72-14-33
matrix2Sage N
#latticePoints polar P
#latticePoints P

Sigma = normalFan P
V = normalToricVariety Sigma

loadPackage "defMethods"
matrix2Sage N
--transpose mingens toricIdeal sub(N',ZZ)
N' = N || matrix{{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}}
NN = makeCone toricIdeal sub(N',ZZ)
VV = normalToricVariety normalFan convexHull transpose((transpose NN)_{0..8})
polar polytope sum toList apply(0..11, i-> VV_i)
convexHull transpose((transpose NN)_{0..8})
N
vertices polar polytope sum toList apply(0..11, i-> VV_i)
vertices polar polytope sum toList apply(0..11, i-> V_i)
fan VV

N
isReflexive P



--- try to find crepant resolution
latticePoints P

---
vertices polar P
BLV = blowup({0,1,2,3,4,5},V)
BBLV = blowup({6,7,8,9,10,11},BLV)

fromCDivToPic BBLV
transpose matrix rays BBLV
--decompose dual monomialIdeal BBLV
transpose matrix degrees ring BBLV
ring BBLV
KV = -toricDivisor BBLV
--volume polytope KV



PP17 = projectiveSpace 5
rays oo
BBB = blowup({3,4,5},blowup({0,1,2},PP17))
pic BBB
fan PP17
fan BBB
----
-- for P^1xP^1xP^1 ** P^2 x P^2:
PP1 = projectiveSpace 1;
PP2 = projectiveSpace 2;
P1 = PP1 ** PP1 ** PP1
M1  = transpose matrix rays (P1)
P1 =vertices polar convexHull M1
P2 = PP2 ** PP2
M2  = transpose matrix rays (P2)
P2 =vertices polar convexHull M2
Z1 = 0*random(ZZ^4,ZZ^8)
Z2 = 0*random(ZZ^3,ZZ^9)
N = (2*P1 | Z2) || (Z1 | 2*P2) || matrix{{1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1}}
matrix2Sage N
---
-- for P^1 x P^1 x P^1
PP1 = projectiveSpace 1;
P1 = PP1 ** PP1 ** PP1
M  = transpose matrix rays (P1)
P1 =vertices polar convexHull M

Z = 0*random(ZZ^3,ZZ^8)
N = (2*P1 | Z) || (Z | 2*P1) || matrix{{1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1}}

matrix2Sage N
P = convexHull N
Sigma = normalFan P
V = normalToricVariety Sigma

transpose matrix apply(latticePoints polar P, l -> flatten entries flatten l)
vertices polar P
apply(cones(7, fan V), rays)

VV = blowup({0,1,2,3,4,5},V)
VVV = blowup({6,7,8,9,10,11},VV)
K = -toricDivisor VVV
polytope K == P
polytope fan VVV

isCartier K
-K+VVV_12+VVV_13
isAmple (-(-K+VVV_12+VVV_13))
matrix apply(rays fan VVV,r ->  flatten entries flatten r)
convexHull transpose oo
polar oo

interiorLatticePoints polytope K
apply(faces(1,P),f -> #interiorLatticePoints f)
latticePoints P
numgens ideal VVV
ring VVV
DM = transpose matrix degrees oo
rrr = rays fan VVV
DM_0
numgens ideal VVV
apply(rrr, r -> r*(transpose DM))

M = transpose matrix{{1,0},{1,1},{0,1},{-1,0},{-1,-1},{0,-1}}
P1 =vertices polar convexHull M
tex( (P1 | Z )||(Z | P1))
Z = 0*random(ZZ^2,ZZ^6)
N = (2*P1 | Z) || (Z | 2*P1) || matrix{{1,1,1,1,1,1,-1,-1,-1,-1,-1,-1}}
matrix2Sage N
P = convexHull N
Sigma = normalFan P
V = normalToricVariety Sigma

cl V
pic V
VV = makeSmooth V

cotangentSheaf VV


transpose matrix rays VV
fan V
(cones(5,fan V))#0
rays oo


V = PP2 ** PP2
isFano V
-toricDivisor V
isAmple oo
isAmple (V_0+V_3)
polytope(V_0+V_3)
vertices oo
----
--- finn crepant for polar: dp6xdp6

---

M = transpose matrix{{1,0},{1,1},{0,1},{-1,0},{-1,-1},{0,-1}}
	tex M
dP6P = polar convexHull M
P = dP6P * dP6P
Sigma = faceFan P
Y = normalToricVariety Sigma
--K = -toricDivisor Y
--vertices polytope K
transpose matrix rays Y
vertices polar P

----
restart
loadPackage "defMethods"
R = QQ[x_1..x_9]
minimalPresentation ((minors(2,sub(genericMatrix(R,3,3),x_9 => x_5))) + ideal(x_9))
M = makeCone oo

matrix2Sage transpose ((transpose M)_{0..2})
