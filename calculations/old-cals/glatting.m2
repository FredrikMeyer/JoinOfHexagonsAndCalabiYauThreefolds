restart
R = ZZ/31[x_1..x_6,z_1..z_6,y_0,y_1,  Weights => {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1}]

I = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)
f = map(R,R,{z_1,z_2,z_3,z_4,z_5,z_6,x_1,x_2,x_3,x_4,x_5,x_6,y_1,-y_0})
J = f I + I

T = R[t_1..t_24]
L = value get "./defdata66";
f1 = L#0; r1 = L#1; g1 = L#2; c1 = L#3;

G1 = ideal mingens ideal g1
loadPackage "MinimalPrimes"
installMinprimes()
Glist = minprimes G1
intersect Glist == G1
apply(Glist, i -> dim i - 14)
transpose mingens Glist#0

use ring sum g1
opts =  {t_1 => 1, t_2 => 3, t_3 => 4, t_4 => 5, t_9 => 6, t_10 => 7,
	t_11 => 2, t_12 => 6, t_5 => 8, t_6 => 10, t_7 => 12, t_8 => 14,
	t_13 => 1,t_14 => 3, t_15 => 4, t_16 => 5, t_21 => 6, t_22 => 7,
	t_23 => 2, t_24 => 6, t_17 => 8, t_18 => 10, t_19 => 12, t_20 => 14
	}
    
    
f1g = sub(transpose sum f1, opts)
Fg = ideal f1g

dim ideal transpose mingens ideal sum g1
28-14



betti res J

------------------- TEST UNDER -------------

k1 = sub(Fg, x_1 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k1) --glatt i kartet

k2 = sub(Fg, x_2 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k2)
dim oo  --svar 2, så sing lok har dim 1

k3 = sub(Fg, x_3 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k3) 
dim oo  --svar 2, så sing lok har dim 1

k4 = sub(Fg, x_4 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k4) 
dim oo --svar 2, så sing lok har dim 1

k5 = sub(Fg, x_5 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k5) 
dim oo ---- samme

k6 = sub(Fg, x_6 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k6) 
dim oo --- samme

ky1 = sub(Fg, y_1 => -1)
transpose mingens ideal minimalPresentation(R/ky1)
----- umulig å regne ut (per nå) (men går bra)

-- samme for y_2

kz1 = sub(Fg, z_1 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation(R/kz1)
dim oo
--- dim 2 - 1 OK

kz2 = sub(Fg, z_2 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation(R/kz2)
dim oo
--- dim 2 - 1 OK

kz3 = sub(Fg, z_3 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation(R/kz3)
dim oo
--- dim 2 - 1 OK

kz4 = sub(Fg, z_4 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation(R/kz4)
dim oo
-- ok

kz5 = sub(Fg, z_5 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation(R/kz5)
dim oo
-- ok

kz6 = sub(Fg, z_6 => -1)
radical ideal mingens ideal singularLocus ideal minimalPresentation(R/kz6)
dim oo
-- ok
----------------- konklusjon: singularitetene i Y glattet har dim 1 overalt


-------------- sjekk struktur f1g --------

f1g
transpose mingens ideal f1g
toString f1g
g1g = matrix {{z_4*z_6+3*x_1*z_5+z_5*y_0},
      {z_3*z_6-30*x_1^2-34*x_1*y_0-8*y_0^2}, {z_2*z_6+5*x_1*z_1+4*z_1*y_0},
      {z_3*z_5+10*x_1*z_4+8*z_4*y_0}, {z_2*z_5-70*x_1^2-116*x_1*y_0-48*y_0^2},
      {z_1*z_5+14*x_1*z_6+12*z_6*y_0}, {z_2*z_4+7*x_1*z_3+6*z_3*y_0},
      {z_1*z_4-42*x_1^2-50*x_1*y_0-12*y_0^2}, {z_1*z_3+6*x_1*z_2+2*z_2*y_0},
      {x_4*x_6+x_5*z_1+3*x_5*y_1}, {x_3*x_6-8*z_1^2-34*z_1*y_1-30*y_1^2},
      {x_2*x_6+4*x_1*z_1+5*x_1*y_1}, {x_3*x_5+8*x_4*z_1+10*x_4*y_1},
      {x_2*x_5-48*z_1^2-116*z_1*y_1-70*y_1^2}, {x_1*x_5+12*x_6*z_1+14*x_6*y_1},
      {x_2*x_4+6*x_3*z_1+7*x_3*y_1}, {x_1*x_4-12*z_1^2-50*z_1*y_1-42*y_1^2},
      {x_1*x_3+2*x_2*z_1+6*x_2*y_1}}
  
  
g1g = matrix {{z_4*z_6+3*y_1*z_5+z_5*y_0},
      {z_3*z_6-30*y_1^2-34*y_1*y_0-8*y_0^2}, {z_2*z_6+5*y_1*z_1+4*z_1*y_0},
      {z_3*z_5+10*y_1*z_4+8*z_4*y_0}, {z_2*z_5-70*y_1^2-116*y_1*y_0-48*y_0^2},
      {z_1*z_5+14*y_1*z_6+12*z_6*y_0}, {z_2*z_4+7*y_1*z_3+6*z_3*y_0},
      {z_1*z_4-42*y_1^2-50*y_1*y_0-12*y_0^2}, {z_1*z_3+6*y_1*z_2+2*z_2*y_0},
      {x_4*x_6+x_5*y_0+3*x_5*y_1}, {x_3*x_6-8*y_0^2-34*y_0*y_1-30*y_1^2},
      {x_2*x_6+4*x_1*y_0+5*x_1*y_1}, {x_3*x_5+8*x_4*y_0+10*x_4*y_1},
      {x_2*x_5-48*y_0^2-116*y_0*y_1-70*y_1^2}, {x_1*x_5+12*x_6*y_0+14*x_6*y_1},
      {x_2*x_4+6*x_3*y_0+7*x_3*y_1}, {x_1*x_4-12*y_0^2-50*y_0*y_1-42*y_1^2},
      {x_1*x_3+2*x_2*y_0+6*x_2*y_1}}  

toString f1g

hilbertPolynomial ideal g1g
hilbertPolynomial ideal f1g


singlist = {}
for i from 1 to 6 do (
    k1 := sub(ideal g1g, x_i => -1);
    k2 := sub(ideal g1g, z_i => -1);
    s1 := radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k1);
    s2 :=  radical ideal mingens ideal singularLocus ideal minimalPresentation (R/k2);
    print dim s1;
    print dim s2;
    singlist = singlist | {sub(s1,R)} | {sub(s2,R)};
    print "..";
    )



minprimes ideal mingens (intersect singlist + ideal g1g)
#oo



loadPackage "VersalDeformations"
CT^2(0, gens ideal f1g)
CT^2(0, gens ideal g1g)

k2 = sub(ideal g1g, y_0 => -1)
transpose gens ideal minimalPresentation (R/k2)


----

I = Fg
Y = Proj(R/I)
J = ideal mingens( Fg + ideal(random(1,R),random(1,R)))


M = presentation (J/J^2)
M;

loadPackage "BGG"

E = QQ[e_1..e_14, SkewCommutative => true]
regularity coker M
tateResolution(M,E,-2,3)

------------
use R
IC = (J + random(1,R) + random(1,R))

I2 =  minimalPresentation(sub(IC, x_2 => -1)+ideal(x_2))
S = ring oo
I2S = sub(I2,S)

ideal singularLocus I2S
ISSS = oo;
radical ideal (gens ISSS)_{0..10}

---
transpose mingens minimalPresentation  sub(Fg, z_1 => 1)
torus = ideal apply(flatten apply(apply(apply(flatten entries oo, monomials), v ->  flatten entries v), j -> subsets(j,2)),
    s -> s_0-s_1)
tt = (decompose torus)#10
minimalPresentation tt
mingens ideal singularLocus ideal oo
transpose gens Fg


----
loadPackage "VersalDeformations"
IX0 = J + random(1,R) + random(1,R)

CT^1(0, gens IX0)
gens gb IX0
IX02 =  minimalPresentation ideal gens gb IX0

f = random(1,R)
g = random(1,R)
A = R/J
res sub(ideal(f,g),A)
betti oo
matrix{{f,g}} ** A

X = Proj A
f = random(1,R)
g = random(1,R)
IX = J + f + g
B = R/IX
XX = Proj 

apply(0..5, i-> HH^i OO_X(-2))
euler OO_X(-2)
CN= sheaf (J/J^2)
euler CN
euler CN(-1)
euler CN(-2)
28-196

use R
J0 = sub(J, {y_0 => 0, y_1 => 0})
IX0 = ideal mingens(J + ideal(y_0,y_1))

B0 = R/IX0
prune ((IX0/(IX0^2))/(J0/(J0^2)))
sheaf ((prune (IX0/J0)) ** B0)
sheaf prune ((IX0^2/J0^2) ** B0)
sheaf(oo ** B0)

X0 = Proj B0
euler (OO_X0(-1))
prune ((IX/J) ** B)
((IX * IX)/(J * J) ** B)
------
loadPackage("defMethods", Reload => true)
loadPackage("NormalToricVarieties", Reload => true)

transpose ((transpose makeCone J)_{0..4})
P = convexHull oo
fVector P
fVector polar  polar P

interiorLatticePoints P
interiorLatticePoints (2*P)
isNormal 
fVector polytope KX
#latticePoints polar polytope KX
#latticePoints polytope KX
V1 = normalToricVariety normalFan polytope KX
V2 = normalToricVariety normalFan polar polytope KX
prune( cl V2 )

VV = toricIdeal transpose (transpose (lift((vertices polar polytope KX | transpose matrix{{0,0,0,0,0}}),ZZ)) | transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1,1}})
VV2 = toricIdeal transpose (transpose (lift((vertices polar polytope KX),ZZ)) | transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1}})
transpose mingens VV2
VV3 = toricIdeal transpose( matrix apply(select(latticePoints polar polytope KX, v -> v != 0), l -> flatten entries l) | transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1,1,1}})
toricIdeal (lift( transpose(transpose (vertices polytope KX) | transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1}}),ZZ) | transpose matrix{{0,0,0,0,1,1},{0,0,0,0,-1,1}})

viewHelp polytope

betti res VV3
betti res J


V = normalToricVariety P
VV = normalToricVariety polytope KX
 --- fan V == fan VV :)
vertices polytope KX
#latticePoints polytope KX
#latticePoints  P
KX = sum toList apply(0..11, i-> V_i)
KXVV = sum toList apply(0..11, i-> VV_i)
isCartier KX
isAmple KX
cl V
dim V
normalFan polytope(KX) == normalFan  P -- true
vertices polytope(KX)
isReflexive polytope(KX) -- samme normalvifte
isReflexive(P)

fVector P
fVector polar P
vertices (2*P)


projEmb(KX)
-first degrees OO KX
viewHelp basis
rays normalFan P
transpose basis(-first degrees OO KX, ring V)
KV = ring V
S =  QQ[t_0..t_13]

KX
v= matrix{{0,0,0,0,0,1}}
(fromCDivToWDiv V) * transpose v
fromCDivToPic V
M = (fromWDivToCl V)
KXW = transpose matrix {{1,1,1,1,1,1,1,1,1,1,1,1}}
M* KXW * matrix{{1/2}}

coeffs = -flatten entries transpose(N*transpose v) --- vår embedding
D = sum toList apply(0..11, i-> coeffs#i * V_i)
isAmple D
isCartier D
isNef D
projEmb D

transpose basis(-first degrees OO D, ring V)
basis(-apply(first degrees OO KX, i-> lift(i/2,ZZ)), ring V)
map(KV,S,basis(-apply(first degrees OO KX, i-> lift(i/2,ZZ)), ring V))
transpose mingens ker oo
vertices polytope KX
latticePoints KX
degrees ring V
fromWDivToCl V


VP = normalToricVariety lift(vertices polar polytope KX,ZZ)

transpose mingens radical IVP
IVP = toricIdeal lift(transpose(transpose vertices polar polytope KX | transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1}}),ZZ)

transpose mingens IVP
hilbertPolynomial(IVP, Projective => false)
hilbertPolynomial(IVP + random(1,ring IVP) + random(1,ring IVP), Projective => false)

hilbertPolynomial(J + random(1,R) + random(1,R), Projective => false)
apply(gens ring IVP,x -> dim ideal mingens ideal singularLocus minimalPresentation (IVP + ideal(x-1)))

(f1,r1,g1,c1) = versalDeformation(gens IVP, CT^1(0, gens IVP), CT^2(0, gens IVP));

gens ring sum f1

IVPgl  = ideal transpose sub(sum f1, apply(gens ring sum f1, t -> (t => -1)))

IVPgl2 = ideal transpose sub(sum f1, toList apply(1..23, i-> t_i => 0) | {t_24 => 1})

transpose mingens IVPgl



ideal transpose mingens minimalPresentation(IVPgl + ideal(i-1))
transpose mingens oo
betti res IVPgl == betti res IVP
euler sheaf module (IVP*IVP)
euler sheaf module (J*J)

findT2(ring IVP,IVPgl)

ggens = IVPgl_*
ff = last ggens
((ff-(last (IVP_*)))/d - d*f)/b

fff = lift((ff +d^2*e-d^2*f)/b, ring ff)
fff
(fff-b*c)/d
(lift((ggens#4-j*k^2+j*k*l)/h, ring fff)+h*l+j*l)/k
ggens#4
((((ggens#4-(IVP_*)#4)+j*k*l)/h)+j*l)/k
(IVP_*)#4
#terms oo
((ggens#5 - (IVP_*)#5)/l  + j*l)/k
(IVP_*)#5

ggens#10

time mingens(minors(2,jacobian IVPgl) + IVPgl)
time mingens(minors(3,jacobian IVPgl) + IVPgl)

hilbertPolynomial(J + random(2,R) + random(2,R), Projective => false)

transpose gens IVP
latticePoints polar polytope KX

transpose gens IVP
IVP1 = ideal (IVP_*)_{2,3,11..17}

IVP1
B = QQ[a,b,c,d,e,f]
IVP1 = sub(IVP1,B)
X = Proj (B/IVP1)
transpose mingens IVP1
CT^1(0, gens IVP1)
(f2,r2,g2,c1) = versalDeformation(gens IVP1, CT^1(0, gens IVP1), CT^2(0, gens IVP1))

transpose sum f2
makeCone IVP1
IVP1
transpose gens IVP1

----------
loadPackage "VersalDeformations"

S = ZZ/31[x_1..x_6,z_1..z_6]--,  Weights => {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4}]
JJ = sub(eliminate(eliminate(J,y_0),y_1),S)

(f1,r1,g1,c1) = versalDeformation(gens JJ, CT^1(0, gens JJ), CT^2(0, gens JJ));

loadPackage "defMethods"
makeCone J
makeCone JJ
loadPackage "NormalToricVarieties"
viewHelp normalToricVariety
transpose ((transpose makeCone J)_{0..4})
transpose ((transpose makeCone JJ)_{0..4})
matrix2Sage transpose ((transpose makeCone JJ)_{0..4})
matrix2Sage transpose ((transpose makeCone J)_{0..4})

P1 = convexHull transpose ((transpose makeCone J)_{0..4})
P2 = convexHull transpose ((transpose makeCone JJ)_{0..4})
V1 = normalToricVariety P1
V2 = normalToricVariety P2

--- sage utregning: JJ og J er samme toriske!!

CT^2(0, gens JJ);
CT^2(0, gens J);



A = S/JJ
C = res J;
N = Hom(image C.dd_1,A)
dI = sub(transpose jacobian C.dd_1,A)
T1 = N/image dI
T1S = prune sheaf T1

HH^0(T1S)

CT^1(0, gens J)
CT^1(0, gens JJ)
 --- merk: T^1_0 ikke det samme som H^0(T^1) alltid... krever Gorenstein eller noe slikt

loadPackage "Depth" 
depth(JJ,ring JJ)
depth(J,ring J)

---------
loadPackage "SimplicialComplexes"
S =  simplicialComplex monomialIdeal (ideal leadTerm J)
SS =  simplicialComplex monomialIdeal (ideal leadTerm J + ideal(x_1*z_1,x_2*z_2))
fVector S
fVector SS


