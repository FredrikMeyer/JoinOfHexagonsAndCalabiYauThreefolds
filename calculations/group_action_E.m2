restart

{*
Skal lage gruppevirkning av S_3 på E=k^3 via permutasjonsvirkning.
Da splitter E**E som U^2+U'+W^3, hvor U er den trivielle representasjonen,
U' er sign-rep, og W er standard-representasjonen. Jeg vil ha en basis for W^3.
*}

P = transpose matrix{{0,1,0},{0,0,1},{1,0,0}}
S = transpose matrix{{1,0,0},{0,0,1},{0,1,0}}

-- sjekker at dette er rep av S_3
P^3 -- id
S^2 -- id
S*P*S*P -- id

v = transpose matrix{{1,1,1}} -- basis for trivelle rep

-- basis for W, st.rep
v1 = transpose matrix {{1,-1,0}}
v2 = transpose matrix{{1,0,-1}}

v11 = v1 * transpose v1
v12 = v1 * transpose v2
v21 = v2 * transpose v1
v22 = v2 * transpose v2

{*
Bruker projeksjonsformel til å finne basis for den isotypiske 
komponenten til W.
*}

G =  {P^3,P,P^2,S,S*P, S*P^2}
chiTriv = {1,1,1,1,1,1}
chiST = {2,-1,-1,0,0,0}
chiSign = {1,1,1,-1,-1,-1}

-- standard basis for E
basisE = set apply({matrix{{1,0,0}},matrix{{0,1,0}},matrix{{0,0,1}}},transpose)
-- basis for E tensor E
basisEE = apply(toList(basisE ** basisE), S -> S#0 * transpose S#1)

-- Basis for projection onto isotrivial (?) factors
basisW3 = apply(toList (basisE ** basisE), S -> sum toList apply(0..5, i-> chiST#i * ((G#i * S#0) * (transpose (G#i * S#1))))) -- spenner W^3
basisSign = apply(toList (basisE ** basisE), S -> sum toList apply(0..5, i-> chiSign#i * ((G#i * S#0) * (transpose (G#i * S#1))))) -- spenner sign
basisTrivial = apply(toList (basisE ** basisE), S -> sum toList apply(0..5, i-> chiTriv#i * ((G#i * S#0) * (transpose (G#i * S#1))))) -- spenner U


basisW3 = basisW3_{0..5} --- holder å velge de første seks
 -- spenner triv-delen i W**W:
tt = sum toList apply(0..5, i-> chiTriv#i * ((G#i * v1) * (transpose (G#i * v2))))

Z = matrix{{0,0,0},{0,0,0},{0,0,0}}

b1 = apply(basisW3, M -> (M | Z) || (Z | Z))
b2 = apply(basisW3, M -> (Z | Z) || (Z | M))

basis12 = b1 | b2  --- dette vil gi en join med singlokus av dim > 0. 
--rank matrix toList apply(oo, M -> flatten entries M) 12

w =  {6, 7, 11, 6, 7, 11, 9, 10, 1, 12, 7, 2}
T = ZZ/1009[x_0..x_11, Weights => w] -- ring for P^11--- KAR 1009
T = QQ[x_0..x_11, Weights => w] -- ring for P^11--- KAR 1009
rms1 = apply(0..11, i -> random(T^3,T^3))
rms2 = apply(0..11, i -> random(T^3,T^3))

e1 = (toList basisE)#2
e2 = (toList basisE)#1
e3 = (toList basisE)#0

M_1 = e1 * transpose e2 -- + (e3 * transpose e3)
M_2 = e1 * transpose e3 --+ (e1 * transpose e1)
M_3 = e2 * transpose e3 --+ (e2 * transpose e2)
M_4 = e2 * transpose e1 --+ (e3 * transpose e3)
M_5 = e3 * transpose e1 --+ (e1 * transpose e1)
M_6 = e3 * transpose e2 --+ (e2 * transpose e2)
D_1 = e3 * transpose e3
D_2 = e2 * transpose e2
D_3 = e1 * transpose e1
D_4 = e3 * transpose e3
D_5 = e2 * transpose e2
D_6 = e1 * transpose e1
--blocks = apply(1..6, i-> (M_i | Z) || (Z | Z)) |apply(1..6, i-> (Z | Z) || (Z | M_i))
blocks = apply(1..6, i-> (M_i | Z) || (Z | 2*D_i)) |apply(1..6, i-> (2*D_i | Z) || (Z | M_i)) -- også invariant
--------------
--- når får vi rang 12-rom??
e1 = sub(e1,T)
e2 = sub(e2,T)
e3 = sub(e3,T)

f1 = x_0 *  e2 * transpose e1   + x_1 *e1 * transpose e2 + x_2* e3 * transpose e3
f2 = x_0 *  e3 * transpose e2   + x_1 *e2 * transpose e3 + x_2* e1 * transpose e1
f3 = x_0 *  e1 * transpose e3   + x_1 *e3 * transpose e1 + x_2* e2 * transpose e2
f4 = x_0 *  e1 * transpose e2   + x_1 *e2 * transpose e1 + x_2* e3 * transpose e3
f5 = x_0 *  e3 * transpose e1   + x_1 *e1 * transpose e3 + x_2* e2 * transpose e2
f6 = x_0 *  e2 * transpose e3   + x_1 *e3 * transpose e2 + x_2* e1 * transpose e1
f11 = x_3 *  e2 * transpose e1   + x_4 *e1 * transpose e2 + x_5* e3 * transpose e3
f22 = x_3 *  e3 * transpose e2   + x_4 *e2 * transpose e3 + x_5* e1 * transpose e1
f33 = x_3 *  e1 * transpose e3   + x_4 *e3 * transpose e1 + x_5* e2 * transpose e2
f44 = x_3 *  e1 * transpose e2   + x_4 *e2 * transpose e1 + x_5* e3 * transpose e3
f55 = x_3 *  e3 * transpose e1   + x_4 *e1 * transpose e3 + x_5* e2 * transpose e2
f66 = x_3 *  e2 * transpose e3   + x_4 *e3 * transpose e2 + x_5* e1 * transpose e1
bb  = {f1,f2,f3,f4,f5,f6}
bb2  = {f11,f22,f33,f44,f55,f66}

b1 = apply(0..5, i-> (bb#i | Z) || (Z | bb2#i))

S = (Z | id_(ZZ^3)) || (id_(ZZ^3) | Z)
S * ff1 * S
b2 = apply(b1 , M -> S * M* S)

blocks = toList (b1 | b2)
subl = toList apply(0..5, i-> x_i => random(0,T))
blocks = apply(blocks, b -> sub(b,subl))

----------------

q1 = matrix toList apply(apply(1..6, i-> (M_i | Z) || (Z | Z)) | apply(1..3, i-> (D_i | Z) || (Z | Z)), M -> flatten entries M)
q2 = matrix toList apply(blocks, M -> flatten entries M)
ker transpose (q1 || q2) --- H cap U er 3-dim 

f12 = (M_1 | Z) || (Z | D_1)
f21 = (M_4 | Z) || (Z | D_1)
f12-f21 --- rank 2


b11 = apply(basisW3, M -> M ** R)
b22 = apply(basisW3, M -> M ** R)
vr = sub(random(T^3,T^1),ZZ)
--vr = transpose sub(matrix{{0,1,0}},ZZ)
--vr' = transpose sub(matrix{{0,0,1}},ZZ)
vr' = sub(random(T^3,T^1),ZZ)
vrr =sub(random(T^3,T^1),ZZ)
--vrr' = transpose sub(matrix{{0,0,1}},ZZ)
--vrr = transpose sub(matrix{{1,0,0}},ZZ)
vrr' =sub(random(T^3,T^1),ZZ)
lll = apply(G, g -> (g * vr) *transpose (g * vrr))
lll' = apply(G, g -> (g * vr') *transpose (g * vrr'))
LL = apply(lll, M -> sub(M,T))
LL' = apply(lll', M -> sub(M,T))
ZR = Z ** T
block1 = apply(LL, M -> (M | ZR)|| (ZR | 0))

blocks = apply(0..5, i -> (LL#i | Z) || (Z | LL'#i)) | apply(0..5, i -> (LL'#i | Z) || (Z | LL#i))

--blocks = block1 | block2 
--blocks  = apply(0..5, i -> (b11#i | Z) || (Z | Z)) | apply(0..5, i -> (Z | Z) || (Z | b11#i))
--blocks =  apply(0..5, i-> (LL#i | Z) || (Z | LL'#i)) | apply(0..5, i-> (LL'#i | Z) || (Z | LL#i))
--blocks =  apply(0..5, i-> (b11#i | Z) || (Z | b11#i)) | apply(0..5, i-> (b11#i | Z) || (Z | b11#i))
matrix toList apply(blocks, M -> flatten entries M) ||matrix toList apply(apply(1..6, i-> (M_i | Z) || (Z | Z)) | apply(1..3, i-> (D_i | Z) || (Z | Z)), M -> flatten entries M)

apply(blocks, rank)
--rank matrix toList apply(blocks, M -> flatten entries M)

--randombasis = apply(0..11, j -> sum toList apply(0..11, i-> random(0,R)*blocks#i))
-- sjekk rang 12
--matrix toList apply(randombasis, M -> flatten entries M)
--blocks = apply(0..11, i-> ((random(T^3,T^1)*random(T^1,T^3)| Z) || (Z | random(T^3,T^1)*random(T^1,T^3))))


blocks = apply(blocks, b -> sub(b,T))
AA = sum toList apply(0..11, i-> x_i *blocks#i)
--AAnew = sum toList apply(0..11, i-> x_i *blocksnew#i)
IX =  (minors(2,AA_{0..2}^{0..2}) + minors(2,AA_{3..5}^{3..5}))
IX =  (minors(2,AAnew_{0..2}^{0..2}) + minors(2,AAnew_{3..5}^{3..5}))
AA
decompose radical ideal mingens(IX + ideal(x_6..x_11))--- så den invariante treffer
--- sing lok til M 6*2 ganger. Totalt gir dette 12 del Pezzo sings. 

IX01 = IX + ideal(x_0-1)
IX0 = ideal mingens minimalPresentation(IX01)
ISING = time ideal mingens ideal singularLocus IX0; -- laaang tid: 1200 sek= 20 min !!! (t=1), brukte 283 sek=4 min for t=2.
---- brukte 78 min over Q !!
ISING = time radical ISING ; -- 30 sek, (1.2 sek for t=2), 16 sek over Q

--"singsIX_invariant" << toString ISING << endl << close
--


loadPackage "MinimalPrimes"
installMinprimes()
time decompose ISING -- 13 sings (totalt 12+12= 24 ??) -- 4 sek
---- ser ut som 24+12=36 sings for t=1.
-----  12+12 sings for t=2
LL = decompose ISING

MP = first LL
MP = LL#1
MP
transpose mingens (MP^3 + ISING)
---

lift(jacobian(ISING) ** (ring ISING/MP),ZZ/1009)

use ring IX0
IX0sub = sub(IX0, {x_10 => x_10+261, x_7 => x_7+261, x_5 => x_5-375})
transpose mingens minimalPresentation ideal transpose mingens minimalPresentation sub(IX0sub,(ring MP)/MP^3)
--- opp til orden 2 ser det ut som om origo er et dobbeltpunkt!!
-----
loadPackage "LocalRings"
setMaxIdeal ideal gens ring IX0sub
setMaxIdeal ideal gens ring IX0
ideal transpose localMingens  gens IX0sub
ideal transpose localMingens  gens IX0
W = oo_*  --- her regner jeg med IX0sub !!
f = W#0
(f - (f % x_8))/x_8
g = (f - ((f % x_8)))/x_8
W = apply(W, k -> numerator sub(k, x_8 => -(f % x_8)/g))
f = W#1
(f - (f % x_11))/x_11
g = (f - ((f % x_11)))/x_11
W = apply(W, k -> numerator sub(k, x_11 => -(f % x_11)/g))
f = W#3
(f - (f % x_5))/x_5
g = (f - ((f % x_5)))/x_5
W = apply(W, k -> numerator sub(k, x_5 => -(f % x_5)/g))
Iloka= minimalPresentation(ideal mingens ideal W  + ideal(x_8,x_5,x_11))
ideal transpose localMingens  gens Iloka
Iloka
tangentCone Iloka
decompose oo
M = matrix{{-27*2, 126  , -375, 317},
	   {0    ,-2*411, -1  , 1 },
	   {0    ,    0, -27*2,  0},
	   {0    ,   0 , 0 ,    0}}
M = sub(M + transpose M, ring Iloka)

---- regner i origo
W = oo_*
f = W#3
use ring f
(f - (f % x_6))/x_6
g = (f - ((f % x_6)))/x_6
W = apply(W, k -> numerator sub(k, x_6 => -(f % x_6)/g))
W#3
f =W#0
(f - (f % x_1))/x_1
g = (f - ((f % x_1)))/x_1
W = apply(W, k -> numerator sub(k, x_1 => -(f % x_1)/g))
f = W#1
f = (ideal mingens ideal W )_*#1
(f - (f % x_5))/x_5

IW = ideal transpose mingens minimalPresentation ideal mingens( ideal apply(W,k -> numerator sub(k, x_1 => -(f % x_1)/g)) + x_1)
W = oo_*
f = W#0
use ring f
g = (f - ((f % x_6)))/x_6
IW = ideal transpose mingens minimalPresentation ideal mingens( ideal apply(W,k -> numerator sub(k, x_6 => -(f % x_6)/g)) + x_6)
W = oo_*
f = W#0
use ring f
g = (f - ((f % x_10)))/x_10
IW = ideal transpose mingens minimalPresentation ideal mingens( ideal apply(W,k -> numerator sub(k, x_10 => -(f % x_10)/g)) + x_10)

IW = ideal transpose mingens minimalPresentation ideal mingens( ideal apply(W,k -> numerator sub(k, x_1 => -(f % x_1)/g)) + x_1)
W = IW_*
f = W#0
use ring f
g = (f - ((f % x_6)))/x_6qqq
IW = ideal transpose mingens minimalPresentation ideal mingens( ideal apply(W,k -> numerator sub(k, x_6 => -(f % x_6)/g)) + x_6)


W = IX0_*
use ring IX0
sub(W#2, x_6 => -126*x_7*x_11/(x_7*x_11+x_7*x_8+x_10*x_11+x_8*x_10-252))
W1 = apply(W, f -> numerator sub(f, x_6 => -126*x_7*x_11/(x_7*x_11+x_7*x_8+x_10*x_11+x_8*x_10-252)))
W = (ideal mingens ideal W1)_*
f = W#4
g = (f - ((f % x_1)))/x_1
f == (f % x_1) + g*x_1 -- true
IMM = ideal mingens minimalPresentation ideal mingens (ideal apply(W,k -> numerator sub(k, x_1 => -(f%x_1)/g)) + ideal(x_6,x_1));
W = IMM_*
f = W#5
use ring f
g = (f - (f % x_5))/x_5
f == (f % x_5) + g*x_5 -- true
IMM = ideal mingens minimalPresentation ideal mingens(ideal apply(W, k -> numerator sub(k, x_5 => -(f%x_5)/g)) + ideal(x_5))
W = IMM_*

radical ideal mingens ideal singularLocus ideal W
-----
decompose ideal oo
MP = first oo -- max ideal til en av de potensielle dobbelpktene
MP = ideal(x_1..x_11)
rank sub(lift((jacobian IX) ** (T/MP),ring IX), x_11 => 1) -- 7
rank sub(lift((jacobian IX) ** (T/MP), ring IX), x_0 => 1)-- 7
rank sub(lift((jacobian IX) ** (T/MP),ring IX), x_8 => 1)




sub(jacobian IX01, {x_1 => 0, x_5 => 0, x_6 => 0, x_7 => 0, x_8 => 0, x_10 => 0, x_11 => 0})
time tangentCone IX01


M = mingens ideal mingens ideal singularLocus minimalPresentation(Ip + ideal(x_0-1))

loadPackage "MinimalPrimes"
installMinprimes()

L = decompose ideal M
apply(L, degree)

TT  = ZZ/1009[x_0..x_5]
IX2 = ideal mingens sub(Ip,TT)
radical ideal mingens ideal singularLocus IX2
singIX2 = oo
L = decompose singIX2
apply(L, degree)

toList apply(0..9, j -> toList (j,toList select(0..9, i-> dim (L#j + L#i)==1)))
apply(flatten apply(oo, p -> apply(p#1, i -> {i,p#0})), set)
apply(unique oo, toList)
apply(oo, p -> x_(p#0) * x_(p#1))
loadPackage "SimplicialComplexes"
R = QQ[x_0..x_9]
S = simplicialComplex oo


mingens ideal singularLocus last L
singlastl = ideal oo


 
radical ideal oo

MM = (ideal M)_*
radical ideal MM_{0..20}



---

sub(AA, {x_0 => 0})
IX

pt = sub(AA, {x_0 => 0, x_1 => 0, x_2 => 0, x_3 => 0, x_4 => 0, x_5 => 0, x_6 => 0, x_7 => 0, x_8 => 0, x_9 => 0, x_10 => 0, x_11 => 1})

pt = sub(AA, {x_0 => -374, x_1 => 0, x_2 => 0, x_3 => -374, x_4 => 0, x_5 => 0, x_6 => 0, x_7 => -375, x_8 => 1, x_9 => 0, x_10 => 0, x_11 => 0})

-- ^ dette er et singulært punkt på X
matrix toList apply(blocks, M -> flatten entries M)

sings =unique flatten apply(unique singlist, I -> decompose I)
sings#1

loadPackage "gfanInterface"
gs =  IX_*
mons = {x_0*x_3, x_0*x_4, x_3*x_5,x_1*x_3,x_1*x_4,x_2*x_4,x_0*x_2,x_1*x_5,x_2*x_5,x_6*x_9,x_6*x_10,x_9*x_11,x_7*x_9,x_7*x_10,
    x_8*x_10, x_6*x_8, x_7*x_11, x_8*x_11}
w = weightVector(mons_{0..13}, gs_{0..13})



loadPackage "MinimalPrimes"
installMinprimes()
decompose IX
mingens ideal singularLocus ideal  mingens  minimalPresentation(IX + ideal(x_0-1))
ideal mingens ideal singularLocus ideal (mingens IX)_{0..8}

CT^1( ideal transpose mingens minimalPresentation (IX + ideal(x_0-1)))

PP = oo
decompose radical ideal transpose mingens  minimalPresentation ideal  mingens ideal PP
transpose mingens IX

mingens ideal singularLocus ideal transpose mingens minimalPresentation (IX + ideal(x_0-1))
(dim ideal oo, degree radical ideal oo)


singlist = {}
for i from 0 to 11 do {
    sz = sub(IX, x_i => 1);
    singz = time radical ideal mingens ideal singularLocus  minimalPresentation sz;
    print (dim singz, degree singz);
    singsz = time (decompose singz);
    invz = sz.cache.minimalPresentationMap;
    singlist = singlist | apply(singsz, I -> saturate(homogenize(preimage(invz, singz),x_i)));
    print i;
    print ".";
    }


sings  = unique flatten apply(unique singlist, decompose)
unique flatten apply(singlist, decompose)
sings = oo
#sings --13 * 12 singulariteter!!

Ix0 = minimalPresentation(IX + (x_0-1))
transpose mingens tangentCone(Ix0)

IX + ideal(x_11+x_8, x_7+x_10, x_6+x_9)
decompose oo
apply(oo, dim)
AA

decompose ideal oo
dlist = oo
apply(dlist, degree)
dlist#1

decompose ideal mingens tangentCone(sub(IX, {x_0 => x_0+75*x_10, x_5 => x_5 + 127*x_10-375*x_0, x_7 => x_7+x_10+148*x_0}) + x_10)


---- 48 isolerte! for t=5
sings#29
Iloc = ideal transpose mingens minimalPresentation(IX + ideal (x_11-1))
sub(sings#29, {x_6 => x_6+x_9, x_1 => x_1-x_9, x_2 => x_2-x_9})
Iloc2 = ideal transpose mingens minimalPresentation(sub(IX, {x_6 => x_6+x_9, x_1 => x_1-x_9, x_2 => x_2-x_9})+ideal(x_9-1))
transpose mingens Iloc2


transpose mingens Iloc
AAA = sub(sub(sub(sub(AA, x_11 => 1), x_0 => x_0-x_3), x_2 => x_2-x_5), x_1 => x_1-x_4)
I11 =  (minors(2,AAA_{0..2}^{0..2}) + minors(2,AAA_{3..5}^{3..5}))

G = I11_*
G#2
transpose matrix {apply(G, f -> numerator sub(f, x_4 => x_3*x_5/(x_8+1)))}
GG = (ideal transpose gens minimalPresentation ideal mingens  ideal oo)_*
f = GG#2

ll = (f - (f % x_10))/x_10
GGG = (ideal transpose mingens ideal apply(GG, f -> numerator sub(f, x_10 => -(f % x_10)/ll)))_*
f = GGG#1
ll = numerator ((f - (f % x_6))/x_6)
transpose mingens ideal apply(GGG, f -> numerator sub(f, x_6 =>  -(f% x_6)/ll))
mingens ideal singularLocus ideal oo
decompose radical ideal oo
transpose mingens minimalPresentation ideal transpose mingens I11

loadPackage "VersalDeformations"
CT^1(0, gens IX);

blocks 

torus = ideal apply(flatten apply(apply(apply(flatten entries gens IX, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)
radical ideal mingens torus
loadPackage "Binomials"
toruskomps = BPD torus
toruskomps = select(toruskomps, I -> dim I == 1) --- har en Z_3-virkning!!! (altså den med f_ij)

AA
GG = apply(G, g -> (g | Z) || (Z | g))
AA
GG#1 * AA * (inverse GG#1)
AA
=======
IX = ideal transpose mingens IX
J = ideal (IX_*)_{0..8}
decompose J
JJ = sub(J, {x_10 => x_10 + x_11, x_9 => x_9+x_11, x_7 => x_7+x_8, x_6 => x_6+x_8})
decompose JJ
JJ = sub(JJ, {x_8 => x_8 + x_10})
decompose JJ
JJ = sub(JJ, {x_11 => (x_11-x_8)/374})
decompose J


mingens ideal flatten apply(decompose ideal (IX_*)_{0..8}, m -> flatten entries gens m)
T = ZZ/1009[z_0..z_23]
f = map(R,T, flatten apply(decompose ideal (IX_*)_{0..8}, m -> flatten entries gens m))
ker f
apply(oo, I->preimage(f,I))



loadPackage "MinimalPrimes"
minimalPrimes(IX)



GG = transpose mingens  minimalPresentation(IX + ideal(x_0-1))
GG_0_0
use ring GG
transpose mingens eliminate(ideal sub(GG, x_10 => (1+x_10)/504),x_10)

dim radical( (ideal AA_{0..2}^{0..2})  + IX) --
matrix toList apply(blocks, M -> flatten entries M)
--radical((ideal AA_{0..2}^{0..2}) + minors(2,AA_{3..5}^{3..5})) <-- vi må konkludere med at generelle X ikke snitter Sing(M)

loadPackage "defMethods"
--matrix2Sage transpose oo

apply(oo, I -> (degree I, dim I))

loadPackage "VersalDeformations"

CT^1(0, gens IX);


sub(transpose mingens IX, {x_8 => x_8 + x_10, x_7 => x_7+x_9, x_1 => x_1+x_3})

d1 = b11#0
d2 = b11#5
d1' = (d1-d2) * (1/3)
d2' = d1'+d2
four = b11_{1,2,3,4}

b1' = (four#1+2*four#3)*(1/3)
b2' = four#3 - b1'

b3' = (2*four#0 + four#2)*(1/3)
four
b4' = (b3' - four#2)*(1/2)

B11 = {d1',d2',b1',b2',b3',b4'}

rank matrix apply(B11, M -> flatten entries M) -- fremdeles rang 6


minimalPresentation (ideal oo + ideal(x_8-1))

loadPackage "MinimalPrimes"
installMinprimes()

minimalPrimes sub(ideal ((IX_*)_{0..9}), {x_2 => x_2+x_4, x_1 => x_1+x_4, x_0 => x_0+x_4+x_5})
minimalPrimes ideal ((IX_*)_{0..9})


mingens ideal singularLocus (minimalPresentation ((ideal mingens IX) + (x_5-1)))
transpose mingens IX

CT^1(0, gens IX)


mingens ideal singularLocus minimalPresentation ideal transpose mingens(IX + ideal(x_0..x_5))

P
J = minimalPresentation ideal transpose mingens(IX + ideal(x_0..x_5))
S = ring oo
C = Proj(S/J)
CT^1(0, gens J)
dim singularLocus  J



----------
N = matrix toList apply(blocks, M -> flatten entries M)

rotateMatrix = (M) -> (
    r := rank source M;
    c := rank target M;
    return matrix table(r,c, (i,j) -> M_(c-j-1,r-i-1));
    )

a = gens gb rotateMatrix N
l = rotateMatrix leadTerm a
a = rotateMatrix a


blocksnew = apply(entries a, ff -> matrix toList apply(0..5, i-> ff_{6*i..6*i+5}))




-----------

P
S
PP
I = id_(ZZ^3)
PP = (P | Z) || (Z | P)
SS = (S | Z) || (Z | S)
S2 = (Z | I) || (I | Z)

R = PP^2 * S2

-- da er D_6 generert av R og SS:
SS*R*SS == R^5  -- TRUE
SS^2 == id_(ZZ^6) -- TRUE

D6 = toList (apply(0..5, i-> R^i) | apply(0..5, i-> R^i*SS)) -- som rep på E+E.

---
-- vil ha D6 som rep på E**E + E**E
MMM = (basisW3#0 | Z) || (Z | basisW3#0)
apply(D6, g -> g * MMM * (inverse g))

D6#5

I2 = 2*I

M  = (I2 | Z) || ( Z | I)

M = (random(T^3,T^3) | Z) || (Z | (random(T^3,T^3)))

bane = apply(D6, g -> g * M * (inverse g))
rank matrix toList apply(bane, m -> flatten entries m)

CT^1(0, mingens IX)
--------
--- regne med revlex

mingens ideal singularLocus minimalPresentation (IX + ideal(x_0-1))
sing0 = radical ideal oo
decompose sing0


Iloc = minimalPresentation (IX + ideal(x_0-1))
R = ring Iloc
TT = newRing(R, MonomialOrder => RevLex, Global => false)
transpose mingens gb sub(Iloc,TT)

Iloc : saturate(Iloc)

(decompose sing0)#0
singlocus = radical ideal mingens ideal singularLocus Iloc
(decompose singlocus)#0

describe ring Iloc

tangentCone Iloc
----


---
restart
R = ZZ/1009[x_1..x_9,y_1..y_9]
N = random(R^12,R^18)
a
minors(13,a || genericMatrix(R,1,18))
II = ideal mingens o10

M = (genericMatrix(R,3,3) | Z) || (Z | genericMatrix(R,y_1,3,3))

A = R/II

MA = sub(M,A)
IX = minors(2,MA_{0..2}^{0..2}) + minors(2,MA_{3..5}^{3..5})

mingens ideal singularLocus minimalPresentation (IX + ideal(y_1-1))


---- prove aa se lokalt på sings til IX01 i origo
use ring IX0
M = IX0_*
g = numerator( ( M#2 - (M#2 % x_6))/x_6)
 M#2 % x_ 6 + x_6 * g 
sub(M#2, x_6 => -(M#2 % x_6)/g) == 0
M = apply(M, k -> numerator sub(k, x_6 => -(M#2 % x_6)/g));
M#5
g = numerator( ( M#5 - (M#5 % x_1))/x_1)
sub(M#5, x_1 => -(M#5 % x_1)/g) == 0
M = apply(M, k -> numerator sub(k, x_1 => -(M#5 % x_1)/g));

g = numerator( ( M#8 - (M#8 % x_5))/x_5)
M = apply(M, k -> numerator sub(k, x_5 => -(M#8 % x_5)/g));

J = minimalPresentation (ideal select(M, f -> f!=0)  +ideal(x_5,x_1,x_6))
----- !!! har at radical(J) er dobbelpunkt opp til høyere orden
M = IX0sub_*
M#0
g = numerator( ( M#0 - (M#0 % x_6))/x_6)
sub(M#0, x_6 => -(M#0 % x_6)/g) == 0
M = apply(M, k -> numerator sub(k, x_6 => -(M#0 % x_6)/g));
g = numerator ((M#3 - (M#3 % x_1))/x_1)
sub(M#3, x_1 => -(M#3 % x_1)/g)
M = apply(M, k -> numerator sub(k, x_1 => -(M#3 % x_1)/g));

g = numerator ((M#8- (M#8 % x_5))/x_5)
M = apply(M, k -> numerator sub(k, x_5 => -(M#8 % x_5)/g))

sub(M#7, x_5 => -(M#8 % x_5
	
	)/g)
