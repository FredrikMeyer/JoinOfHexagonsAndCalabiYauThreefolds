restart
loadPackage "MinimalPrimes"
loadPackage "VersalDeformations"
loadPackage "SimplicialComplexes"
installMinprimes()
kk = ZZ/101
R = kk[x_1..x_6,z_1..z_6,y_0,y_1,  Weights => {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1}]
T = tensor(R,kk[t_1..t_24], DegreeRank => 1, Degrees => {38:1})
(f1,r1,g1,c1) =  value get "./nyGlatting";

transpose mingens ideal sum  g1
comps = decompose ideal mingens ideal g1;
apply(comps, i -> dim i - dim R)

G0 = comps#0
G1 = last comps
loadPackage "defMethods"


opts0 = {t_17 => 1, t_18 => 2, t_20 => 3, t_19 => 4, t_23 => 5,
    t_15 => 1, t_16 => 2, t_22 => 3, t_21 => 4, t_13 => 5,
    t_2 => 1, t_1 => 2, t_10 => 3, t_9 => 4, t_4 => 5, t_3 => 6,
    t_12 => 1, t_11 => 2, t_8 => 3, t_7 => 4, t_6 => 5, t_5 => 6,
    t_14 => 6, t_24 => 6  
    }
opts = {t_17 => 1, t_18 => 2*t_1, t_20 => 3*t_1, t_19 => 4, t_23 => 5,
    t_15 => 1, t_16 => 2*t_1, t_22 => 3*t_1, t_21 => 4, t_13 => 5,
    t_2 => 1*t_1, t_1 => 2, t_10 => 3*t_1, t_9 => 4, t_4 => 5*t_1, t_3 => 6,
    t_12 => 1*t_1, t_11 => 2, t_8 => 3*t_1, t_7 => 4, t_6 => 5*t_1, t_5 => 6,
    t_14 => 6*t_1, t_24 => 6*t_1  
    }

transpose sum f1
apply(comps, I -> sub(I,opts) == 0)
F1 = sub(transpose sum f1, opts)
f1g = ideal F1
f1g = ideal sub(transpose gens f1g, R)
sub(transpose gens f1g, t_1 => 0)
Y = ideal oo
singlist = {}
for i from 1 to 6 do {
    sz = sub(f1g, z_i => 1);
    sx = sub(f1g, x_i => 1);
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

simplicialComplex monomialIdeal intersect singlist -- <- sing.lokus er to kjeder av 4 P^1.
prune HH oo

ISS = intersect singlist;


r1 = sub(random(1,R),T)
r2 = sub(random(1,R),T)

use ring f1g
f1g1 =  sub(f1g, x_3 => 1)
minf1g1 = minimalPresentation f1g1

transpose gens minf1g1
gens ring minf1g1
use ring minf1g1
sub(minf1g1, t_1 => 1)
gens ring oo

U = kk[x_1,x_3,z_1..z_6,y_0]
U = kk[x_2,x_4,z_1..z_6,y_0]
FU = sub(sub(minf1g1, t_1 => 1),U)
transpose gens FU
FU
radical ideal mingens ideal singularLocus FU
dim oo
I0 = ideal mingens sub(ISS, t_1 => 0);
isslist = decompose I0
apply(isslist, degree)
sum apply(isslist, degree)
apply(isslist, I -> minimalPresentation ideal mingens(I + ideal(sub(random(1,R),T),sub(random(1,R),T))))
isslist#0
transpose gens minimalPresentation  (isslist#0 + ideal(sub(random(1,R),T),sub(random(1,R),T)))
decompose ideal irk
irk = o103

use ring isslist#0
sub(isslist#0, {z_2 => z_2+2*y_0, z_1 => z_1+42*y_0, z_6 => z_6-25*y_0,z_5 =>  z_5 + 12*y_0})

------------------------


W = kk[x_1..x_6,z_1..z_6,y_0,y_1,t_1,  Weights => {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1,1}]
W = kk[x_1..x_6,z_1..z_6,y_0,y_1,  Weights => {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1}]
IW = sub(f1g, W)

r1 = random(1,W)
r2 = random(1,W)

sw = sub(IW, x_2 => -1) + ideal (x_2)
IWM = minimalPresentation (sw) ---, Exclude => {14})
f = sw.cache.minimalPresentationMap
CY = IWM + (ideal f (r1+1)) + (ideal f (r2+1))
CYM = minimalPresentation(CY)---, Exclude => {10})

dim CYM --- dim 3 :p

source gens gb CYM


A = ring CYM

(gens CYM)_6
gens A

AA = kk[x_1,z_6,y_1,z_2,z_3,x_3,y_0, Weights => {1,8,3,8,7,3,7}]
CYMAA = sub(CYM,AA)
source gens gb CYMAA
------ klart å finne passe liten grøbner-base. men ligner ikke på noen def av C(dP_6)!!
loadPackage "VersalDeformations"
CYMAA = ideal mingens CYMAA
res CYMAA
initialIdeal CYMAA

ideal mingens ideal singularLocus CYMAA
"cyikart1" << toString CYMAA << endl << close;

---
CYF =  f1g + random(1,R)  + random(1,R)

torus = ideal apply(flatten apply(apply(apply(flatten entries gens CYF, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)

ideal mingens oo
dim oo
-----------
transpose gens Y
G = sub(transpose gens Y, z_5 => z_5/2)
G = sub(G, z_1 => z_1/6)
G = sub(G, z_3 => z_3/3)
G = sub(G, y_0 => y_0/35)
G = transpose mingens ideal G
---
sub(Y,R)
torus = ideal apply(flatten apply(apply(apply(flatten entries gens sub(Y,R), monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)

--------
loadPackage "defMethods"
loadPackage "NormalToricVarieties"
ttt = transpose ((transpose makeCone torus)_{0..4})
P = convexHull ttt
Q = 2*P
vertices Q
v = (interiorLatticePoints (2*P))_0
Z = affineImage(Q,-v)
isReflexive Z
V = normalToricVariety (normalFan Z)
apply(apply(cones(5,normalFan Z), rays), m -> minors(5,m))

---

I1 = minimalPresentation (sub(f1g, x_1 => 1) + ideal x_1 + ideal(t_2..t_24))

I1

transpose gens oo
(numgens ring oo) - 24

