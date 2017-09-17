restart
--- mål: study the Z/3-action on X_t
A = QQ[w]/ideal(w^2+w+1)  --- the action is given by multiplication by a third root of unity
R = A[x_1..x_18]          --- in PP^17.
S = R[lambda]

--- these are the coordinats of omega * (...x_i...).
AB = transpose genericMatrix(R,3,3) | transpose genericMatrix(R,x_10,3,3)  -- a pair of 3x3 matrices
gAB = matrix{{x_1,   w*x_2,   w^2*x_3, x_10,  w*x_11,  w^2*x_12},   --- multiplied by g in Z/3.
             {w*x_4, w^2*x_5, x_6,    w*x_12, w^2*x_14, x_15},
             {w^2*x_7, x_8,    w*x_9, w^2*x_16, x_17, w*x_18}}
I = ideal (sub(AB,S) - lambda * sub(gAB,S))
fikspkts = saturate(saturate(I, lambda), sub(ideal(x_1..x_18),S))
LL = decompose fikspkts

IH =  ideal(x_9-2*x_11-2*x_13,2*x_6+2*x_8-x_10,x_5-2*x_12-2*x_16,2*x_3+2*x_7-x_14,2*x_2+2*x_4-x_18,x_1-2*x_15-2*x_17) -- idealet til span(fij)
IM = minors(2,AB_{0..2}) + minors(2,AB_{3..5})
L2 = decompose ideal mingens (fikspkts + IH + IM)
P = QQ[x_1..x_18]
apply(L2, I -> sub(ideal (I_*)_{0..16},P))




sub(Iinv,S)
apply(decompose fikspkts, I -> ideal mingens(sub(Iinv,S) + I))
apply(oo, dim)

-----
use P
M1 = matrix {{2*x_15+2*x_17,0,0},{0,0,x_6},{0,x_8,0}}
M2 = matrix{{2*x_6+2*x_8,0,0},{0,0,x_15},{0,x_17,0}}
decompose (minors(2,M1)+ minors(2,M2))

--- finn idealet til invariant plan
W = A[t_1..t_12]
fija = (i,j,a) -> (
    Eij  :=  (id_(W^3))_{i} * transpose (id_(W^3))_{j};
    Eij' :=  (id_(W^3))_{(-i-j) % 3} * transpose (id_(W^3))_{(-i-j) % 3};
    if (a == 0) then (
      Eij | 2 * Eij'
    )
    else (
      2 * Eij' | Eij
    )
  )
WW = R ** W
use W
spennplan = sub(t_1*fija(0,1,0) + t_2*fija(0,2,0) + t_3*fija(1,0,0) + t_4*fija(1,2,0) + t_5*fija(2,0,0) + 
    t_6*fija(2,1,0) + t_7*fija(0,1,1) + t_8*fija(0,2,1) + t_9*fija(1,0,1) + t_10*fija(1,2,1) + t_11*fija(2,0,1) + t_12*fija(2,1,1),WW)

Q = QQ[t_1..t_12,x_1..x_18]
--- idealet til invariante planet:
Iinv = eliminate(toList(t_1..t_12),sub(ideal (spennplan - sub(AB,WW)),Q))
dim Iinv
dim ring Iinv
gens ring Iinv
sub(Iinv, R)
use R

P = QQ[x_1..x_18]
sub(Iinv,P)
toString Iinv

---: regne på X_G
restart
kk = ZZ/3001
Z = QQ[x_1..x_12]

pars = {2,3,5}
pars = {1,1,1}
fija = (i,j,a) -> (
    Eij := (id_(Z^3))_{i} * transpose (id_(Z^3))_{j};
    Eij' :=  (id_(Z^3))_{(-i-j) % 3} * transpose (id_(Z^3))_{(-i-j) % 3};
    if (a == 0) then (
    Eij | pars#((i+j)%3) * Eij'
    )
    else (
    pars#((i+j)%3) * Eij' | Eij
    )
    )

MG = x_1*fija(0,1,0) + x_2*fija(0,2,0) + x_3*fija(1,0,0) + x_4*fija(1,2,0) + x_5*fija(2,0,0) + 
    x_6*fija(2,1,0) + x_7*fija(0,1,1) + x_8*fija(0,2,1) + x_9*fija(1,0,1) + x_10*fija(1,2,1) + x_11*fija(2,0,1) + x_12*fija(2,1,1)

--- virkning på P^11
Zw = Z[w,lambda]/ideal(w^2+w+1)

Zv = matrix{{x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,                    x_10,   x_11,x_12}}
gZv = matrix{{w*x_1,w^2*x_2,w*x_3,x_4,w^2*x_5,x_6,w*x_7,w^2*x_8,w*x_9,x_10,w^2*x_11,x_12}}

INVLIST = decompose saturate(ideal(Zv - lambda*gZv), sub(ideal(x_1..x_12),Zw))
G1 = x_4*fija(1,2,0)+x_6*fija(2,1,0)+x_10*fija(1,2,1) + x_12*fija(2,1,1)
G2 = x_1*fija(0,1,0)+x_3*fija(1,0,0)+x_7*fija(0,1,1) + x_9*fija(1,0,1)
G3 = x_2*fija(0,2,0)+x_5*fija(2,0,0)+x_8*fija(0,2,1) + x_11*fija(2,0,1)

decompose (minors(2,G1_{0..2}) + minors(2,G1_{3..5}) + ideal ((INVLIST#0)_*)_{0..7})
decompose (minors(2,G2_{0..2}) + minors(2,G2_{3..5}) + ideal ((INVLIST#1)_*)_{0..7})
decompose (minors(2,G3_{0..2}) + minors(2,G3_{3..5}) + ideal ((INVLIST#2)_*)_{0..7})
--- fra over: konkluderer med at fikspunktene til gruppa er kun de tolv punktene (0....:1:...).

---
--MG =  matrix {{2*x_10+2*x_12, x_1, x_2, 2*x_4+2*x_6, x_7, x_8}, {x_3,
--        2*x_8+2*x_11, x_4, x_9, 2*x_2+2*x_5, x_10}, {x_5, x_6, 2*x_7+2*x_9,
--        x_11, x_12, 2*x_1+2*x_3}}
IX = minors(2,MG_{0,1,2}) + minors(2,MG_{3,4,5})
IX = ideal mingens IX
time gens gb IX;

loadPackage("DefMethods", Reload => true)
fastSingularities IX

minimalPresentation(IX + ideal(x_2-1))
time singularLocus IX;

ISING = time mingens ideal singularLocus ideal mingens minimalPresentation(IX + ideal(x_2-1)); -- KAN REGNE MED EN TIME MINST (40 min for t=1)
ISING = radical ideal ISING



--"singsIX_t2" << toString ISING << endl << close -- SKRIV TIL FIL

loadPackage "MinimalPrimes"
installMinprimes()

decompose ISING

fastSingularities = method()
fastSingularities(Ideal) := I -> (
    R := ring I;
    n := numgens R;
    gensR := gens R;
    singlist := {};
    for i from 0 to 2 do { ---(n-1) do {
	affineChart := I + ideal(gensR_i - 1);
	print(i);
	print("Computing singular locus...");
	sing := ideal singularLocus ideal mingens minimalPresentation affineChart;
	print("computing radical of singular locus...");
	sing        = radical ideal mingens sing;
	inv         := affineChart.cache.minimalPresentationMap;
	singlist = singlist | {(homogenize(preimage(inv,sing),gensR_i))};
	};
    saturate intersect(singlist)
    )

fastSingularities(IX)
decompose oo
apply(oo, degree)
)
--ISING = value get "singsIX_t2";
---
loadPackage "LocalRings"
I1 =  minimalPresentation(IX + ideal(x_1-1))
setMaxIdeal ideal gens ring I1
LL = flatten entries localMingens gens I1
f = (LL)#2
use ring f
(f % x_9)
numerator((f - (f % x_9))/x_9)
LL' = apply(LL, g -> numerator sub(g, x_9 => x_10*x_11/numerator((f - (f % x_9))/x_9)))
f' = LL'#3
f' % x_7
numerator((f' - (f' % x_7))/x_7)
LL'' = apply(LL', g -> numerator sub(g, x_7 => x_8*x_12/numerator((f' - (f' % x_7))/x_7)))
f'' = LL''#0
f'' % x_6
numerator((f'' - (f'' % x_6))/x_6)
LL''' = apply(LL'', g -> numerator sub(g, x_6 => -(2*x_8*x_10*x_12+2*x_10*x_11*x_12-x_8*x_11)/numerator((f'' - (f'' % x_6))/x_6)))
LL'''#1

---- fire stk! så lurer på mo det ikke er C(dP6)???

---

--toString ideal  mingens( IX + ideal(x_12,x_11+2,x_10,x_9,x_8+2,x_6-1,x_2))
ideal(x_12,x_11+2,x_10,x_9,x_8+2,x_7,x_6-1,x_5,x_4,x_3,x_2,x_1-1)

PT = sub(MG,{x_12 => 0, x_11 => -2, x_10 => 0, x_9 => 0, x_8 => -2, x_7 => 0, x_6 => 1, x_5 => 0, x_4 => 0, x_3 => 0, x_2 => 0, x_1 => 1})
(decompose ISING)#1
--
(decompose ISING)#1
minIX = minimalPresentation(IX + ideal(x_1-1))
use ring minIX
IX0 = sub(minIX,  {x_10  => 1/2+x_10, x_8 => x_8 + 1/2})

loadPackage "LocalRings"
setMaxIdeal ideal gens ring minIX
fff = flatten entries localMingens gens IX0
(fff#1 - (fff#1 % x_6))/x_6
f2 = apply(fff, f -> numerator sub(f, x_6 => (-((fff#1 % x_6))/((fff#1 - (fff#1 % x_6))/x_6))))

(f2#0 - (f2#0 % x_9))/ x_9
f3 = apply(f2, f -> numerator sub(f, x_9 => -(f2#0 % x_9)/((f2#0 - (f2#0 % x_9))/x_9)))

Iloc = ideal f3
last decompose Iloc
transpose localMingens gens oo
-------
IX0' = ideal mingens ideal apply(0..8, i-> numerator sub((IX0_*)#i, x_2 => -(x_6*x_10-x_6+(1/2)*x_9-(1/2)*x_8*x_9)/(x_8*x_10+2*x_10*x_11-x_8-x_10-2*x_11+1)))
f = (IX0')_*#0
IX0'' = ideal mingens ideal apply(0..5, i-> numerator sub((IX0_*)#i, x_9 => (x_10*x_11-x_11)/numerator((f - (f % x_9))/x_9)))
f = (IX0'')_*#1
ideal mingens (IX0'' + ideal(x_11))
((ideal mingens (IX0'' + ideal(x_11)))_*)#1
decompose ideal oo
last oo

--
EPS = QQ[x_2,x_6,x_8,x_9,x_10,x_11,x_12][eps]/ideal(eps^2)
sub(sub(minimalPresentation(IX + ideal(x_1-1)),EPS), {x_10 => (1/2)*x_10-(1/2), x_8 => (1/2)*x_8-(1/2)})
map(EPS,EPS, {eps*x_2, eps*x_6, eps*x_8, eps*x_9, eps*x_10, eps*x_11,eps*x_12,eps})
oo sub(minimalPresentation(IX + ideal(x_1-1)), EPS)
oo sub(sub(minimalPresentation(IX + ideal(x_1-1)),EPS), {x_10 => (1/2)*x_10-(1/2), x_8 => (1/2)*x_8-(1/2)})
mingens oo
ideal oo

---
torus = ideal apply(flatten apply(apply(apply(flatten entries gens IX, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)
radical ideal mingens torus
loadPackage "Binomials"
toruskomps = BPD torus
toruskomps = select(toruskomps, I -> dim I == 1) --- har en Z_3-virkning!!! (altså den med f_ij)


