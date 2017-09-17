restart
R1 = QQ[x_1..x_6,y_0]
R2 = QQ[z_1..z_6,y_1]

loadPackage "SimplicialComplexes"
loadPackage "VersalDeformations"
loadPackage "Depth"

S1 = simplicialComplex {x_1*x_2*x_3, x_1*y_0*x_3, x_3*x_4*y_0, x_6*x_1*y_0, x_6*x_4*y_0, x_6*x_5*x_4}
S2 = simplicialComplex {z_1*z_2*z_3, z_1*y_1*z_3, z_3*z_4*y_1, z_6*z_1*y_1, z_6*z_4*y_1, z_6*z_5*z_4}
I1 = ideal S1
I2 = ideal S2

{*
Experimenting with another triangulation. This is a join of two non-hexagons.
*}

R = QQ[first entries vars R1 | first entries vars R2]
I = sub(I1,R) + sub(I2,R)
depth (R/I) -- 6, so depth=dim, so C-M
{*
Below: the second T11 is obtained by decomposing the second order terms of the relation
ideal of the versal deformation of I. Code for this is below.

IMPORTANT: do not change term order for this to work. The last T11 is a submatrix of the
second T11. So to get the last, compute nr 1 and 2 first.

The second T11 gives a 44-dim family over the dual numbers, which can be lifted over higher k[t]/t^n.
I havent been able to lift it to k[t]. 

The third T11 contains paramters lifting to a 8-dim family of toric varieties. This makes sense, because
the Batyrev-Borisov construction gives a toric with CY with Hodge numbers (44,8). 
*}
T11 = CT^1(0, gens I); -- 58-dim
T11 = submatrix'(CT^1(0, gens I),{26,25,24,23,22,21,20,53,52,51,50,49,48,47}) --- svarer til den største komponenten til defrommet
--T111 = submatrix'(T11,{10,11,13,14,15,16,17,18,19,30,31,32,33,34,35,36,37,38})
T11 = submatrix(T11, {2,4,6,8,12,14,18,21,25,29,31,25,39,40,41,42,43}) --- try to get Y_0
T11 = submatrix(T11, {43,41,2,29,42,40,12,39}) --- gir del Pezzo**del Pezzo
T11 = submatrix'(T11,{2,12,29,39}) -- fjerne kvadrater

(f1,r1,g1,c1) = versalDeformation(gens I, T11, CT^2(0, gens I), HighestOrder => 2, SmartLift => false);
(f1,r1,g1,c1) = versalDeformation(gens I, T11, CT^2(0, gens I),HighestOrder =>5, Verbose => 10);--- , SmartLift => false);


ideal sum g1
loadPackage "MinimalPrimes"
loadPackage "Binomials"
installMinprimes()
decompose ideal sum g1
BPD ideal sum g1
--(transpose sum f1)_0_18
transpose ( (sum f1)*(sum r1))+(sum c1)*sum(g1)==0  --:D


defI = ideal sub(sum f1, toList apply(1..17, i-> t_i => -1))
minimalPresentation(defI + ideal(x_6-1))
ideal mingens ideal singularLocus oo
dim oo
-- not smooth. But possibly the same deformation as the one found in
-- another file (which one?)



f1' = apply(f1, f -> sub(f, toList apply(1..44, i-> t_i => t_1)))
r1'= apply(r1, f -> sub(f, toList apply(1..44, i-> t_i => t_1)));
g1' = apply(g1, f -> sub(f, toList apply(1..44, i-> t_i => t_1)));
c1' = apply(c1, f -> sub(f, toList apply(1..44, i-> t_i => t_1)));


(sub(transpose sum f1,{t_20 => 0, t_19 => 0, t_18 => 0, t_17 => 0, t_16 => 0, t_15 => 0, t_14 => 0, t_11 =>0,t_12 => 0,t_39 => 0, t_38 => 0,t_37 => 0, t_36 => 0, t_35 => 0, t_34 => 0,t_33 => 0, t_32 => 0, t_31 => 0}))_0_17

(transpose sum f1)
defI =  ideal transpose  mingens ideal sub((transpose sum f1), toList apply(1..8, i -> t_i => -i))


transpose gens defI

mingens ideal singularLocus minimalPresentation(defI + ideal(x_3 -1))
SingX3 = oo
dim ideal oo


transpose sum f1

loadPackage "MinimalPrimes"

sub(transpose sum f1, apply({1,2,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,31,32,33,34,35,36,37,38,39,
	    21,22,23,24,25,26,27,28,29}, i -> (t_i => 0)))

loadPackage "Binomials"
IG = ideal sum g1

IG1 = ideal ((IG_*)_{0..37})
IG2 = ideal ((IG_*)_{38..75})

T = QQ[t_1..t_27]
IGT1 = sub(IG1,T)
kList = BPD IGT1
apply(kList, dim)
select(kList, I -> I == radical I)
apply(oo, dim)

first kList
select(kList, I -> (I == radical I))

T2 = QQ[t_28..t_58]
IGT2 = ideal mingens sub(IG,T2)
BPD IGT2
kList2 = oo
select(kList2, I -> I == radical I)
apply(oo, dim)

TT = QQ[t_1..t_58]
kListTT = apply(kList, I -> sub(I,TT))
kList2TT = apply(kList2, I -> sub(I,TT))

IGTT  = ideal sub(sum g1,TT)
time intersect(intersect kListTT, intersect(kList2TT)) ==  IGTT --- long time;



comps = apply(toList ((set kListTT) ** (set kList2TT)), S -> ideal mingens (S#0 +S#1));
apply(comps, dim)
maxComp = select(comps, I -> dim I == 44)
toString oo
----{ideal(t_54,t_53,t_52,t_51,t_50,t_49,t_48,t_27,t_26,t_25,t_24,t_23,t_22,t_21)}

T11 = submatrix'(CT^1(0, gens I),{53,52,51,50,49,48,47,26,25,24,23,22,21,20,0,2,11,1,10,12,34,35,36,44,45,46,54,55,56,57}) -- fjerne frie variable
T11 = submatrix'(CT^1(0, gens I),{53,52,51,50,49,48,47,26,25,24,23,22,21,20,0,2,11,1,10,35,36,44,46,54,55,56}) -- fjerne frie variable
T11 = submatrix'(CT^1(0, gens I),{53,52,51,50,49,48,47,26,25,24,23,22,21,20})

T11 = submatrix(T11, {43,41,2,29,42,40,12,39}) --- gir del Pezzo**del Pezzo

(f1,r1,g1,c1) = versalDeformation(gens I, T11, CT^2(0, gens I), HighestOrder => 5, SmartLift => true);



CT^2(0,gens sub(ideal transpose sum f1, toList apply(1..20, i -> t_i => -1) | toList apply(21..32, i-> t_i => 0)))
CT^1(0,gens sub(ideal transpose sum f1, toList apply(1..8, i -> t_i => -1)))















S = simplicialComplex monomialIdeal I

flatten entries vars R
apply(flatten entries vars R, v -> #select(flatten entries faces(1,S), f -> (ideal f) : v != ideal f))
sum oo

X0 = Proj(R/((ideal S) + (y_0+y_1+x_1), (x_1+z_1+z_2)))

dim ideal mingens (ideal X0)


----
restart
R = QQ[x_0..x_13]

L = {(0, 1, 3, 5, 7, 9), (6, 7, 8, 9, 10, 12), (2, 3, 6, 7, 10, 12), (2, 3, 4, 6, 9, 10), (0, 2, 3, 6, 7, 12), (0, 2, 3, 6, 9, 12), (0, 1, 3, 7, 9, 12), (3, 4, 6, 9, 10, 11), (2, 7, 8, 9, 10, 12), (0, 2, 6, 7, 8, 12), (2, 3, 5, 7, 10, 12), (2, 3, 5, 9, 10, 12), (0, 1, 6, 7, 9, 12), (1, 3, 6, 7, 9, 12), (2, 5, 7, 9, 10, 12), (0, 1, 6, 7, 8, 9), (0, 6, 7, 8, 9, 12), (0, 2, 6, 8, 9, 12), (0, 2, 3, 4, 6, 9), (3, 6, 7, 10, 11, 12), (5, 7, 9, 10, 11, 12), (3, 5, 7, 9, 11, 12), (2, 3, 6, 9, 10, 12), (3, 5, 7, 10, 11, 12), (0, 2, 7, 8, 9, 12), (0, 2, 3, 5, 7, 12), (0, 2, 5, 7, 9, 12), (6, 7, 8, 9, 10, 11), (3, 6, 7, 9, 11, 12), (3, 6, 9, 10, 11, 12), (6, 7, 9, 10, 11, 12), (0, 1, 3, 6, 7, 12), (3, 5, 9, 10, 11, 12), (2, 6, 7, 8, 10, 12), (2, 6, 8, 9, 10, 12), (0, 3, 5, 7, 9, 12), (0, 2, 3, 5, 9, 12), (0, 1, 3, 6, 9, 12), (0, 1, 3, 4, 6, 9)}

L  = {(2, 3, 5), (1, 3, 4), (2, 3, 6), (1, 3, 6), (1, 2, 6), (0, 1, 2)}
R = QQ[x_0..x_6]



ls = {{(2, 3, 5), (1, 3, 4), (2, 3, 6), (1, 3, 6), (1, 2, 6), (0, 1, 2)}, {(1, 3, 4), (2, 5, 6), (1, 3, 6), (1, 2, 6), (0, 1, 2), (3, 5, 6)}, {(3, 4, 6), (2, 3, 6), (0, 1, 2), (2, 3, 5), (1, 2, 6), (1, 4, 6)}, {(2, 3, 5), (1, 3, 4), (2, 3, 6), (0, 2, 6), (1, 3, 6), (0, 1, 6)}, {(3, 4, 6), (2, 3, 6), (0, 2, 6), (2, 3, 5), (1, 4, 6), (0, 1, 6)}, {(1, 3, 4), (0, 2, 6), (2, 5, 6), (1, 3, 6), (3, 5, 6), (0, 1, 6)}, {(3, 4, 6), (0, 1, 2), (2, 5, 6), (1, 2, 6), (1, 4, 6), (3, 5, 6)}, {(3, 4, 6), (0, 2, 6), (2, 5, 6), (1, 4, 6), (3, 5, 6), (0, 1, 6)}, {(3, 4, 6), (2, 3, 6), (0, 2, 6), (2, 3, 5), (0, 1, 4), (0, 4, 6)}, {(0, 1, 2), (2, 5, 6), (3, 4, 5), (1, 2, 6), (1, 4, 6), (4, 5, 6)}, {(1, 3, 4), (1, 3, 6), (0, 5, 6), (3, 5, 6), (0, 1, 6), (0, 2, 5)}, {(1, 4, 6), (0, 2, 6), (2, 5, 6), (3, 4, 5), (0, 1, 6), (4, 5, 6)}, {(3, 4, 6), (0, 2, 6), (2, 5, 6), (0, 1, 4), (0, 4, 6), (3, 5, 6)}, {(3, 4, 6), (0, 1, 6), (0, 5, 6), (1, 4, 6), (3, 5, 6), (0, 2, 5)}, {(0, 2, 6), (2, 5, 6), (3, 4, 5), (0, 1, 4), (0, 4, 6), (4, 5, 6)}, {(3, 4, 6), (0, 1, 4), (0, 4, 6), (0, 5, 6), (3, 5, 6), (0, 2, 5)}, {(0, 1, 6), (0, 5, 6), (3, 4, 5), (1, 4, 6), (4, 5, 6), (0, 2, 5)}, {(0, 5, 6), (3, 4, 5), (0, 1, 4), (0, 4, 6), (0, 2, 5), (4, 5, 6)}} --- <- alle regulære trianguleringer av sekskanten 

I = ideal simplicialComplex apply(L, f -> product toList apply(f, v -> x_v))
Is = apply(ls, l -> ideal simplicialComplex apply(l, f -> product toList apply(f, v -> x_v)))

tr = {{0,1,2,3,4,5}, {0,1,2,4,5,6}, {0,1,2,5,6,7}, {0,1,2,6,7,8}, {1,2,3,4,5,12}, {1,2,3,4,9,12}, {1,2,3,5,9,12}, {1,2,4,5,6,12}, {1,2,4,6,9,12}, {1,2,5,6,7,12}, {1,2,5,7,9,12}, {1,2,6,7,8,12}, {1,2,6,8,9,12}, {1,2,7,8,9,12}, {1,3,4,5,9,12}, {1,4,5,6,9,12}, {1,5,6,7,9,12}, {1,6,7,8,9,12}, {2,3,4,5,10,12}, {2,3,4,9,10,12}, {2,3,5,9,10,12}, {2,4,5,6,10,12}, {2,4,6,9,10,12}, {2,5,6,7,10,12}, {2,5,7,9,10,12}, {2,6,7,8,10,12}, {2,6,8,9,10,12}, {2,7,8,9,10,12}, {3,4,5,9,10,11}, {3,4,5,9,10,12}, {4,5,6,9,11,12}, {4,5,6,10,11,12}, {4,5,9,10,11,12}, {4,6,9,10,11,12}, {5,6,7,9,11,12}, {5,6,7,10,11,12}, {5,7,9,10,11,12}, {6,7,8,9,11,12}, {6,7,8,10,11,12}, {6,8,9,10,11,12}, {7,8,9,10,11,12}}
I = ideal simplicialComplex apply(tr, f -> product toList apply(f, v -> x_v))
hilbertPolynomial I


unique apply(Is, I -> hilbertPolynomial I)

apply(Is, I -> source CT^1(0, gens I))
apply(Is, I -> source CT^2(0, gens I))

apply(select(Is, I -> source CT^2(0, gens I) == 0), I -> source CT^1(0, gens I))

----
loadPackage "defMethods"
toricIdeal matrix{{0, 0, 0, 0, 0, 0, -2, 2, 0, -2, 0, 2}, {0, 0, 0, 0, 0, 0, 0, 0, -2, -2, 2, 2}, {-2, 2, 0, -2, 0, 2, 0, 0, 0, 0, 0,0}, {0, 0, -2, -2, 2, 2, 0, 0, 0, 0, 0, 0}, {-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1},{1,1,1,1,1,1,1,1,1,1,1,1}}
transpose mingens oo
I


----
-- def av del pezzo ** 2

IdP = ideal sub(sum f1, toList apply(1..8, i -> t_i => -1))

T1dP6 = submatrix(CT^1(0, gens IdP),{1,5,8,12,22,26,29,33})

(f1,r1,g1,c1) = versalDeformation(gens IdP, T1dP6, CT^2(0, gens IdP));

decompose ideal mingens ideal sum g1
dP6def = ideal sub(transpose sum f1, {t_6 => 1, t_8 => 1, t_5 => 2, t_7 => 2, t_2 => 3, t_4 => 3, t_1 => 4, t_3 => 4})
dP6def = ideal sub(transpose sum f1, {t_6 => 1, t_8 => 1, t_5 => 1, t_7 => 1, t_2 => 1, t_4 => 1, t_1 => 1, t_3 => 1})

for i from 1 to 6 do (
    print dim ideal mingens ideal singularLocus minimalPresentation(dP6def + ideal(x_i -1)) --- osv 
)

transpose gens dP6def

--gens gb dP6def
loadPackage "CharacteristicClasses"
CSMClass dP6def
CSMClass IdP

mingens ideal singularLocus dP6def

singlist = {}
for i from 1 to 6 do {
    sz = sub(dP6def, z_i => 1);
    sx = sub(dP6def, x_i => 1);
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

sings = intersect singlist
sing1 = first decompose sings

singlist = decompose sings
dim(singlist#2 + singlist#5)

subsets({0,1,2,3,4,5},2)
apply(subsets({0,1,2,3,4,5},2), S -> dim(singlist#(S#0) + singlist#(S#1)))
--- snittgraf: 123, 054
mingens(singlist#1 + singlist#3) == 
mingens(singlist#1 + singlist#2) == mingens (singlist#2 + singlist#3)

h1 = x_1+x_2+x_3+x_4+x_5+x_6+z_1+z_2+z_3+z_4+z_5+z_6
h2 = x_1-x_2+x_3-x_4+x_5-x_6+z_1-z_2+z_3-z_4+z_5-z_6


sigma = map(R,R,{x_2,x_1,x_6,x_5,x_4,x_3,y_0,z_2,z_1,z_6,z_5,z_4,z_3,y_1})
tau   = map(R,R,{x_5,x_4,x_3,x_2,x_1,x_6,y_0,z_5,z_4,z_3,z_2,z_1,z_6,y_1})
alfa  = map(R,R,{z_1,z_2,z_3,z_4,z_5,z_6,y_1,x_1,x_2,x_3,x_4,x_5,x_6,y_0})

alfa dP6def == dP6def

CT^2(0, gens dP6def)
--------

transpose mingens ideal sum g1
transpose mingens IG1
