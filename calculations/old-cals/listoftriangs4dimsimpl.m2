restart

loadPackage "SimplicialComplexes" 

R = QQ[a_0..a_12]
L = value get "alltriangsm2.txt";


listOfIdeals = apply(L, T->  ideal simplicialComplex apply(T, f -> product apply(f, i -> a_i)));

h1 = sum sum entries vars R

loadPackage "VersalDeformations"

T1s = apply(listOfIdeals, I-> rank source CT^1(0, gens I));
#listOfIdeals
T1s = {}
for i from 0 to 4440 do (
    T1s = T1s | {rank source CT^1(0, gens listOfIdeals#i)};
    if (i % 100 == 0) then (
	print(i);
	);
    );

min T1s -- 60

T2s = {}
for i from 0 to 4440 do (
    T2s = T2s | {rank source CT^2(0, gens listOfIdeals#i)};
    if (i % 100 == 0) then (
	print(i);
	print min T2s;
	);
    );
min T2s -- 

apply(listOfIdeals, I-> rank source CT^2(0, gens I))
positions(T2s, i -> i == 53) == positions(T1s, i -> i == 60)
min oo

-------------
Ik = listOfIdeals#100
(f1, r1, g1, c1) = versalDeformation(gens Ik, CT^1(0, gens Ik), CT^2(0, gens Ik), HighestOrder => 2, SmartLift => false);

loadPackage "MinimalPrimes"
installMinprimes()

decompose  ideal mingens ideal g1
--fVector simplicialComplex monomialIdeal listOfIdeals#100    

loadPackage "Binomials"
decompose ideal select((ideal g1)_*, f -> #terms f <= 2)

------------

Ik = listOfIdeals#100
T1 = CT^1(0, gens Ik);
T2 = CT^2(0, gens Ik);

(f1,r1,g1,c1) = versalDeformation(gens Ik, T1, T2, HighestOrder => 2, SmartLift => false);

IG = (ideal sum g1)
transpose mingens IG
loadPackage "Binomials"

lms =  BPD ideal select(IG_*, f -> #terms f <= 2)

toString IG



