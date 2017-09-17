restart
loadPackage "SimplicialComplexes"
loadPackage "VersalDeformations"
kk = ZZ/109
R = kk[x_1..x_8]

stanleyReisner = L -> (
    facetList := apply(L, l -> product apply(l, i -> x_i));
    S := simplicialComplex facetList;
    ideal S
    )

T1T2dims = I -> (
    T1 := CT^1(0, gens I);
    T2 := CT^2(0, gens I);
    (rank source T1, rank source T2)
    )

triangs  = value get "grunbaum_list.m2";
noObstructions = select(triangs, T -> (T1T2dims stanleyReisner T)_1 == 0)

apply(noObstructions, T -> T1T2dims stanleyReisner T)
apply(noObstructions, T -> betti res  stanleyReisner T)
apply(noObstructions, T -> prune HH simplicialComplex monomialIdeal stanleyReisner T)

T1 = noObstructions#0
T2 = noObstructions#1 -- nr 16
T3 = noObstructions#2
T4 = noObstructions#3
T5 = noObstructions#4
T6 = noObstructions#5


I1 = stanleyReisner T1
I2 = stanleyReisner T2
I3 = stanleyReisner T3
I4 = stanleyReisner T4
I5 = stanleyReisner T5
I6 = stanleyReisner T6--





hilbertPolynomial(I6, Projective => false)
16*12/3
T1T2dims stanleyReisner T5

T1 = CT^1(0, gens I5)
T1 = submatrix'(T1, {20,21,26,27,28,38,39,40})
T1 = T1 * transpose random(R^1,R^62, Density => 0.3)
(f1,r1,g1,c1) = versalDeformation(gens I5, T1, CT^2(0, gens I5), HighestOrder => 20, SmartLift => false, PolynomialCheck => false);
I2t = ideal sub(transpose sum f1, toList apply(1..1, i-> (t_i => random(0,R))))
--I2t = ideal sub(transpose sum f1, toList apply(1..72, i-> (t_i => random(0,R))))

decompose I2t

M = mingens ideal singularLocus I2t;
dim ideal M
radical ideal M
T6

I1t = sub(ideal sum f1, toList apply(1..83, i-> (t_i => random(0,R))))

X = Proj(R/I1t)
CT^1(0, gens I1t)

CT = cotangentSheaf X;


betti res I1t
betti res prune (I1t/I1t^2)

NX = prune sheaf (I1t/I1t^2);
HH^2(NX)

prune sheaf Hom(I1t/I1t^2,R/I1t)
HH^0(oo)
HH^1(o68) --- Så Hodge tall til nr 1 er (1,76)...

dim ideal mingens ideal singularLocus I1
---
t2 = triangs#1
I2 = stanleyReisner t2

T1 = CT^1(0, gens I2)
T1 = submatrix'(T1,{94,93,90,89})
T1 = submatrix(T1, {0,3,4,5,10,15,20,25,30,31,32,33,34,35,40,45,50,51,52,53,54,55,60,65,70,75,80,85,90,93,92,91})
T2 = CT^2(0, gens I2)

(f1,r1,g1,c1) = versalDeformation(gens I2, T1, T2, HighestOrder => 200);

It = ideal sub((sum f1, toList apply(1..32, i-> t_i => random(0,R))))

ideal singularLocus It;




loadPackage "MinimalPrimes"
installMinprimes()

decompose oo

