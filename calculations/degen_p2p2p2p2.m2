restart
{*
Idea: try to deform the ideal of X_00 inside the ideal of
P^2 * P^2 instead of directly. Maybe the equations turn out easier.

Answer is "yes", but not easy enough.
   *}

R = QQ[x_1..x_18]

M1 = genericMatrix(R,3,3)
M2 = genericMatrix(R,x_10,3,3)

I = minors(2,M1) + minors(2,M2)
IX = ideal(x_1,x_5,x_9,x_10,x_14,x_18) + I 

loadPackage "VersalDeformations"

T1 = CT^1(0, gens IX);
T2 = CT^2(0, gens IX);


(f1,r1,g1,c1) = versalDeformation(gens IX, T1,T2, HighestOrder => 2, SmartLift => false);
transpose sum f1
sum g1
loadPackage "MinimalPrimes"
installMinprimes()

I1 = ideal ((transpose sum g1)_{0..35})
ideal mingens I1
decompose oo
decompose ideal sum g1
first I1_*
varsused = unique flatten apply(I1_*, s-> flatten apply(decompose ideal (monomials s), I -> I_*))

TT1 = QQ[varsused]
I1TT1 = sub(I1, TT1)
L1 = decompose oo

