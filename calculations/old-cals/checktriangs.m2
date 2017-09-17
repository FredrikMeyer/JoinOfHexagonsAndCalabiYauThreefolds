restart
loadPackage("VersalDeformations")
loadPackage("SimplicialComplexes")

{*
The file listoftriangs33.txt contains all regular, fine, connected triangulations of
the polytope (PxP)^dual. There are 4441 of them.

In this file I compute all their T^1 and T^2 modules.
*}

listoftriangs = value get "listoftriangs33.txt"; --- 4441 triangulations

R = QQ[x_0..x_12]

-- f takes a triangulation and produces the corresponding
-- Stanley-Reisner ideal.
f = method();
f List := t -> ideal simplicialComplex apply(t, s -> product apply(s, v -> x_v));

-- It has the correct/expected Hilbert polynomial.
hilbertPolynomial(f listoftriangs#0)

-- I compute all the T1's and T2's.
T1s = apply(listoftriangs, t -> rank source CT^1(0, gens f t))
T2s = apply(listoftriangs, t -> rank source CT^2(0, gens f t))
min T1s -- 60
min T2s -- 53

core = method();
apply(listoftriangs,k -> #fold(s -> (set s) * (set k#0), k))
lll = (listoftriangs#0)
I = f lll
A = R/I

select({0,1,2,3,4,5,6,7,8,9,10,11,12}, i ->  all(lll, F -> member(i,F)))

apply(listoftriangs, lll ->  #select({0,1,2,3,4,5,6,7,8,9,10,11,12}, i ->  all(lll, F -> member(i,F))))

B = select(listoftriangs, lll ->  #select({0,1,2,3,4,5,6,7,8,9,10,11,12}, i ->  all(lll, F -> member(i,F))) > 0)
lll = oo
select(listoftriangs, t -> betti res f t == B)
--- bare én har samme res som den med kjegle
isGorenstein = method(); -- method to check if a triangulation is Gorenstein (assumes S=Core(S), i.e S is not a cone)
K =  simplicialComplex monomialIdeal f lll#0
apply(listoftriangs, t -> { simplicialComplex monomialIdeal f t;
	all(0..4, i-> all(flatten entries faces(3, K), F -> 
		rank HH_i(link(K,F)) == i));
	})

A = R/ideal(x_12)
I = f lll#0
J = sub(I,A)
CT^1(0, gens J)
CT^2(0, gens J)

use R
betti res f listoftriangs#345

flatten entries faces(3, simplicialComplex monomialIdeal f lll)
rank HH_4(simplicialComplex monomialIdeal f lll)

betti res f listoftriangs#0
select(listoftriangs, t -> CT^1(-2,f t) != 0)
-- Note: See relationt1t2.py and t1vst2.png for a plot of these values.
-- There seem to be a linear relation.

-- The minimal T^2-dimension is 53. I select every triangulation with
-- this obstruction space.
minIndices = select(0..(#listoftriangs-1), i -> T2s#i == 53); -- 144 stk
minT2triangs = apply(minIndices, i -> f listoftriangs#i);

-- I compute the T1's of these.
minT1triangs = apply(minT2triangs,  t -> rank source CT^1(0, gens t))
-- All of these have T^1=60, which is the minimum dimension.

--- TODO: find Gorenstein trians



-- I choose one of these.
tr = listoftriangs#5  -- number 5, 7
J = f trp
T1 = CT^1(0, gens J);
T2 = CT^2(0, gens J);

-- I create a degegenerated Calabi-Yau.
IX0 = ideal gens(J + ideal sum toList apply(0..12, i -> x_i))
IX0 = minimalPresentation IX0
T1X0 = CT^1(0, gens IX0); -- 53 (8 in degree -1)
T2X0 = CT^2(0, gens IX0); -- 28 (21 in degree -1)

-- I deform it to second order
(f1,r1,g1,c1) = versalDeformation(gens IX0, T1X0, T2X0, HighestOrder => 2, Verbose => 5, SanityCheck => false,
    SmartLift => false);
--- g1 is veeeery complicated (less so without usings mingens in the def of IX0)
(transpose sum f1)_0_2

---
ideal sum g1
loadPackage "minimalPrimes"
installMinprimes()
decompose ideal sum g1
LL = oo
--

-- Instead deform the original triang
(f1,r1,g1,c1) = versalDeformation(gens J, T1, T2, HighestOrder => 5, Verbose => 3, SmartLift => false, SanityCheck => false);
transpose gens ideal sum g1 -- Not too complicated

loadPackage "MinimalPrimes"
installMinprimes()

transpose gens ideal sum g1 --- Too complited to do something with
