restart
Z = QQ[x_1..x_12]

pars = {2,3,5}
fija = (i,j,a) -> (
    Eij := (id_(Z^3))_{i} * transpose (id_(Z^3))_{j};
    Eij' :=  (id_(Z^3))_{(-i-j) % 3} * transpose (id_(Z^3))_{(-i-j) % 3};
    if (a == 0) then (
    Eij | pars#((i-j)%3) * Eij'
    )
    else (
    pars#((i-j)%3) * Eij' | Eij
    )
    )

MG = x_1*fija(0,1,0) + x_2*fija(0,2,0) + x_3*fija(1,0,0) + x_4*fija(1,2,0) + x_5*fija(2,0,0) + x_6*fija(2,1,0) +
     x_7*fija(0,1,1) + x_8*fija(0,2,1) + x_9*fija(1,0,1) + x_10*fija(1,2,1) + x_11*fija(2,0,1) + x_12*fija(2,1,1)

IX = minors(2,MG_{0..2}) + minors(2,MG_{3..5})
time gens gb IX;

-- automorphisms of X
sigma = map(Z,Z,{x_4,x_3,x_6,x_5,x_1,x_2,x_10,x_9,x_12,x_11,x_7,x_8}) -- Z/3-virkning på E, e_i -> e_(i+1)
tau   = map(Z,Z,{x_7,x_8,x_9,x_10,x_11,x_12,x_1,x_2,x_3,x_4,x_5,x_6}) -- bytter om E*E-faktoreneb
mu    = map(Z,Z,{x_3,x_5,x_1,x_6,x_2,x_4,x_9,x_11,x_7,x_12,x_8,x_10}) -- bytter om på e_ij -> e_ji
sigma IX == IX -- true
tau IX == IX -- true

loadPackage "MinimalPrimes"
installMinprimes()

IX1 = IX + (x_1-1)
sing1 = time ideal singularLocus ideal mingens minimalPresentation IX1;
sing1       = time radical ideal mingens sing1;
inv1         = time IX1.cache.minimalPresentationMap;
preim1 = time saturate homogenize(preimage(inv,sing),x_1)
IX2 = IX + (x_2-1);
sing2 = time ideal singularLocus ideal mingens minimalPresentation IX2;
sing2       = time radical ideal mingens sing2;
inv2         = time IX2.cache.minimalPresentationMap;
preim2 = time saturate homogenize(preimage(inv2,sing2),x_2);
IX3 = IX + (x_3-1);
sing3 = time ideal singularLocus ideal mingens minimalPresentation IX3;
sing3       = time radical ideal mingens sing3;
inv3         = time IX3.cache.minimalPresentationMap;
preim3 = time saturate homogenize(preimage(inv3,sing3),x_3);
IX4 = IX + (x_4-1);
sing4 = time ideal singularLocus ideal mingens minimalPresentation IX4;
sing4       = time radical ideal mingens sing4;
inv4         = time IX4.cache.minimalPresentationMap;
spreim4 = time saturate homogenize(preimage(inv4,sing4),x_4);
IX5 = IX + (x_5-1);
sing5 = time ideal singularLocus ideal mingens minimalPresentation IX5;
sing5       = time radical ideal mingens sing5;
inv5         = time IX5.cache.minimalPresentationMap;
preim5 = time saturate homogenize(preimage(inv5,sing5),x_5);

decompose intersect {preim,tau preim, mu preim, tau mu preim, preim2,tau preim2, mu preim2, tau mu preim2,
    preim3, tau preim3, mu preim3, tau mu preim3, preim4, tau preim4, mu preim4, tau mu preim4,
    preim5, tau preim5, mu preim5, tau mu preim5}
apply(oo, degree)
sum oo
--- supert: betyr at vi har 48 sings!!-

singsSS = intersect decompose intersect {preim,tau preim, mu preim, tau mu preim, preim2,tau preim2, mu preim2, tau mu preim2,
    preim3, tau preim3, mu preim3, tau mu preim3, preim4, tau preim4, mu preim4, tau mu preim4}


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

f = (I1_*)#2
use ring f
(f % x_9)
numerator((f - (f % x_9))/x_9)
LL = apply(I1_*, g -> numerator sub(g, x_9 => x_10*x_11/numerator((f - (f % x_9))/x_9)))

localMingens matrix{LL}

-- automorphisms of X
--
loadPackage "MinimalPrimes"
installMinprimes()


divlist = decompose ideal mingens (IX  + ideal(x_1))
D = last oo

degree D
intersect divlist == IX + ideal(x_1)

hilbertPolynomial divlist#0
hilbertPolynomial D

D
minimalPresentation D
radical oo
betti res oo
betti res first divlist

mingens ideal singularLocus D

apply(divlist, degree)

decompose(IX + random(1,Z))


complist = {}
for i from 1 to 12 do (
    complist = complist | decompose (IX + ideal(x_i));
    print(i);
    )
#complist
#unique complist
complist

complist = unique complist
apply(complist, degree)

#complist

apply(complist, I -> dim singularLocus I)

last complist
betti res oo
betti res minimalPresentation last complist


-- finne fikspunkter til sigma
sigma

W = Z[lambda]


sub(IX,W) + ideal (matrix{ gens Z} - lambda * sub(matrix sigma,W))
decompose oo
decompose ideal (matrix{ gens Z} - lambda * sub(matrix sigma,W))
