restart
R = QQ[a_1..a_5]
S = QQ[x,y,z]


f = map(S,R,{x*z^2,y*z^2,x*y*z,y*x^2,x*y^2})
I = ker f

loadPackage "VersalDeformations"

sub1 = sub(I, a_1 => 1)
J1 = minimalPresentation(sub1)
g = sub1.cache.minimalPresentationMap

T1 = CT^1(-2, gens J1)
T2 = CT^2(-2, gens J1);
(f1,r1,g1,c1) = versalDeformation(gens J1, T1, T2);
sum f1
I1 = ideal(sub(sum f1, t_1 => 1))
fkand = transpose mingens homogenize(preimage(g, I1),a_1)

fk = ideal fkand
hilbertPolynomial fk ==hilbertPolynomial I

X = Proj(R/I)
X1 = Proj(R/fk)

apply(0..3, i-> HH^i(tangentSheaf X1))
apply(0..3, i-> HH^i(tangentSheaf X))
apply(0..3, i-> HH^i(tangentSheaf Y))
apply(0..3, i-> HH^i(cotangentSheaf Y))

loadPackage "NormalToricVarieties"

V = normalToricVariety matrix {{1,0,-1,0},{0,1,0,-1}}
W = makeSmooth V

OO (-sum toList apply(0..7, i-> W_i))
OO ( (-6*W_1+14*W_2+4*W_3+5*W_5+10*W_7))

(inverse fromCDivToPic W)*v
C =  posHull fromCDivToPic W
rays dualCone C
v = 12*(interiorPoint convexHull rays dualCone C)

projEmb = (D) -> (X = variety D;
    L = OO D;
    m = rank HH^0(X,L);
    S = ring X;
    R = QQ[y_0..y_(m-1)];
    phi = map(S,R,basis(-first degrees L, S));
    kernel phi)
projEmb     


decompose ideal singularLocus X1

apply(1..5, i -> decompose ideal singularLocus minimalPresentation (fk + ideal(a_i - 1)))

fglatt = ideal sub((sum first versalDeformation(gens I, CT^1(0, gens I), CT^2(0, gens I))), {t_1 => 2, t_2 => 2, t_3 => 2, t_4 => 2})

Y = Proj(R/fglatt)


radical ideal  singularLocus fglatt

sub2 = sub(fk, a_2 => 1)
J1 = minimalPresentation(sub1)
g = sub1.cache.minimalPresentationMap

T1 = CT^1(0, gens J1)
T2 = CT^2(0, gens J1);
(f1,r1,g1,c1) = versalDeformation(gens J1, T1, T2);

I1 = ideal(sub(sum f1, t_1 => 1))
fkand = transpose mingens homogenize(preimage(g, I1),a_1)

fk = ideal fkand
hilbertPolynomial fk ==hilbertPolynomial I

X = Proj(R/I)
X1 = Proj(R/fk)

apply(0..3, i-> HH^i(tangentSheaf X1))
apply(0..3, i-> HH^i(tangentSheaf X))

decompose ideal singularLocus X1

apply(1..5, i -> decompose ideal singularLocus minimalPresentation (fk + ideal(a_i - 1)))

CT^1(gens fk)
