restart
loadPackage "Depth"
loadPackage "VersalDeformations"
loadPackage "SimplicialComplexes"
-- Kuhnels 9 vertex
T = {(2, 5, 6, 7, 8), (1, 3, 5, 7, 9), (2, 4, 5, 8, 9), (1, 2, 5, 6, 8), (1, 5, 6, 8, 9),
    (1, 2, 3, 7, 8), (4, 6, 7, 8, 9), (2, 3, 6, 7, 9), (1, 2, 4, 5, 6), (1, 2, 4, 6, 7),
    (1, 3, 4, 5, 7), (2, 4, 6, 7, 9), (1, 2, 4, 5, 9), (3, 4, 5, 7, 8), (3, 5, 6, 7, 9),
    (1, 2, 3, 8, 9), (2, 3, 4, 8, 9), (2, 3, 4, 6, 9), (2, 3, 5, 7, 8), (2, 3, 5, 6, 7),
    (2, 3, 4, 5, 6), (2, 3, 4, 5, 8), (1, 3, 4, 5, 6), (3, 4, 6, 8, 9), (1, 2, 6, 7, 8),
    (1, 3, 6, 8, 9), (1, 2, 3, 7, 9), (1, 2, 5, 8, 9), (5, 6, 7, 8, 9), (1, 3, 5, 6, 9),
    (1, 4, 6, 7, 8), (4, 5, 7, 8, 9), (1, 2, 4, 7, 9), (1, 4, 5, 7, 9), (1, 3, 4, 7, 8), (1, 3, 4, 6, 8)}

R = ZZ/1009[x_1..x_9]
mons = apply(T, t -> product toList apply(t, i -> x_i))
S =  simplicialComplex mons

IX = ideal S
depth(R/IX) --- 3

T1 = CT^1(0, gens IX);
--T1 = normalMatrix(0, gens IX);
T2 = CT^2(0, gens IX);

b1 = transpose matrix{{1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1}}
b2 = transpose matrix{{0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0}}
T1 = (T1 * b1) | (T1 * b2)
T1 = (T1 * b2)

T1 = T1_{0,1,3,5,7,9,11,13,14,16,18,20} -- baseromligninger = 0
T1 = T1_{0,2,4,6,7,8,10,12,15,17,19,20}  -- samme her
T1  = T1_{0,2,4,5,6,7,9,10,11,13,16,17,18,19,20}

(f1,r1,g1,c1) = versalDeformation(gens IX, T1,T2, HighestOrder => 50, SmartLift => false, PolynomialCheck => false);
(f1,r1,g1,c1) = versalDeformation(gens IX, T1,T2, HighestOrder => 5, SmartLift => false, PolynomialCheck => false);


ideal mingens ideal sum g1
loadPackage "Binomials"
BPD ideal mingens ideal sum g1
TR = ZZ/1009[gens ring sum f1]
BPD sub(ideal sum g1, TR)

apply(oo, dim)

decompose ideal  mingens ideal sum g1
last oo
comp3 = oo
loadPackage "defMethods"

minimalPresentation(comp3 + ideal(t_2-1,t_13-1,t_6-1,t_11-1))
(f1,r1,g1,c1) = versalDeformation(gens IX, T1,T2, HighestOrder => 20);

IXd = ideal sub(transpose sum f1, toList apply(1..21, i-> (t_i => -i)))

mingens(last L + ideal(t_2-1,t_9 - 1,t_3-1,t_18-1,t_15-1,t_4 -2,t_16-3,t_10-4,t_20+504))
mingens ideal sub(mingens ideal sum g1, {t_20 => -504, t_19 => -84, t_18=> 1, t_17 => 337, t_3 => 1, t_18 => 1, t_10 => 4, t_15 => 1,
	t_2 => 1, t_4 => -1/504, t_16 => 3, t_12 => 2, t_13 => -126, t_14 => 1, t_11 => -84, t_9 => 1, t_7 => -504, t_6 => -336, t_5 => -252})
IXd = ideal sub(transpose sum f1, {t_20 => -504, t_19 => -84, t_18=> 1, t_17 => 337, t_3 => 1, t_18 => 1, t_10 => 4, t_15 => 1,
	t_2 => 1, t_4 => -1/504, t_16 => 3, t_12 => 2, t_13 => -126, t_14 => 1, t_11 => -84, t_9 => 1, t_7 => -504, t_6 => -336, t_5 => -252,
	t_21 => 1,t_1 => 1,t_8 => 1})

IXd = ideal transpose mingens IXd

L = decompose ideal sum g1
mingens(last L + ideal(t_2-1,t_9 - 1,t_3-1,t_18-1,t_15-1,t_4 -2,t_16-3,t_10-4))


--- the equilbrium triangulation

S = ZZ/1009[x_0..x_6,X,Y,Z]
Xf = {{0,1,3,4},{0,1,2,3},{1,2,4,5},{1,2,3,4},{2,3,5,6},{2,3,4,5},{3,4,6,0},{3,4,5,6},{4,5,0,1},{4,5,6,0},{5,6,1,2},{5,6,0,1},{6,0,2,3},{6,0,2,1}}
Yf = {{0,1,2,3},{0,2,4,6},{1,2,3,4},{1,3,5,0},{2,3,4,5},{2,4,6,1},{3,4,5,6},{3,5,0,2},{4,5,6,0},{4,6,1,3},{5,6,0,1},{5,0,2,4},{6,0,1,2}, {6,1,3,5}}
Zf = {{0,2,4,6},{0,1,3,4},{1,3,5,0},{1,2,4,5},{2,4,6,1},{2,3,5,6},{3,5,0,2},{3,4,6,0},{4,6,1,3},{4,5,0,1},{5,0,2,4},{5,6,1,2},{6,1,3,5},{6,0,2,3}}
mons = apply(Xf, t -> product toList apply(t, i -> x_i)*X) | apply(Yf, t -> product toList apply(t, i -> x_i)*Y) | apply(Zf, t -> product toList apply(t, i -> x_i)*Z)

K = simplicialComplex mons

IX= ideal K
T1 = CT^1(0, gens IX);
T1 = T1 * transpose matrix {toList apply(1..42, i -> 1)}
T2 = CT^2(0, gens IX);
nnormalMatrix(0, gens IX);

transpose (f1#0 + f1#1)
betti res IX

--- her løfter faktisk den verselle familien!!
(f1,r1,g1,c1) = versalDeformation(gens IX, T1, T2);

IXd = ideal sub(transpose sum f1, toList apply(1..42, i-> (t_i => random(0,S))))
IXd = ideal sub(transpose sum f1, toList apply(1..1, i-> (t_i => random(0,S))))

IXd =  ideal  mingens IXd
hilbertPolynomial IXd
hilbertPolynomial IX

loadPackage "MinimalPrimes"
installMinprimes()
LL = decompose IXd
#oo
--- så generiske komponenten er union av tre stk (
hilbertPolynomial LL#0
betti res  minimalPresentation(LL#0) --- codim 3, så er en Pfaffian resolusjon
use ring sum f1
IXt = ideal sub(transpose sum f1, toList apply(1..42, i-> (t_i => t_1)))

---
loadPackage "NormalToricVarieties"
P = convexHull transpose((transpose makeCone comp3)_{0..9})
normalFan P
F  = oo
V = normalToricVariety F

polytope sum toList apply(0..49, i-> V_i)
fVector P
isFano V


---
-- Bagschi/Datta 10 vertex triang

R = QQ[x_11,x_12,x_13,x_14,x_22,x_23,x_24,x_33,x_34,x_44]
basicFacets = {x_11*x_22*x_33*x_12*x_13,
    x_11*x_22*x_12*x_14*x_34, x_11*x_22*x_14*x_24*x_34,
    x_11*x_22*x_12*x_13*x_24, x_11*x_22*x_13*x_24*x_34}

gens R
alpha = map(R,R,{x_22,x_23,x_12,x_24,x_33,x_13,x_34,x_11,x_14,x_44})
beta  = map(R,R,{x_22,x_12,x_24,x_23,x_11,x_14,x_13,x_44,x_34,x_33})

G = {alpha * alpha * alpha, alpha, alpha * alpha,
    beta, alpha*beta, (alpha*beta)^2,
    beta * alpha * beta, (beta * alpha * beta)^2,
    alpha * alpha *beta, (alpha * alpha *beta)^2,
    beta * alpha * alpha *beta, (beta * alpha * alpha *beta)^2}
G = unique apply(subsets(G,2), P -> P#0 * P#1)

topFaces = unique flatten apply(G, g -> apply(basicFacets, f -> g f))

S = simplicialComplex topFaces
fVector S

betti res ideal S
I = ideal S

T1 = CT^1(0, gens I);
T2 = CT^2(0, gens I);
normalMatrix(0, gens I)

T1 = submatrix'(T1,{0,5,8,11,12,13,21,23,32,33,37,38})
T1 = submatrix'(T1,{0,1,2,4,5,7,8,12,13,15,21,23,27,32})
(f1,r1,g1,c1) = versalDeformation(gens I, T1, T2, HighestOrder => 4);
(f1,r1,g1,c1) = versalDeformation(gens I, T1, T2, HighestOrder => 50);



decompose ideal mingens ideal sum g1
apply(oo, dim)
sub(transpose sum f1, toList apply(1..27, i-> (t_i => -1)))
It = ideal oo
decompose It
CT^1(0, gens It)
loadPackage "MinimalPrimes"
installMinprimes()

LL = decompose ideal mingens ideal sum g1
apply(LL, i -> (dim i)-dim R)

last LL
dim ring ideal sum g1



