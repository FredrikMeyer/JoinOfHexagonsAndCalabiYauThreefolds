restart
loadPackage "VersalDeformations"
kk = ZZ/1009
R = kk[x_1..x_12]

generateX2 = () -> (
    ts = apply(1..12, i-> entries random(R^3,R^3));
    ts' = apply(1..12, i-> entries random(R^3,R^3));
    N = transpose matrix toList apply(ts, t -> flatten  flatten t);
    N' = transpose matrix toList apply(ts', t -> flatten flatten t);
    K = N' || N;
    rotateMatrix = (M) -> (
    	r := rank source M;
    	c := rank target M;
    	return matrix table(r,c, (i,j) -> M_(c-j-1,r-i-1));
    	);
    a = transpose gens gb rotateMatrix transpose K;
    b = entries a;
    b = apply(0..11, i-> apply(b#i, z -> z*x_(i+1)));
    bb = sum toList b;
    bb1 = bb_{0..8};
    bb2 = bb_{9..17};
    I1 = minors(2, matrix toList apply(0..2, i-> toList apply(0..2, j-> bb2#(3*i+j))));
    I2 = minors(2, matrix toList apply(0..2, i-> toList apply(0..2, j-> bb1#(3*i+j))));
    I1+I2
    )

IX = generateX2()
time gens gb IX;

loadPackage "VersalDeformations";
T1 = time CT^1(0, gens IX);
source T1


loadPackage "MinimalPrimes"
installMinprimes()

decompose ideal mingens (IX + ideal(x_1))
last subsets(1..12,3)
apply(subsets(1..12, 3), S -> dim(IX + ideal(x_(S#0),x_(S#1),x_(S#2))))
dim ideal mingens(IX + ideal(x_1,x_2,x_))

decompose minimalPresentation ideal mingens(IX + ideal(x_10,x_11,x_12))

