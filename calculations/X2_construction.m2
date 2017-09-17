restart 
R = ZZ/7[x_1..x_12]

L = {{{x_1,x_2},{x_3,x_4}},{{x_5,x_6},{x_7,x_8}}}

minors222tensor = (L) -> ( -- L is a list of lists of lists
    eqs = {L#0#0#0*L#1#0#1-L#0#0#1*L#1#0#0,
	L#1#0#0*L#1#1#1-L#1#1#0*L#1#0#1,
	L#1#1#0*L#0#1#1-L#1#1#1*L#0#1#0,
	L#0#1#0*L#0#0#1 - L#0#1#1*L#0#0#0,
	L#1#0#1*L#0#1#1-L#1#1#1*L#0#0#1,
	L#1#0#0*L#0#1#0-L#1#1#0*L#0#0#0};
    eqs = eqs | {L#0#0#0 * L#1#1#1 - L#0#0#1*L#1#1#0,
	L#1#0#0*L#0#1#1-L#1#0#1*L#0#1#0,
	L#0#0#1*L#1#1#0 - L#1#0#1*L#0#1#0};
    ideal eqs
    )

randomTensor = () -> (
    toList apply(0..1, j -> toList apply(0..1, k -> toList (random(0,R),random(0,R))))
    )

generateX2 = () -> (
    ts = apply(1..12, i-> randomTensor());
    ts' = apply(1..12, i-> randomTensor());
    N = transpose matrix toList apply(ts, t -> flatten flatten flatten t);
    N' = transpose matrix toList apply(ts', t -> flatten flatten flatten t);
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
    bb1 = bb_{0..7};
    bb2 = bb_{8..15};
    I1 = minors222tensor {{{bb1#0,bb1#1},{bb1#2,bb1#3}},{{bb1#4,bb1#5},{bb1#6,bb1#7}}};
    I2 = minors222tensor {{{bb2#0,bb2#1},{bb2#2,bb2#3}},{{bb2#4,bb2#5},{bb2#6,bb2#7}}};
    I1+I2
    )

IX = generateX2()

time gens gb IX;


--mingens ideal singularLocus minimalPresentation(IX + ideal(x_11-1))

loadPackage "VersalDeformations"

T1 = time CT^1(0, gens IX);

hilbertPolynomial(IX, Projective => false)
---T1 --- 29-  (over Z/11 og Z/13)
