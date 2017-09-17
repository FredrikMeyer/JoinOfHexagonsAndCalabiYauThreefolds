restart 
loadPackage "VersalDeformations"
kk = ZZ/11
R = kk[x_1..x_12]

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

generateX2 = () -> (
    K = random(R^16,R^12);
    a = transpose gens gb K;
    b = entries transpose K;
    b = entries a;
    b = apply(0..11, i-> apply(b#i, z -> z*x_(i+1)));
    bb = sum toList b;
    bb1 = bb_{0..7};
    bb2 = bb_{8..15};
    I1 = minors222tensor {{{bb1#0,bb1#1},{bb1#2,bb1#3}},{{bb1#4,bb1#5},{bb1#6,bb1#7}}};
    I2 = minors222tensor {{{bb2#0,bb2#1},{bb2#2,bb2#3}},{{bb2#4,bb2#5},{bb2#6,bb2#7}}};
    I1+I2
    )

IX = generateX2();
time gens gb IX;


apply(1..12, i-> rank source gens gb minimalPresentation(IX + ideal(x_i-1)))
--mingens ideal singularLocus minimalPresentation(IX + ideal(x_4-1))
GG = oo;
dim ideal GG -- 0
degree ideal GG -- 4 
IS = radical ideal GG
IS = oo --ideal(-x_11,-x_10+2,-x_8-1,-x_4,-x_3,-x_2+1,-x_1+2)




T1 = time CT^1(0, gens IX);
rank source T1

CT = time cotangentSheaf Proj(R/IX);


hilbertPolynomial(IX, Projective => false)
---T1 --- 29-  (over Z/11 og Z/13)
--T1 -- 29 igjen over Z/23
--- fikk 29 også over Z/1009....

--- nå prøvde jeg med random koeffisenter uten å endre noe
--- og fikk T^1=29 etter 4.5 dager....