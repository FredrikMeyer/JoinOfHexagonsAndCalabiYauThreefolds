restart

R = ZZ/101[x_0,x_1,x_2,y_0,y_1,y_2,z_0,z_1,z_2,w_0,w_1,w_2,r,s]

M1 = matrix{{x_0,x_1,x_2}}
M2 = matrix{{y_0,y_1,y_2}}
M3 = matrix{{z_0,z_1,z_2}}
M4 = matrix{{w_0,w_1,w_2}}

N = r* transpose M1 * M2
N' = s*transpose M3 * M4


mons = flatten entries N |  flatten entries N'

L = {{0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 1, 0, -5, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, -5, 0, 0, 0, 0, -5, 0, 0, 0, 0}, {0,
      0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {-1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0}, {0, -5, 0, 0, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}
h1 = L#4
matrix {mons} * transpose matrix{h1}


--I = ideal flatten flatten apply(0..5, i-> entries( matrix {mons} * random(R^18,R^1)))
I = ideal flatten flatten apply(0..5, i-> entries( matrix {mons} * transpose matrix{L#i}))


I = saturate(I,ideal(x_0,x_1,x_2))
I = saturate(I,ideal(y_0,y_1,y_2))    
I = saturate(I,ideal(z_1,z_2,z_2))
I = saturate(I,ideal(w_0,w_1,w_2))
I = saturate(I,ideal(r,s))

minimalPresentation ideal mingens (I +ideal(s-1)+ideal(x_0-1) + ideal(z_1-1)+ideal(y_0-1) + ideal(w_0-1));
mingens oo;
mingens ideal singularLocus ideal oo
Ielim = eliminate(eliminate(eliminate(eliminate(eliminate(eliminate(eliminate(eliminate(I,w_0),z_1),w_2),w_1),z_2),z_0),r),s) + ideal(w_0,w_1,w_2,z_0,z_1,z_2,r,s)


dim ideal mingens ideal singularLocus minimalPresentation Ielim







for i from 0 to 2 do (
    for j from 0 to 2 do (
	for k from 0 to 2 do (
	    for l from 0 to 2 do (
		print degree ideal mingens minimalPresentation(I + ideal(x_i-1) + ideal(y_j-1) + ideal(z_k-1)+ideal(w_l-1)+ideal(r-1))
		);
	    );
	);
    );


ideal mingens minimalPresentation(I + ideal(x_2-1) + ideal(y_0-1) + ideal(z_0-1)+ideal(w_0-1)+ideal(r-1))
degree oo



--Iproj = eliminate(eliminate(eliminate(eliminate(eliminate(eliminate(I,z_0),z_1),z_2),w_0),w_1),w_2) + ideal(z_0..z_2,w_0..w_2)




for i from 0 to 2 do (
    for j from 0 to 2 do (
	for k from 0 to 2 do (
	    
    	    )
	)
    )

minimalPresentation ((IX + ideal(x_2-1) + ideal(z_1-1) + ideal(w_2-1)+ideal(r-1) + ideal(y_1-1)))
mingens ideal singularLocus ideal mingens oo;
dim ideal oo

ideal mingens (IX + ideal(r))




loadPackage "MinimalPrimes"
installMinprimes()
decompose IX


--------------------
restart

basisE = set apply({matrix{{1,0,0}},matrix{{0,1,0}},matrix{{0,0,1}}},transpose)
basisEE = apply(toList(basisE ** basisE), S -> S#0 * transpose S#1)
e1 = (toList basisE)#2
e2 = (toList basisE)#1
e3 = (toList basisE)#0

M_1 = e1 * transpose e2 -- + (e3 * transpose e3)
M_2 = e1 * transpose e3 --+ (e1 * transpose e1)
M_3 = e2 * transpose e3 --+ (e2 * transpose e2)
M_4 = e2 * transpose e1 --+ (e3 * transpose e3)
M_5 = e3 * transpose e1 --+ (e1 * transpose e1)
M_6 = e3 * transpose e2 --+ (e2 * transpose e2)
D_1 = e3 * transpose e3
D_2 = e2 * transpose e2
D_3 = e1 * transpose e1
D_4 = e3 * transpose e3
D_5 = e2 * transpose e2
D_6 = e1 * transpose e1
Z = 0 * id_(ZZ^3)
blocks = apply(1..6, i-> (M_i | 5*D_i)) |apply(1..6, i-> D_i| M_i)



entries transpose (matrix (gens ker matrix toList apply(blocks, M -> flatten entries M)))
toString oo
------------
------------- projeksjoner


restart
L = {{0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 1, 0, -5, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, -5, 0, 0, 0, 0, -5, 0, 0, 0, 0}, {0,
      0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {-1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0}, {0, -5, 0, 0, 0, 0, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}
loadPackage "RationalPoints"
kk = ZZ/109
--kk = QQ
R = kk[x_1..x_18]
A = matrix{{x_1,x_2,x_3},{x_4,x_5,x_6},{x_7,x_8,x_9}}
B = matrix{{x_10,x_11,x_12},{x_13,x_14,x_15},{x_16,x_17,x_18}}

I= (minors(2,A) + minors(2,B))
mons = gens R

IH =  ideal apply(L, l -> sum toList apply(1..18, i-> x_i * l#(i-1)))
IHr = ideal( random(1,R),random(1,R) ,random(1,R) ,random(1,R) , random(1,R),random(1,R))

IX = I + IHr
IX = I + IH

IXmin = minimalPresentation IX
use ring IXmin

IXP = minimalPresentation (eliminate(toList (x_10..x_18), IX) + ideal(x_10..x_18))
omegaX = Ext^(codim IX)((ring IX)^1/IX,(ring IX)^{-19})
prune sheaf omegaX
pruned = oo
apply(0..4, i-> HH^i(pruned))



dualModule = Hom( omegaX, (ring IXP)^1/IXP)
f = homomorphism dualModule_{0}
canGens = f*basis(0,omegaX)
ringX = (ring IXP)/IXP
P36 = kk[w_0..w_35]
idealXcan = trim kernel map(ringX, P36,substitute(matrix canGens,ringX))




betti res IXP
II = ideal (IXP_*)_{0..8}
f =  last (IXP_*)


CT^2(0, gens IXP)

X = Proj(ring IXP/IXP)
HH^3(OO_X) --- 100!!!

use ring IXP
Iloc = ideal transpose mingens minimalPresentation(IXP + ideal(x_1-1))

L = ((ideal mingens ideal jacobian Iloc) + Iloc)_*;

#rationalPoints ideal L
gammap = oo/11^4
apply(1..50, d -> (log(d)-log(gammap))/log(11))
sum toList oo
J = jacobian ideal L;
tally apply(8000, i-> (
	point := random(kk^1,kk^4);
	if sub(ideal L,point) == 0 then
	    rank sub(J,point)
	))

sub(59*11^2/8000,RR)
sub(7*11/8000,RR)

dim ideal L
radical ideal L
degree oo

torus = ideal apply(flatten apply(apply(apply(flatten entries gens IXP, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)

loadPackage "Binomials"
BPD torus
LL = oo
select(LL, I -> dim I == 1)




radical ideal L

--#rationalPoints ideal L
--#rationalPoints Iloc

C = res IXP
betti C
ideal transpose C.dd_5 == ideal C.dd_1 --- gorenstein ?? 


(IXP_*)_{0..8}
f = last IXP_*
I1 = ideal ((IXP_*)_{0..8})

hilbertPolynomial IXP
hilbertPolynomial(IXP,Projective => false)
S = ring IXP
XP = Proj (S/IXP)

    
HH^3(OO_XP)
apply(0..3, i-> HH^i(OO_XP))
viewHelp Ext
A = S/IXP
II = prune sheaf (module IXP ** A)
Ext^2(II, OO_XP)
Hom(II,OO_XP) -- 847!!!
HH^3(II)

HH^0(sheaf Hom(prune ((module IXP) ** A),

use ring IXP
I9 =  minimalPresentation(IXP + ideal(x_9-1))
singi = ideal mingens (I9 + ideal jacobian I9)

degree radical singi


loadPackage "MinimalPrimes"
installMinprimes()



minprimes(IXP, Strategy => "NoBirational")


loadPackage "VersalDeformations"
CT^1(0, gens IXP)

----
restart
R = QQ[x_1..x_9,y_1..y_9]
M1 = genericMatrix(R,3,3)
M2 = genericMatrix(R,y_1,3,3)
I = minors(2,M1) + minors(2,M2) + minors(2,M1+M2)
I = saturate(I,ideal(x_1..x_9)*ideal(y_1..y_9))

loadPackage "MinimalPrimes"
installMinprimes()

I1 = first decompose I
I2 = last decompose I

degree first decompose(I1+I2)

v = flatten entries( random(R^3,R^1) * random(R^1,R^3))
v = flatten entries matrix{{1,0,0},{0_R,0,0},{0,0,0}}


 ideal mingens sub(I, toList apply(1..9, i -> x_i => v#(i-1))) + ideal(x_1..x_9)


 

minimalPresentation oo
decompose oo

--eliminate(toList (y_1..y_9),I1) + ideal(y_1..y_9)
Z = 0 * id_(ZZ^3)
ms = toList apply(1..12, i-> ((random(R^3,R^1)*random(R^1,R^3) | Z) || (Z | random(R^3,R^1)*random(R^1,R^3))))
Pt = ms#0
A = Pt_{0..2}^{0..2}
M = matrix toList apply(0..11, i-> flatten entries ms#i)

AX = (A | Z) || (Z | matrix{{x_1,x_2,x_3},{x_4,x_5,x_6},{x_7,x_8,x_9}})
MX = M || matrix{flatten entries AX}
MX' = submatrix'(MX,{3,4,5,9,10,11,15,16,17,18,19,20,24,25,26,30,31,32})

Imins = ideal mingens minors(13,MX');
fiber = Imins + minors(2,matrix{{x_1,x_2,x_3},{x_4,x_5,x_6},{x_7,x_8,x_9}}) + ideal(x_10..x_18)
radical fiber


A = R/(radical fiber)
sub(matrix{{x_1,x_2,x_3},{x_4,x_5,x_6},{x_7,x_8,x_9}},A)


