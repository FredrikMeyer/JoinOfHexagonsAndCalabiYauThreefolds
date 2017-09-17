restart
kk = ZZ/101
R = kk[x_1..x_6,x_7,x_8,x_9,x_10,x_11,y_0]
--R = kk[x_1..x_6,y_0]
S = kk[u,v,w,a_1,a_2,a_3,a_4,a_5]

f = map(S,R,{w^4, u*w^3, u^2*v*w, u^2*v^2, u*v^2*w, v*w^3,u*v*w^2,a_1,a_2,a_3,a_4,a_5})
f = map(S,R,{w^4, u*w^3, u^2*v*w, u^2*v^2, u*v^2*w, v*w^3,u*v*w^2})
I = ker f
A = R/I
Y = Proj A
betti res I
C = res I^2
HH^8(sheaf (image (C.dd_6) ** (R/I)))
betti res (I^2)
HH^7(OO_Y(-6))

HH^5(sheaf prune ((ker C.dd_1) **A))
M2 = prune (image (C.dd_3 ** A));
HH^6(sheaf M2)
C = res I
HH^7(OO_Y(-4))
HH^7(OO_Y(-7))
(rank HH^7(OO_Y(-8)))*30
HH^7(OO_Y(-5))
2070-1330
960-740
220-135
2*85
170-143
x-27 = -
8*(rank HH^7(OO_Y(-9)))  - (rank HH^7(OO_Y(-10))) 
HH^5(sheaf module I^2)

C = res I^2
image (g C.dd_1) ** A;
II = sheaf oo
HH^4(II)

g = map(R,R,matrix{gens R} * random(R^12,R^12));
J = ideal mingens(I + g I)
euler sheaf (J/J^2)
euler sheaf module J
euler sheaf module (J^2)
loadPackage "VersalDeformations"
CT^1(gens I)
CT^2(gens I)

dP6 = Proj(R/I)
R2  = kk[x_1..x_6,y_0,y_1]
I2 = sub(I,R2)
Cdp6 = Proj(R2/I2)
apply(0..3, i-> HH^i(tangentSheaf oo))



apply(0..2, i-> HH^i(tangentSheaf dP6))
apply(0..2, i-> HH^i(cotangentSheaf dP6))
apply(0..2, i-> HH^i((cotangentSheaf dP6)(-2)))

hilbertPolynomial(dP6)

transpose gens I

T = kk[t_1..t_9]

M = matrix{{t_1,t_2,t_3},{t_4,t_5,t_6},{t_7,t_8,t_9}}
Ip2p2 = minors(2,M)
B = T/Ip2p2

use R
phi = map(R,T,{x_4,x_5,y_0,y_0,x_6,x_1,x_3,y_0,x_2})
phi1 = map(R,T,{x_4,x_5,y_0,y_0+random(0,R),x_6,x_1,x_3,y_0+random(0,R),x_2})
psi = map(A,B,{x_4,x_5,y_0,y_0,x_6,x_1,x_3,y_0,x_2})
psi1 = map(A,T,{x_4,x_5,y_0,y_0+random(0,R),x_6,x_1,x_3,y_0+random(0,R),x_2})


psi1  = map(A1,B,matrix phi1 ** A1)
transpose mingens ker psi1

C = B/ker psi1

A1 = R/(minors(2,phi1 M))



-------
restart
R = QQ[x_11,x_12,x_13,x_21,x_22,x_23,x_31,x_32,x_33]

M = 2*genericMatrix(R,3,3)
A = 2*matrix{{x_11, (x_12+x_21)/2, (x_13+x_31)/2},
          {(x_12+x_21)/2, x_22, (x_23+x_32)/2},
	  {(x_13+x_31)/2, (x_23+x_32)/2, x_33}}
B = M-A

P  = transpose (transpose(B | A) | transpose (-A | -B))
If = pfaffians(4,P)

S = R[t]
MS = sub(M,S)
AS = sub(A,S)
BS = MS-AS
PS  = transpose (transpose(BS | AS) | transpose (-AS | -BS*t))
If = pfaffians(4,PS)
transpose mingens If

b = transpose matrix{{B_(1,2),B_(2,0),B_(0,1)}}
IS =  (minors(2,A) + ideal(A*b))

transpose mingens IS

dim radical ideal mingens ideal singularLocus IS

X = Proj(R/IS)
I = minors(2,M)
Y = Proj(R/I)
dim Y
hilbertPolynomial(X,Projective => false)
 hilbertPolynomial(Y,Projective => false)

radical ideal mingens ideal singularLocus X
dim oo

transpose mingens IS

gens R
use R
f = map(T,R,{y_11,y_12, y_13, y_0-y_12,y_22,y_23,y_0-y_13,y_0-y_23,y_33})
T = QQ[y_11,y_12,y_13,y_22,y_23,y_33,y_0,y_1]
transpose mingens (f IS)
IdP = homogenize(ideal oo, y_1)

hilbertPolynomial  IdP

dim radical ideal mingens ideal singularLocus IdP
X = T/IdP
use T
IdP2 = sub(IdP, {y_33 => y_33/2, y_22 => y_22/2, y_11 => y_11/2, y_12 => y_12/2, y_13 => y_13/2, y_23 => y_23/2})

decompose IdP2

loadPackage "VersalDeformations"
T1 = CT^1(0, gens IdP2)
T2 = CT^2(0, gens IdP2)

(F,RL,G,C) = versalDeformation(gens IdP2, T1, T2, Verbose => 4, SmartLift => false, HighestOrder => 5);

ideal sum G
decompose tangentCone ideal sum G
ring oo#0
dim T
apply(oo, dim)

C = res I
M = C.dd_3
Y = Proj(A)
(sheaf (ker M ** A)) ** OO_Y(6)
NY = sheaf Hom(I/I^2,A)
sheaf (ker(M **A)) ** OO_Y(6) -- <- normalknippet

sheafHom(NY ** OO_Y(-7),OO_Y)

gens module N
(prune exteriorPower(4, NY)) ** OO_Y(-7)
---
R =QQ[x_0..x_8]
M = genericMatrix(R,3,3)
I = minors(2,M)
codim I
Y = Proj(R/I)
NY = sheaf Hom(I/I^2, R/I)
(prune exteriorPower(4, NY)) ** OO_Y(-9)

f = sum toList apply(0..8, i -> x_i^3)
J = I +f

X = Proj(R/J)

euler X
hh^(1,1)(X)
hh^(1,2)(X)
---
restart
R = QQ[u,v,r,s,a,b]
S = QQ[x_0..x_7]
f = map(R,S,{u*r*a,u*r*b, u*s*a,u*s*b, v*r*a,v*r*b, v*s*a, v*s*b})
I = ker f
factor 162
---

R = QQ[x_0..x_6]
M = matrix{{x_0,x_1,x_2},{x_3,x_0,x_4},{x_5,x_6,x_0}}

I = minors(2,M)

loadPackage "LocalRings"
setMaxIdeal ideal gens R
localMingens gens I
0
