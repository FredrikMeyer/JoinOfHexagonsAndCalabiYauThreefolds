restart
loadPackage "VersalDeformations"

R = QQ[x_1..x_6,z_1..z_6]
M1 = matrix{{0,x_1,x_2},{x_4,0,x_3},{x_5,x_6,0}}
M2 = matrix{{0,z_1,z_2},{z_4,0,z_3},{z_5,z_6,0}}
I = ideal mingens(minors(2,M1) + minors(2,M2))
-- SR-degenerasjon av M

(CT^1(0, gens I))
(CT^2(-2, gens I))

S = QQ[x_1..x_6,y]
M1 = matrix{{0,x_1,x_2},{x_4,0,x_3},{x_5,x_6,0}}
I = ideal mingens minors(2,M1)

T1 = CT^1(0,gens I)
T2 = CT^2(0,gens I)

T1 = T1_{1,3,5,7,9,11}
(f1,r1,g1,c1) = versalDeformation(gens I, T1, T2, HighestOrder => 20);

transpose sum f1
decompose ideal g1

It = ideal sub(sum f1, {t_1 => 1, t_2 => 1, t_5 => 1,t_6 => 2, t_3 => 2, t_4 => 2})
It = ideal sub(sum f1, {t_1 => 1, t_2 => 1, t_5 => 1,t_6 => 1, t_3 => 1, t_4 => 1})

X = Proj(S/It)
HH^1(tangentSheaf X)
dim I

HH^0(sheaf prune (Hom(It/It^2,S^1/It)))
HH^0(sheaf prune (Hom(It/It^2,S/It)))
--- this is the normal sheaf of dP6, it is 46-dimensional

-----


restart
loadPackage "VersalDeformations"
kk = ZZ/1009
R = kk[x_0..x_8]
R' = kk[x_9..x_17]
S = kk[u,v,a,b]

f = map(S,R,{u^2*a^2,u^2*a*b,u^2*b^2,u*v*a^2, u*v*a*b, u*v*b^2,v^2*a^2,v^2*a*b, v^2*b^2})
f' = map(S,R',{u^2*a^2,u^2*a*b,u^2*b^2,u*v*a^2, u*v*a*b, u*v*b^2,v^2*a^2,v^2*a*b, v^2*b^2})
I = ker f
I
-- kan V(I) glattes?
transpose mingens I
T1 = CT^1(gens I)
T2 = CT^2(gens I)
(f1,r1,g1,c1) = versalDeformation(gens I, T1, T2);
I1 = minimalPresentation ideal sub(transpose sum f1, {t_1 => -1}) --<-- glatt
--ideal mingens I1
---
I' = ker f'

T = kk[x_0..x_17]
II = sub(I,T) + sub(I',T)

--singularLocus minimalPresentation(II + ideal(x_0-1))


--CT^2(0, gens II)
altern = sum toList apply(0..17, i-> (-1)^i * (gens T)_i)
IX=  (II + ideal(sum gens T,random(1,T)))
IX=  (II + ideal(sum gens T,altern))

codim minimalPresentation IX

source gens gb minimalPresentation IX

degree II
T1 = CT^1(0, gens II);
T1 = T1 * transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}}
T2 = CT^2(0, gens II);

(f1,r1,g1,c1) = versalDeformation(gens II, T1, T2, HighestOrder => 8, SmartLift => false, PolynomialCheck => false);

eulers(II/II^2)
euler (II/II^2)
viewHelp eulers
dim ring II
X = Proj(T/IX)


