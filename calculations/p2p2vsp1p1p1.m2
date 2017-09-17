restart
R = QQ[x_0..x_7, Weights => {1,2,2,1,2,2,2,2}]
S = QQ[r,s,u,v,a,b]

M = matrix{{x_0,x_1,x_2},{x_3,x_4,x_5},{x_6,x_7,x_0+x_3+x_5+x_7}}
I1 = minors(2,M)

f = map(S,R,{r*u*a, r*u*b,r*v*a,r*v*b, s*u*a,s*u*b,s*v*a,s*v*b})
I2 = ker f

N1 = sheaf prune Hom(I1/I1^2,R^1/I1);
HH^0(N1) -- 55
N2 = sheaf prune Hom(I2/I2^2,R^1/I2);
HH^0(N2) -- 54


Y = Proj(R/I2)
9*9-1
2*55-80
-36=x-30
-6

--apply(0..7, i-> dim singularLocus minimalPresentation(I1 + (x_i-1)))

loadPackage "VersalDeformations"

IL1 = ideal mingens ideal leadTerm I1
IL2 = ideal leadTerm I2
NL1 = sheaf prune Hom(IL1/IL1^2,R^1/IL1);
NL2 = sheaf prune Hom(IL2/IL2^2,R^1/IL2);
HH^0(NL2)
HH^0(NL1)

T11 = CT^2(0, gens IL1)
T12 = CT^2(0, gens IL2)

X1 = Proj(R/I1)
X2 = Proj(R/I2)

HH^0(tangentSheaf X2) -- 9
HH^0(tangentSheaf X1) -- 8

leadTerm I1
(ideal leadTerm I1)_*
(I1)_*
I2


loadPackage "SimplicialComplexes"
S1 = simplicialComplex monomialIdeal leadTerm I1
S2 = simplicialComplex monomialIdeal leadTerm I2

saturate(ideal faces(3, S1), ideal(x_0*x_7))

faces(3,S2)
(f1,r1,g1,c1) = versalDeformation(gens ideal leadTerm I2,CT^1(0, gens ideal leadTerm I2), CT^2(0, gens ideal leadTerm I2),HighestOrder => 2);
apply(decompose ideal sum g1, i -> dim i - dim R)

---
-- dp6
---
restart
R = QQ[x_0..x_6]
