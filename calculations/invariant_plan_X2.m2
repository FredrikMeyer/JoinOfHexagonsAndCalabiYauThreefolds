restart
kk  = ZZ/3001
R = kk[x_1..x_12]


pars = {2,3,5}
--pars = {1,1,1}
--pars = {2,2,2}
L1 = {{{pars#0 * x_7 + pars#1 * x_8 + pars#2 * x_10,x_1},{x_2,x_3}},{{x_4,x_5},{x_6,pars#2 *x_9 + pars#1 * x_11 + pars#0 * x_12}}}
L2 = {{{pars#0 * x_1 + pars#1 * x_2 + pars#2 * x_4,x_7},{x_8,x_9}},{{x_10,x_11},{x_12,pars#2 *x_3 + pars#1 * x_5 + pars#0 * x_6}}}

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


IX = (minors222tensor L1) + (minors222tensor L2)
time gens gb IX;


--- finne fikspunkter
S  = R[lambda]
M1 =  matrix{{x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10,x_11,x_12}}
M2 =  matrix{{x_1*lambda,x_2*lambda,-x_3*lambda,x_4*lambda,-x_5*lambda,-x_6*lambda,-lambda*x_7,-x_8*lambda,lambda*x_9,-lambda*x_10,lambda*x_11,lambda*x_12}}

Ifiks = saturate(ideal (M1-M2), sub(ideal gens R,S))
decompose(Ifiks + IX)
-- konklusjon: 12 fikspunkter!! de er de eneste

-- automorphisms
wu   = map(R,R, {x_1,x_2,-x_3,x_4,-x_5,-x_6,-x_7,-x_8,x_9,-x_10,x_11,x_12})
sigma = map(R,R, {x_6,x_5,x_4,x_3,x_2,x_1,x_12,x_11,x_10,x_9,x_8,x_7}) -- e_i -> e_i+1
tau   = map(R,R,{x_7,x_8,x_9,x_10,x_11,x_12,x_1,x_2,x_3,x_4,x_5,x_6}) -- bytter om E*E*E-faktorene
sigma IX == IX
tau IX == IX
wu IX == IX


loadPackage "MinimalPrimes"
installMinprimes()

(-48 + 2*12 + 3*x)/2 = 48?

preims = {}
for i from 1 to 3  do(
    IXi = IX + (x_i-1);
    singi = time ideal singularLocus ideal mingens minimalPresentation IXi;
    singi       = time radical ideal mingens singi;
    invi         = time IXi.cache.minimalPresentationMap;
    preimi = time saturate homogenize(preimage(invi,singi),x_i);
    preims = preims | {preimi};
    print(i);
    )

unique (preims | apply(preims, I -> sigma I) | apply(preims, I -> tau I) | apply(preims, I -> sigma tau I))
#oo --- får alle singularietetene !! (koordinatpunktene)

decompose intersect oo


IX1 = IX + (x_1-1)
sing1 = time ideal singularLocus ideal mingens minimalPresentation IX1;
sing1       = time radical ideal mingens sing1;
inv1         = time IX1.cache.minimalPresentationMap;
preim1 = time saturate homogenize(preimage(inv1,sing1),x_1)

IX2 = IX + (x_2-1)
sing2 = time ideal singularLocus ideal mingens minimalPresentation IX2;
sing2       = time radical ideal mingens sing2;
inv2         = time IX2.cache.minimalPresentationMap;
preim2 = time saturate homogenize(preimage(inv2,sing2),x_2)
IX11 = IX + (x_11-1)
sing11 = time ideal singularLocus ideal mingens minimalPresentation IX11;
sing11      = time radical ideal mingens sing11;
inv11         = time IX11.cache.minimalPresentationMap;
preim11 = time saturate homogenize(preimage(inv11,sing11),x_11)
IX12 = IX + (x_12-1)
sing12 = time ideal singularLocus ideal mingens minimalPresentation IX12;
sing12      = time radical ideal mingens sing12;
inv12         = time IX12.cache.minimalPresentationMap;
preim12 = time saturate homogenize(preimage(inv12,sing12),x_12)

 preim1

degree sing1
decompose sing1

IX11 = minimalPresentation(IX + ideal(x_1-1))
loadPackage "LocalRings"



---
IX1 = IX + ideal(x_1-1)
IX11 = minimalPresentation IX1
setMaxIdeal ideal gens ring IX11

S = flatten entries localMingens gens IX11
S#3
use ring oo
S#3 % x_12
f = numerator((S#3 - (S#3 % x_12))/x_12)
S' = apply(S, g -> numerator sub(g, x_12 => -x_8*x_10/f))

S'#2
S'#2 % x_11
f = numerator((S'#2 - (S'#2 % x_11))/x_11)
S'' = apply(S', g -> numerator sub(g, x_11 => x_7*x_10/f))
S3 = flatten entries localMingens matrix{{S''#0, S''#1}}

F = gens minimalPresentation ideal gens ((ideal S3) + ideal(x_11,x_12))

loadPackage "VersalDeformations"
setMaxIdeal ideal gens ring F
localMingens gens ideal F
C = localResolution ideal F
A = ring F/ideal F
N = Hom(image C.dd_1, A)
dI = sub(transpose jacobian C.dd_1, A)
T1 = N/image dI
setMaxIdeal ideal gens ring T1
localPrune T1

minimalPresentation tangentCone ideal S3
(+24 + 3 * 24)/2
---  gnaske overbevist om at (1:0:...) er dobbelpunkt

IX2 = IX + ideal(x_2-1)
IX22 = minimalPresentation IX2
setMaxIdeal ideal gens ring IX22

S = flatten entries localMingens gens IX22
use ring S#0
g = S#3 % x_12
f =numerator( (S#3 - g)/x_12)
S1 = apply(S, h -> numerator sub(h, x_12 => -g/f))


g = S1#2 % x_11
f = numerator((S1#2 - g)/x_11)
S2 = apply(S1, h -> numerator sub(h, x_11 => -g/f))


numerator sub(S2#1, x_3 => -(S2#0 % x_6)/(numerator((S2#0 - (S2#0 % x_6))/x_6)))
numerator sub(S2#0, x_3 => -(S2#0 % x_6)/(numerator((S2#0 - (S2#0 % x_6))/x_6)))

minimalPresentation tangentCone (ideal S2 + ideal(x_11,x_12))
S2#0
S2#1

flatten entries localMingens matrix{S2}
oo#0


--

IX1 = minimalPresentation(IX + ideal(x_1-1))

IX1 + random(1,ring IX1) + random(1,ring IX1) + random(1,ring IX1)
degree(oo : saturate(oo))
 
