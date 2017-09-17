restart

kk = ZZ/109
P1ring = kk[u,v,a,b,r,s,u',v',a',b',r',s']
S = kk[x_0..x_15]    
M = matrix {{u*a*r, u*a*s, u*b*r, u*b*s, v*a*r, v*a*s, v*b*r, v*b*s}}
M' = matrix {{u'*a'*r', u'*a'*s', u'*b'*r', u'*b'*s', v'*a'*r', v'*a'*s', v'*b'*r', v'*b'*s'}}
I = ker map(P1ring,S,M | M')


--betti res I
IH = ideal apply(1..4, i-> random(1,S))
IH = ideal (x_0+x_6+x_8+x_14,x_4+x_7+x_12+x_15,x_5+x_3+x_13+x_11,x_1+x_2+x_6+x_9+x_10+x_14)
IX = I + IH
dim IX

singloc = ideal mingens( I + ideal(x_8..x_15));
singloc = singloc * ideal mingens( I + ideal(x_0..x_7));
dim (IH + singloc)
dim singloc

X = Proj(S/IX)
B = S/I
A = S/IX

HH^2(sheaf prune (((module I) ** B) ** A))

IIX = sheaf prune (IX/IX^2);
betti res (IX/IX^2)


--HH^2(IIX)
HH^2( sheaf ((module I) ** A)) 

betti res I
betti res (I/I^2)
res (module I ** B)
IX/IX^2

minpr = minimalPresentation IX
tally apply(10000000, i-> (
	point := random(kk^1,kk^12);
	if (sub(minpr, point) == 0) then 
	    point else 0
	))
J = jacobian minpr
pts = select(elements TT, s -> s != 0)
tally apply(pts, p -> rank sub(J, p))


loadPackage "MinimalPrimes"
installMinprimes()
apply(decompose (IH + singloc), I -> gens ring minimalPresentation I)
source CT^1(0, gens IX)

--CT^2(0, gens IX)
NX = prune Hom(prune (IX/IX^2),S^1/IX)
HH^0(sheaf NX)
HH^1(sheaf NX)


mingens ideal singularLocus minimalPresentation(IX + ideal(x_0-1))
CT^1(minimalPresentation(IX + ideal(x_0-1)))

