restart
loadPackage "VersalDeformations"

kk = ZZ/1009
S = kk[x_1..x_18]
M1 = genericMatrix(S,3,3)
M2 = genericMatrix(S,x_10,3,3)
I  = minors(2,M1) + minors(2,M2)
H  = ideal(x_1+x_3+x_10+x_12,
    x_2+x_6+x_11+x_17,
    x_4+x_7+x_13+x_15,
    x_6+x_5+x_17+x_18,
    x_1+x_9+x_10+x_14,
    x_3+x_7+x_8+x_12+x_16+x_15+x_9)

--singY = intersect (ideal(x_1..x_9) + minors(2,M2),ideal(x_10..x_18) + minors(2,M1))
-- H snitter ikke singloc(I).
IX = I+H
rank source gens gb IX -- (bare 27!!)
--CT^1(0, gens IX)
--mingens ideal singularLocus  minimalPresentation(IX + ideal(x_1-1))


B = S/I
A = B/sub(H,B)
Y = Proj(B)
CTY = time prune cotangentSheaf Y;
--- 190 sek (på MB Pro)
X = Proj(B/sub(H,B))

--euler sheaf ((module CTY) ** A) (too long)

apply(0..6, i-> euler  CTY(-i))
l = oo
apply(0..6, i-> (-1)^i * binomial(6,i)*l#i)
sum toList oo
---HH^0(OO_X)




