restart

kk = ZZ/101
R = kk[x_0..x_9]

I = pfaffians(4,genericSkewMatrix(R,5))

res (I^2)

f = map(R,R,transpose (random(R^10,R^10) * transpose matrix{gens R}))
J = I + f I

C  =res J
betti C
Y = Proj(R/J)
X = Proj(R/I)
HH^6(OO_X(-6))
HH^6(OO_X(-5))
betti res (I/I^2)
C = res I

A = R/I

JA = sub(J,A)
sheaf (JA/JA^2)
HH^3(oo)



Y0 = Proj(R/ideal leadTerm J)
betti res (I/I^2)
HH^0(OO_Y(2))
HH^0(OO_Y(3))
HH^0(OO_Y(5))
PP9 = Proj R
HH^6(OO_X(-6))
HH^6(OO_X(-5))

C =  res (I/I^2)
betti C
HH^0 sheaf  Hom(I/I^2,R/I)

