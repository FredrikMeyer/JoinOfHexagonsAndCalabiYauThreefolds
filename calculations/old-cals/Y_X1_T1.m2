restart

R = QQ[x_0..x_17]

M1 = genericMatrix(R,3,3)
M2 = genericMatrix(R,x_9,3,3)

I = minors(2,M1) + minors(2,M2)

B = R/I

IX =  ideal apply(1..6, i-> random(1,B))

A = B/IX

Hom(I,A)

