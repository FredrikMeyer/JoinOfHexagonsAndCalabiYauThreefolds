restart

kk = ZZ/3001
R = kk[x_1..x_18]
M1 = genericMatrix(R,3,3)
M2 = genericMatrix(R,x_10,3,3)

I = minors(2,M1) + minors(2,M2)
ran =  random(R^6,R^18)
randomforms = (transpose gens gb transpose ran ) * (transpose matrix {gens R})

A = R/I
IH = sub(ideal randomforms,A)
B = A/IH

C = res(I, LengthLimit => 2)

N = Hom(image C.dd_1, B)
dI = sub(transpose jacobian C.dd_1,B)
T1 = N/image dI
T1 == 0 --- faktisk null!!
prune T1
prune N


Y = Proj(A)
TY = tangentSheaf Y;
HH^0(TY)
X = Proj(B)


apply(0..6, i-> HH^0(TY(-i)))
for i from 0 to 6 do (
    for j from 0 to 6 do (
	print (i,j);
	print time  HH^j(TY(-i));
	)
    )
--- konklusjon H^i(TY(-j))=0 om (i,j) != (0,0).

HH^0(sheaf((module TY)**B))
