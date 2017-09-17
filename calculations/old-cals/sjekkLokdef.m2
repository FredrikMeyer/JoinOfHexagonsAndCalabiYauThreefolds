restart
kk = QQ
T = kk[x_2,x_6]
k_1 = random(2,T)+1
k_2 = random(2,T)+1
R = kk[z_1..z_6,y_1,x_2,x_6]
k1 = sub(k_1,R)
k2 = sub(k_2,R)

k1 = sub(hs',R)
k2 = sub(hs, R)
M = matrix{{y_1,z_1,z_2},{z_4,y_1+k1,z_3},{z_5,z_6,y_1+k2}}
I = minors(2,M)

S = QQ[r_1..r_6,Y_1,Y_2,Y_3]
M2 = matrix{{Y_1,r_1,r_2},{r_4,Y_2,r_3},{r_5,r_6,Y_3}}
J = minors(2,M2)

A = R/I
B = S/J
f = map(A,B,{z_1,z_2,z_3,z_4,z_5,z_6,y_1,y_1+k1,y_1+k2})

ker f
mingens ideal singularLocus I

decompose oo
apply(oo, degree)
apply(ooo,dim)
sings = o15
sings#0
transpose mingens sings#0 | transpose mingens sings#1
k1
k2



S = QQ[y_0,x_2..x_6,y_1,y_2,z_1..z_6]
avb = matrix {{y_0, x_2, x_2*y_1, y_0*y_1, y_0*x_6, x_6, y_1, x_2*x_6}}

PP2 = matrix{{y_0,1,x_2},{x_4,y_1,x_3},{x_5,x_6,y_2}}
IPP2 =  minors(2,PP2)
minimalPresentation IPP2
f = IPP2.cache.minimalPresentationMap

SSx = QQ[x_1..x_6]
h1x = random(1,SSx)
h2x = random(1,SSx)

SSz = QQ[z_1..z_6]
h1z = random(1,SSz)
h2z = random(1,SSz)

H1x = sub(h1x,S)+1
H2x = sub(h2x,S)+1
H1z = sub(h1z,S)+1
H2z = sub(h2z,S)+1


fH1x = f H1x
fH2x = f H2x
fH1z = f H1z
fH2z = f H2z
hs = sub(fH1x, y_0 => y_1+fH1z)
hs' = sub(fH2x, y_0 => y_1+fH1z)
