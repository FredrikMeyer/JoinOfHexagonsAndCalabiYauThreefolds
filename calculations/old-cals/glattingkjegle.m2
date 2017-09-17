restart
R = QQ[x_1..x_6,y_0..y_1,s,t]
M = matrix{{x_1,y_0,x_6},{x_2,x_3,y_0+t*y_1},{y_0+s*y_1,x_4,x_5}}
Ifam = minors(2,M)

S = QQ[x_1..x_6,y_0..y_1]
Ispes = sub(sub(Ifam, {s => 2, t => 1}),S)
dim Proj(S/Ispes)

transpose gens Ispes

SS = QQ[x_1..x_6,y_0,y_1,y_2]
M = matrix{{x_1,y_0,x_6},{x_2,x_3,y_1},{y_2,x_4,x_5}}
I =  minors(2,M) + ideal(y_0+y_2-2*y_1) --<- ideal PP(T_P2)

J = I + ideal(y_0-y_1)
Y = Proj (SS/J)

T = QQ[a,b,c,p,q,r]
I = ideal(p*b+a*r-2*q*c, p*b-q*c)
saturate(saturate(I, ideal(p,q,r)),ideal(a,b,c))
saturate((ideal(c) + I), ideal(a,b,c)*ideal(p,q,r))
decompose 
saturate((ideal(p) + I), ideal(a,b,c)*ideal(p,q,r))
decompose oo
