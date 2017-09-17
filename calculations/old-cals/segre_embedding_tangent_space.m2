restart

S = QQ[a_1..a_9]
R = S[x_1..x_9]
M = genericMatrix(R,3,3)
I = minors(2,M)


h = sum toList apply(1..9, i -> x_i * a_i)

transpose mingens (minors(5,jacobian(I+h)) + (I+h))
M
jacobian I





T = QQ[a_1..a_9,x_1..x_9]
IT = sub(transpose mingens (minors(5,jacobian(I+h)) + (I+h)),T)
eliminate(x_4,eliminate(x_3,eliminate(x_2,eliminate(x_1,ideal IT))))
eliminate(x_5,oo)
eliminate(x_6,oo)
eliminate(x_7,oo)
eliminate(x_8,oo)
eliminate(x_9,oo)
