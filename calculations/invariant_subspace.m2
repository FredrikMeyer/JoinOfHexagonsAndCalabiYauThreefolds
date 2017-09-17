restart

R = QQ[w]/ideal(w^2+w+1)

M = random(R^3,R^3)

r = sub(transpose matrix{{0,1,0},{0,0,1},{1,0,0}},R)
s = transpose matrix{{w,0,0},{0,w*w,0},{0,0,1}}

S1 = set{r,r*r,r*r*r}
S2 = set{s,s*s,s*s*s}

G = apply(toList(S1 ** S2), S -> S#0 * S#1)

--- GG
Z = 0 * r
rr = (r | Z) || (Z | r)
ss = (s | Z) || (Z | s)
SS1 = set{rr,rr*rr,rr*rr*rr}
SS2 = set{ss,ss*ss,ss*ss*ss}

GG = apply(toList(SS1 ** SS2), S -> S#0 * S#1)

--------
M = (random(R^3,R^3) | Z) || (Z | random(R^3,R^3))
sum apply(toList SS1, g -> g*M*g)
-----------
(sum apply(G, g -> g*M*g*g)) * (1/9)

G*G*M*G*w

A1 = matrix{{0,1,0},{0,0,0},{0,0,0}}
A2 = matrix{{0,0,1},{0,0,0},{0,0,0}}
A3 = matrix{{0,0,0},{0,0,1},{0,0,0}}
A4 = matrix{{0,0,0},{1,0,0},{0,0,0}}
A5 = matrix{{0,0,0},{0,0,0},{1,0,0}}
A6 = matrix{{0,0,0},{0,0,0},{0,1,0}}
D = matrix{{1,0,0},{0,1,0},{0,0,1}}

AA = sum{A1,A2,A3,A4,A5,A6,D}

unique apply(G,g -> g *A*g*g)
A1+D

Z = 0*A1


r*A1*r*r
M1 = transpose matrix toList apply(unique apply(G,g -> g *(A1+D)*g*g), N -> flatten entries N)
M2 = transpose matrix toList apply(unique apply(G,g -> g *(A2+D)*g*g), N -> flatten entries N)
rank (M1 | M2)




A = matrix{{1,0,0},{0,0,1},{0,1,0}}
G#0 * A * G#0
A = random(R^3,R^3)

matrix{flatten entries M}
sum apply(toList S2 ,g -> g *A*g*g)

unique oo

D1 = (D | Z) || (Z | Z)
D2 = (Z | Z) || (Z | D)

AAi = apply({A1,A2,A3,A4,A5,A6}, a -> (a | Z) || (Z | D))
AAi' = apply({A1,A2,A3,A4,A5,A6}, a -> (D | Z) || (Z | a))

AA = AAi | AAi'


baner = apply(AA,aa ->  matrix toList apply(apply(GG,g -> g *(aa)*g*g), N -> flatten entries N))
MM  = baner#0
apply(baner, b -> MM = MM | b)



