restart

R = QQ[x_0..x_6,t]

M = matrix{
    {0,   0,  0, 0, x_1, x_2},
    {0,   0,  0, x_4, 0, x_3},
    {0,   0,   0, x_5,  x_6,0},
    {0,   0,   0,   0,  0,  0},
    {0,   0,   0,   0,   0,   0},
    {0,0,0 ,0, 0, 0}}
N = M - transpose(M)

M = matrix{
    {0,   (x_0+x_1+x_2)*t,  (x_0+x_3+x_4)*t, x_0, x_1, x_2},
    {0,   0,  (x_5+x_6)*t^2, x_4, x_0, x_3},
    {0,   0,   0, x_5,  x_6,x_0},
    {0,   0,   0,   0,  0,  0},
    {0,   0,   0,   0,   0,   0},
    {0,0,0 ,0, 0, 0}}
N = M - transpose(M)

pfaffians(4,N)

eliminate(pfaffians(4,N),t)

decompose oo



restart
R = QQ[y,x_1..x_6]
M = matrix{{ y,x_1,x_2},
    {x_4,y,x_3},{x_5,x_6,y}}

I = minors(2,M)

loadPackage "VersalDeformations"

CT^1(-1,gens I)
CT^2(gens I)
B = R/I
Hom(I/I^2,B)

use R
J = ideal mingens sub(I,y => 0)
T1=CT^1(0,gens J)
T1n = T1_{0,2,4,6,8,10}
T2=CT^2(0, gens J)
(f1,r1,g1,c1) = versalDeformation(gens J, T1,T2);

transpose sum f1
T1

