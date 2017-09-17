restart
{*
    Experiments with the projection from the blowup of M.
    *}


kk = ZZ/13
R = kk[x_0..x_17]
M1 = genericMatrix(R,3,3)
M2 = genericMatrix(R,x_9,3,3)

I = minors(2,M1)+minors(2,M2)

T = kk[x_0..x_17,y_0..y_8,z_0..z_8]
N1 = matrix{{y_0,y_1,y_2,y_3,y_4,y_5,y_6,y_7,y_8},
    {x_0,x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8}}
N2 = matrix{{z_0,z_1,z_2,z_3,z_4,z_5,z_6,z_7,z_8},
    {x_9,x_10,x_11,x_12,x_13,x_14,x_15,x_16,x_17}}

J = minors(2,N1) + minors(2,N2)


IB = saturate(sub(I,T) + J, ideal(x_0..x_17))
IB = saturate(saturate(IB, ideal(y_0..y_8)), ideal(z_0..z_8))
IB = sub(I,T) + J 

S = IB_*
IH = ideal apply(1..6, i-> sub(random(1,R),T))
J + IH

use ring J
eliminate(toList first {x_0..x_17},IH + J)



