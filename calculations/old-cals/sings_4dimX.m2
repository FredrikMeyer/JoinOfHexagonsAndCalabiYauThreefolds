restart

R = ZZ/1009[y_0,x_1..x_6,z_1..z_6]

M1 = matrix{{y_0,x_1,x_2},{x_4,y_0,x_3},{x_5,x_6,y_0}}
M2 = matrix{{y_0,z_1,z_2},{z_4,y_0,z_3},{z_5,z_6,y_0}}

I = minors(2,M1)+ minors(2,M2)
sum sum entries vars R
IX = I + sum sum entries vars R

mingens ideal singularLocus IX

mingens ideal singularLocus radical minimalPresentation(I + ideal(x_1-1))

mingens ideal singularLocus I

SING = oo

radical ideal SING

dim ideal SING
SS1 = radical ideal (mingens (radical ideal SING_{2000..3702} + ideal SING + ideal(x_1*x_6*z_6,x_5*x_6*z_6,x_1*z_5*z_6)))_{120..137}
SS2 = ideal mingens (radical ideal SING_{2000..3702} + ideal SING + ideal(x_1*x_6*z_6,x_5*x_6*z_6,x_1*z_5*z_6))
radical ideal mingens(SS2 + SS1)
SINGI = oo
decompose SINGI
LL = oo
#LL
