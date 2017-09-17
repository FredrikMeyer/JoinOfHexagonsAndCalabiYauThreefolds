restart

kk = ZZ/1009
R1 = kk[x_0..x_3]
R2 = kk[z_0..z_3]
R = kk[x_0..x_3,z_0..z_3]

I = ideal(sub(random(3,R1),R),sub(random(3,R2),R))
I = ideal(x_0*x_1*x_2, z_0*z_1*z_2)

loadPackage "VersalDeformations"

T1 = CT^1(0, gens I)
T2 = CT^2(0, gens I)

(f1,r1,g1,c1) = versalDeformation(gens I, T1, T2);

II = ideal sub(sum f1, toList apply(1..194, i-> (t_i => random(0,R))))


Im = minimalPresentation (II + random(1,R) + random(1,R))

X = Proj(ring Im/Im)

euler X
--194-144 = 50

HH^2(cotangentSheaf X)
