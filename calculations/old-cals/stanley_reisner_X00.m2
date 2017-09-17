restart
R = QQ[x_1..x_6,z_1..z_6]
M1 = matrix{{0,x_1,x_2},{x_4,0,x_3},{x_5,x_6,0}}
M2 = matrix{{0,z_1,z_2},{z_4,0,z_3},{z_5,z_6,0}}
I = ideal mingens(minors(2,M1) + minors(2,M2))

loadPackage "VersalDeformations"
T1 = CT^1(0, gens I)
T1n = submatrix'(T1, {6,13,20,27,34,41,42,49,56,63,70,77})
T2 = CT^2(0, gens I);

(f1,r1,g1,c1) = versalDeformation(gens I, T1n, T2, HighestOrder => 2, SmartLift => false);
sum g1

loadPackage "MinimalPrimes"
installMinprimes()

sum g1
decompose ideal (transpose sum g1)_{0..35}

betti res I

----
--- finn polytop til join av del Pezzo
--
restart

R = QQ[y_0,x_1..x_6,y_1,z_1..z_6]
M1 = matrix{{y_0,x_1,x_2},{x_4,y_0,x_3},{x_5,x_6,y_0}}
M2 = matrix{{y_1,z_1,z_2},{z_4,y_1,z_3},{z_5,z_6,y_1}}
I = ideal mingens(minors(2,M1) + minors(2,M2))	

loadPackage "defMethods"
loadPackage "NormalToricVarieties"

C = convexHull makeCone I
C = convexHull transpose ((transpose makeCone I)_{0..4})
interiorLatticePoints (2*C)



C = posHull transpose matrix{{1,2},{1,0}}

dualCone C
rays oo
hilbertBasis dualCone C
