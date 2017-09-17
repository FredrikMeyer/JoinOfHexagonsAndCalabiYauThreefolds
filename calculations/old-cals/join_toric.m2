restart

R = QQ[y_0,x_1..x_6,y_1,z_1..z_6]
M1 = matrix{{y_0,x_1,x_2},{x_4,y_0,x_3},{x_5,x_6,y_0}}
M2 = matrix{{y_1,z_1,z_2},{z_4,y_1,z_3},{z_5,z_6,y_1}}

I = minors(2,M1) + minors(2,M2)

loadPackage "defMethods"

V =  transpose ((transpose makeCone I)_{0..4})

loadPackage "NormalToricVarieties"

P = convexHull V

F = normalFan P
X = normalToricVariety F

D = sum toList apply(0..11, i -> X_i)
vertices polytope D
latticePoints (1/2 * polytope D)
latticePoints D
isReflexive polytope D 
isFano X

