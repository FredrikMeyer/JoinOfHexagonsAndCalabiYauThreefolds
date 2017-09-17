restart
kk = QQ
R = kk[x_1..x_6,z_1..z_6,y_0,y_1,  Weights => {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1}]

I = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)
f = map(R,R,{z_1,z_2,z_3,z_4,z_5,z_6,x_1,x_2,x_3,x_4,x_5,x_6,y_1,-y_0})	
J = f I + I

loadPackage "defMethods"
loadPackage "NormalToricVarieties"

makeCone J
C = (makeCone J)^{0..4}
P' = convexHull C
P = convexHull transpose matrix apply(entries transpose vertices (2*P'), v-> v - first entries transpose first interiorLatticePoints (2*P'))
isReflexive P
V = normalToricVariety normalFan P

W = blowup(toList (0..5), V, {0,0,0,0,1})
W = blowup(toList(6..11), W, {0,0,0,0,-1})

isProjective W
KXV = sum toList apply(0..11, i-> V_i)
KXW = sum toList apply(0..13, i-> W_i)
projEmb KX

M = fromWDivToCl W
M * transpose matrix{toList apply(1..14, i-> 1)}
N = fromWDivToCl V
N * transpose matrix{toList apply(1..12, i-> 1)}

KXV

rays W
isCartier KXV

rays first cones(5, fan V)
contains(first cones(5, fan V), transpose matrix {{0,0,0,0,-1}})
transpose ooo
b  = -transpose matrix {{1,1,1,1,1,1,1,1}} 
H \\ b
transpose H

--loadPackage "ToricCompleteIntersections"
--h11OfCY polar polytope KX
--h21OfCY polar polytope KX

latticePoints polytope KXV
