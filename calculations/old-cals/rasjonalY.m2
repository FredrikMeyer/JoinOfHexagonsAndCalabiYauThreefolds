restart
kk = QQ
R = kk[x_1..x_6,y_0,z_1..z_6,y_1,  Weights => {1, 1, 4, 7, 7, 4,1, 1, 1, 4, 7, 7, 4,  1}]

I = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)
f = map(R,R,{z_1,z_2,z_3,z_4,z_5,z_6,y_1,x_1,x_2,x_3,x_4,x_5,x_6,-y_0})
J = f I + I

decompose ideal matrix f

S = QQ[r,s,t,u,v,w]
radical ideal (t^4,r*t^3,r^2*s*t,r^2*s^2,r*s^2*t,s*t^3,r*s*t^2)
PP2 = QQ[r,s,t]
PP2' = QQ[u,v,w]

use S
g = map(S,R, sub((basis(3,PP2))_{1,2,3,4,5,7,8},S) | sub((basis(3,PP2'))_{1,2,3,4,5,7,8},S))
g = map(S,R, {r*t^2,r^2*t,r^2*s,r*s^2,s^2*t,s*t^2,r*s*t,
	u*w^2,u^2*w,u^2*v,u*v^2,v^2*w,v*w^2,u*v*w})

loadPackage "SimplicialComplexes"
simplicialComplex radical ideal matrix g
radical ideal matrix g

matrix g
g' = matrix g + matrix{{0,0,0,0,0,0,0,0,0,0,0,0,0,0}}
transpose mingens ker map(S,R,g')
hilbertPolynomial ker map(S,R,g')
hilbertPolynomial ker map(S,R,g)
dim ker map(S,R,g)
ideal singularLocus ker map(S,R,g')




BlY = QQ[a_0..a_5,r,s,t,u,v,w]
BY = minors(2,matrix{{a_0,a_1,a_2,a_3,a_4,a_5},{r*s,r*t,s*t,u*v,u*w,v*w}})
B = first decompose BY

rB = reesIdeal radical ideal matrix g

ideal mingens(rB + ideal jacobian rB)


----
-- lek med simplicialkomplesker
restart
R = QQ[x_1,x_2,x_3,r_1,r_2,r_3,s_1,s_2,s_3,t_1,t_2,t_3]
loadPackage "SimplicialComplexes"
S1 = facets boundary simplicialComplex {x_1*r_1*r_2*r_3}
S2 = facets boundary simplicialComplex {x_2*s_1*s_2*s_3}
S3 = facets boundary simplicialComplex {x_3*t_1*t_2*t_3}
S = simplicialComplex (flatten entries (S1 | S2 | S3) | {x_1*x_2, x_2*x_3, x_3*x_1})
prune homology S


