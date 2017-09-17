restart
{*
I study deformations of C(dP_6).
*}
S = QQ[x_1..x_6,y_0]
I = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)

loadPackage "VersalDeformations"

(F,R,G,C) = versalDeformation(gens I);

-- Making a common ring
T = QQ[x_1..x_6,y_0,t_1,t_2,t_3,s_1,s_2,s_3];
gsub = sub(sub(ideal mingens ideal G,T), {t_2 => s_2+s_3-s_1, t_1 => s_3, t_3 => s_3-s_1})
fsub = transpose  sub(sub(sum F,T), {t_2 => s_2+s_3-s_1, t_1 => s_3, t_3 => s_3-s_1})

decompose gsub
--- The two different smoothing of C(dP6).
fsub1 = sub(fsub, s_1 => 0)
fsub2 = sub(fsub, {s_3 => 0, s_2 => 0})


I2 = sub(sub(fsub2, s_1 => 1), S)
I1 = sub(sub(fsub1, {s_2 => 1, s_3 => 1}),S)

-- both are smooth:
radical ideal mingens ideal singularLocus ideal I2
radical ideal mingens ideal singularLocus ideal I1

-- I compute if there is a 
-- torus action on I2:
torus = ideal apply(flatten apply(apply(apply(flatten entries I2, monomials), v ->  flatten entries v), j -> subsets(j,2)), s -> s_0-s_1)
-- Radical = true 
radical torus == torus

loadPackage "defMethods"
makeCone first decompose torus
-- the answer is the matrix of points in 
-- the polytope of dP6!