restart
loadPackage "VersalDeformations"

{*
I study the degree zero-deformations of C(dP_6)xA^2. Taking their closure in
Y_0 should give a deformation of Y_0.

I find a deformation of Y_0 by looking at C(dP_6)xA^2 and taking the closure
in PP^13. It turns out that this smoothing has singular locus a union of P^1's.

(so a CI Calabi-Yau inside it is smooth.)
*}

S = QQ[x_1..x_6,y_0,z_1,z_2]
I = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)


T1 = CT^1(0, gens I); -- 6-dim
T2 = CT^2(0, gens I); --6-dim
(F,R,G,C) = versalDeformation(gens I, T1, T2); 
comps = primaryDecomposition ideal sum G
apply(comps, i -> (dim i) - dim S)
--(4,2,0). Last one is nilpotent with support at the origin

-- the two degree zero deformations
I1 = ideal sub(transpose sum F, {t_1 => 0, t_2 => 0, t_3 => 2, t_4 => 1, t_5 => 2, t_6 => 1})
I2 = ideal sub(transpose sum F, {t_1 => 1, t_2 => 2, t_3 => 3, t_4 => 7, t_5 => 1, t_6 => 2})

dim radical ideal mingens ideal singularLocus I1 -- 1
dim radical ideal mingens ideal singularLocus I2 -- 0

-- Making the ring and the ideal of dP_6**dP6
T = QQ[x_1..x_6,z_1..z_6,y_0,y_1,  Weights => {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1}]
Ik = ideal(x_1*x_3-x_2*y_0,
    	x_2*x_4-x_3*y_0,
	x_3*x_5-x_4*y_0,
	x_4*x_6-x_5*y_0,
	x_5*x_1-x_6*y_0,
	x_6*x_2-x_1*y_0,
	x_1*x_4-y_0^2,
	x_2*x_5-y_0^2,
	x_3*x_6-y_0^2)
f = map(T,T,{z_1,z_2,z_3,z_4,z_5,z_6,x_1,x_2,x_3,x_4,x_5,x_6,y_1,-y_0})
J = f Ik + Ik

-- Considering the affine chart where x_1=1.
sub1 = sub(J, x_1 => 1)
J1 = minimalPresentation(sub1) -- eliminating all x_i but x_1,x_2.

-- The transition map from the affine chart.
g = sub1.cache.minimalPresentationMap


-- Deformations in the affine chart (?)
T1 = CT^1(0, gens J1);
T2 = CT^2(0, gens J1);
(F,R,G,C) = versalDeformation(gens J1, T1, T2);

comps' = primaryDecomposition ideal sum G
apply(comps', I -> dim I - dim ring J1) -- Really too many since we have not eliminated x_1.

I1 = ideal sub(transpose sum F, {t_1 => 0, t_2 => 0, t_3 => 0, t_4 => 0, t_7 => 0, t_5 => 2, t_8 => 2, t_6 => 3, t_9 => 3})
I2 = ideal sub(transpose sum F, {t_1 => 0, t_4 => 0, t_7 => 0, t_2 => 1, t_5 => 1, t_3 => 2, t_6 => 2, t_8 => 1, t_9 => 1})

dim radical ideal mingens ideal singularLocus I1 -- 2
dim radical ideal mingens ideal singularLocus I2 -- 1


-- We pull back this deformation to PP^13.
fkand = transpose mingens homogenize(preimage(g, I2),x_1)
dim radical ideal mingens ideal singularLocus minimalPresentation (sub(ideal fkand, x_1 => 1) + x_1) -- 0
-- So fkand has an isolated singularity in the chart x_1!=0.


-- We check in the other charts
apply(1..6, i -> (dim radical ideal mingens ideal singularLocus minimalPresentation  (sub(ideal fkand, x_i => 1) +x_i)))
-- {0,1,1,1,1,-1}
apply(1..6, i -> (dim radical ideal mingens ideal singularLocus minimalPresentation  (sub(ideal fkand, z_i => 1) + z_i)))
-- {2,2,2,2,2,2}

-- So in the x_i != 0-charts we have 1-dim sings.
-- We repeat the procedure in the z_i != 0-charts.

fk = ideal fkand
JJ1 = sub(fk, z_1 => 1)
JJmin = minimalPresentation JJ1  -- this is the chart z_1 != 0 of fkand

gg = JJ1.cache.minimalPresentationMap -- the change of chart map

T1 = CT^1(0, gens JJmin); -- 9 
T2 = CT^2(0, gens JJmin); -- 12
(F,R,G,C) = versalDeformation(gens JJmin, T1, T2);

comps'' = primaryDecomposition ideal sum G
apply(comps'', I -> dim I - dim ring JJmin) -- 3, 6, 0

-- We choose the corresponding deformation (same as last time)
I2 = ideal sub(transpose sum F, {t_1 => 0, t_2 => 1, t_3 => 2, t_4 => 0, t_5 => 0, t_6 => 0, t_7 => 0, t_8 => 3, t_9 => 4})

-- We pull back to PP^13
fkand2 = transpose mingens homogenize(preimage(gg, I2),z_1)
-- (get something quite complicated)

--- We check the singular locus
apply(1..6, i -> (dim radical ideal mingens ideal singularLocus minimalPresentation  (sub(ideal fkand2, x_i => 1) +x_i)))
apply(1..6, i -> (dim radical ideal mingens ideal singularLocus minimalPresentation  (sub(ideal fkand2, z_i => 1) + z_i)))
-- looong time (about an hour)
-- (0,1,1,1,1,-1) (in the x_i charts)
-- (0,1,1,1,1,-1) (in the z_i charts)
-- So sing lokus is a union of one-dim curves.
