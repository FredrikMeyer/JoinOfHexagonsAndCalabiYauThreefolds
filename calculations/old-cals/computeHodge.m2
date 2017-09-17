restart
loadPackage "Polyhedra"
loadPackage "NormalToricVarieties"

V = projectiveSpace(4)
KX= fold(plus,toList apply(0..4, i-> V_i))

P =  polytope KX
toString vertices P
vs =  matrix {{-1, 4, -1, -1, -1}, {-1, -1, 4, -1, -1}, {-1, -1, -1, 4, -1}, {-1, -1, -1, -1, 4}}
P = convexHull vs
#interiorLatticePoints(2*P)
#interiorLatticePoints(P)
126-5*1-20 -- :D riktig

faceOf(P)
F = first faces(2,P)
polarFace F
code(polarFace,Polyhedron)

polarFace2 = method();
polarFace2(Polyhedron,Polyhedron) := (F,P) -> (
	       V := transpose vertices F;
	       R := transpose rays F;
	       P0 := P;
	       P0d := polar P0;
	       codimensionPd := dim F - P0#"dimension of lineality space" + 1;
	       L := faces(codimensionPd,P0d);
	       Pd = first select(1,L, l -> all(flatten entries(V*(vertices l)),e -> e == -1) and V*(rays l) == 0 and R*(vertices l | rays l) == 0);
	       Pd)

(#latticePoints P) - 5 - sum(faces(1,P), F -> #interiorLatticePoints F) + sum(faces(2,P), F -> (#interiorLatticePoints F)*(#interiorLatticePoints polarFace2(F,P)))


apply(faces(3,P), F ->  vertices polarFace F)
sum apply(faces(1,P), F -> #interiorLatticePoints F)


M = matrix {{1, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -1, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0,
      1, 0, -1, -1, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 1, 1, 0, 0}, {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0}}
P =  convexHull transpose ((transpose M)_{0..3})
(#latticePoints P) - 5 - sum(faces(1,P), F -> #interiorLatticePoints F) + sum(faces(2,P), F -> (#interiorLatticePoints F)*(#interiorLatticePoints polarFace2(F,P)))

V = normalToricVariety P
VV = normalToricVariety polar P

vertices polar P
makeSmooth V
Vsm = oo
fan Vsm 

toricIdeal rays first cones(4,fan V)
lift( transpose((transpose vertices P) |  transpose matrix{{1,1,1,1,1,1,1,1,1,1,1,1}}),ZZ) | transpose matrix{{0,0,0,0,1}}
TI = toricIdeal oo

transpose mingens  minimalPresentation(TI + ideal(m-1))
transpose mingens  minimalPresentation(TI + ideal(a-1))
TI1 = ideal oo
(transpose (transpose makeCone TI1)_{0..2})
P1 = convexHull ((transpose (transpose makeCone TI1)_
	{0..2}))

vertices polar P1

V1  = normalToricVariety P1
fan makeSmooth V1
fan V1
8*12

hilbertPolynomial TI
TI
IX0=  (TI + random(1,ring TI))

X00 = Proj(ring TI/(TI+random(1,ring TI)))

euler sheaf prune (IX0/IX0^2) -- -180
180-156
T1X = CT^1(0, gens ideal X00)
T2X = CT^2(0, gens ideal X00);
(f1,r1,g1,c1) = versalDeformation(gens ideal X00, T1X_{0..10},T2X, HighestOrder => 2);




betti res TI
betti res TI1

loadPackage "VersalDeformations"
apply(cones(1,fan V), C -> isSmooth C)

CayleyPoly = convexHull transpose  (matrix ({{0, 0, 1, 0, -1, 1, 0}, {0, 0, 0, 1, -1, 1, 0}, {0, 0, 1, -1, -1, 1, 0}, {0, 0, 0, -1, -1, 1, 0}, {0, 0, -1, 1, -1, 1, 0}, {0, 0, -1, 0, -1, 1, 0}, {0, 0, 0, 0, 0, 1, 0}, {0, -1, 0, 0, 1, 0, 1}, {-1, 0, 0, 0, 1, 0, 1}, {1, -1, 0, 0, 1, 0, 1}, {-1, 1, 0, 0, 1, 0, 1}, {0, 1, 0, 0, 1, 0, 1}, {1, 0, 0, 0, 1, 0, 1}, {0, 0, 0, 0, 0, 0, 1}}))_{0..5}

CayleyPoly =  convexHull transpose matrix apply(entries transpose vertices CayleyPoly, v -> v + {0,0,0,0,0,-1/2})

isReflexive (2*CayleyPoly)
interiorLatticePoints polar (2*CayleyPoly)

latticePoints polar P

#latticePoints CayleyPoly
1 - 7 - 0 + 0 + ____ - 0 - 0 + 0

P = polar CayleyPoly

sum apply(faces(6,P), F -> #interiorLatticePoints(2*polarFace2(F,P)))

sum apply(faces(1,CayleyPoly), F -> #interiorLatticePoints F)

sum apply(faces(5,P), F -> #interiorLatticePoints(F)* #(interiorLatticePoints(2 * polarFace2(F,P))))

sum flatten apply(faces(3,P), y -> apply(faces(1,y), x -> #interiorLatticePoints(x) * #interiorLatticePoints(polarFace2(y,P))))
sum apply(faces(4,P), F -> (#interiorLatticePoints(2*F))*#(interiorLatticePoints(polarFace2(F,P))))
sum apply(faces(3,P), F -> (#interiorLatticePoints(2*F))*#(interiorLatticePoints(polarFace2(F,P))))
