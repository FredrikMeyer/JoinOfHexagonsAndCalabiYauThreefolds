restart
R = QQ[x_0..x_4]
IX = ideal (sum apply(gens R,m -> m^5) + product gens R)

loadPackage "Binomials"
torus = ideal apply(flatten apply(apply(apply(flatten entries gens IX, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)
toruskomps = BPD torus
toruskomps = select(toruskomps, I -> dim I == 1)



