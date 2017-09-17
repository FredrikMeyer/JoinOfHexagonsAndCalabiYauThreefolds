restart
kk = ZZ/3001
Z = QQ[x_1..x_12]

pars = {2,3,5}
pars = {1,1,1}
fija = (i,j,a) -> (
    Eij := (id_(Z^3))_{i} * transpose (id_(Z^3))_{j};
    Eij' :=  (id_(Z^3))_{(-i-j) % 3} * transpose (id_(Z^3))_{(-i-j) % 3};
    if (a == 0) then (
    Eij | pars#((i+j)%3) * Eij'
    )
    else (
    pars#((i+j)%3) * Eij' | Eij
    )
    )

MG = x_1*fija(0,1,0) + x_2*fija(0,2,0) + x_3*fija(1,0,0) + x_4*fija(1,2,0) + x_5*fija(2,0,0) + 
    x_6*fija(2,1,0) + x_7*fija(0,1,1) + x_8*fija(0,2,1) + x_9*fija(1,0,1) + x_10*fija(1,2,1) + x_11*fija(2,0,1) + x_12*fija(2,1,1)

IX = minors(2,MG_{0,1,2}) + minors(2,MG_{3,4,5})
time gens gb IX;

loadPackage "Binomials"
torus = ideal apply(flatten apply(apply(apply(flatten entries gens IX, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)
toruskomps = BPD torus
toruskomps = select(toruskomps, I -> dim I == 1) --- har en Z_3-virkning!!! (altså den med f_ij)




loadPackage "MinimalPrimes"
installMinprimes()

IX + ideal(x_2)
decompose oo

--- try to find local equations for (1:0:0:0...)
loadPackage "LocalRings"

I1 = minimalPresentation(IX + ideal(x_1-1))
S1 = I1_*

f1 =S1#7
use ring f1
u1 = (f1 - (f1 % x_9))/x_9

S2 = apply(S1, h -> numerator sub(h, x_9 => -(f1 % x_9)/u1))


f2 = S1#8
u2 = (f2 - (f2 % x_2))/x_2

S3 = apply(S2, h -> numerator sub(h, x_2 => -(f2 % x_2)/u2))

f3 = S3#4
u3 = (f3 - (f3 % x_6))/x_6

S4 = apply(S3, h -> numerator sub(h, x_6 => -(f3  % x_6)/u3))


setMaxIdeal ideal gens ring u3

matrix{select(S4, f -> f != 0)}

setMaxIdeal ideal gens ring oo
localMingens o157
