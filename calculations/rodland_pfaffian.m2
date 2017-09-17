restart

R = ZZ/1009[x_0..x_6]
S = ZZ/1009[y_1,y_2,y_3]

T = R ** S

M1 = matrix{{0,x_1*y_1, x_2*y_2, x_3*y_3, -x_4*y_3, -x_5*y_2, -x_6*y_1},
           {0,      0 ,x_3*y_1, x_4*y_2, x_5*y_3, -x_6*y_3,  -x_0*y_2},
	   {0,      0 ,      0, x_5*y_1, x_6*y_2,  x_0*y_3,   -x_1*y_3},
	   {0,      0 ,      0,       0, x_0*y_1,  x_1*y_2,   x_2*y_3},
	   {0,      0 ,      0,       0,     0 ,   x_2*y_1,   x_3*y_2},
	   {0,      0 ,      0,       0,     0 ,         0,   x_4*y_1},
	   {0,      0 ,      0,       0,     0 ,         0,       0}}

M  = M1 - transpose M1
use T
m_0 = lift(sub(M, {x_0 => 1}) ** T/ideal(x_0..x_6),T)
apply(0..6, i-> lift(sub(M, {x_i => 1}) ** T/ideal(x_0..x_6),T))

Iy = pfaffians(6,M)

suby = {y_1 => 0, y_2 => sub(random(0,S),T), y_3 => sub(random(0,S),T)}
suby = {y_1 => 0, y_2 => 0, y_3 => 1_T}
suby = {y_1 => sub(random(0,S),T), y_2 => sub(random(0,S),T), y_3 => sub(random(0,S),T)}
suby = {y_1 => 1_T, y_2 => 1_T, y_3 => 1_T} ---- Y = (1:1:1), gir X med 49 sings
sub(M, suby)
generisk = sub(sub(M, suby),R)

--use R
p1 = sub(generisk, {x_0 => 0, x_1 => 1, x_2 => 1, x_3 => 1, x_4 => -1, x_5 => -1, x_6 => -1})

s = transpose matrix{{0,1,0,0,0,0,0},{0,0,1,0,0,0,0},{0,0,0,1,0,0,0},
    {0,0,0,0,1,0,0},{0,0,0,0,0,1,0},{0,0,0,0,0,0,1},{1,0,0,0,0,0,0}}
w = (select(1..1009, i-> i_T^7 == 1))#1
t = sub(transpose matrix{{1,0,0,0,0,0,0},{0,w,0,0,0,0,0},{0,0,w^2,0,0,0,0},{0,0,0,w^3,0,0,0},
        {0,0,0,0,w^4,0,0},{0,0,0,0,0,w^5,0},{0,0,0,0,0,0,w^6}},T)

G = flatten toList apply(0..6, i -> toList flatten apply(0..6, j-> (inverse (s^i*t^j)*p1*s^i*t^j)))

p1
IX = sub(pfaffians(6,generisk),R)

mingens ideal singularLocus IX
degree ideal oo
decompose ideal oo

--select(0..1009, i-> i_T^7 == 1)
NX = sheaf prune Hom(IX/IX^2,R^1/IX)
rank HH^0(NX) - 49-1

X = Proj(R/IX)


singi = ideal mingens ideal singularLocus IX;
rs = radical singi

(decompose rs)#30
use ring generisk
sub(generisk, {x_4 => 0, x_6 => 1, x_5 => 1, x_3 => -1, x_2 => -1, x_1 => -1, x_0 => 1}) * -1






loadPackage "MinimalPrimes"
installMinprimes()
#decompose singi --- <-49, så like mange som Rødland hevder




torus = ideal apply(flatten apply(apply(apply(flatten entries gens IX, monomials), v ->  flatten entries v), j -> subsets(j,2)),    s -> s_0-s_1)

loadPackage "Binomials"
tt = (select(BPD torus, I -> dim I == 1))#1
minimalPresentation tt

-- generisk pfaffian
restart
kk = QQ
R = kk[x_0..x_20]
M = genericSkewMatrix(R,7)

ms  = apply(0..6, i-> sub(M,random(R^1,R^21)))

S = kk[a_0..a_6]

IX = pfaffians(6, sum toList apply(0..6, i -> sub(ms#i,S)*a_i))


loadPackage "VersalDeformations"
--CT^1(0, gens IX)

IIX = sheaf prune (IX/IX^2);
HH^3(IIX)
98-49+1 == 50 -- :D
98-48 == 50 -- :D
X = Proj(S/IX)

prune Hom(prune (IX/IX^2),S^1/IX)

ring IX
KM = (sheaf prune (Ext^3(S^1/IX, S^1) ** (S/IX)))
apply(0..3, i-> HH^i(KM))
betti res IX
X = Proj(S/IX)
(sheaf (Ext^3(S^1/IX, S^1) ** (S/IX)) ** OO_X(-7))
isFreeModule oo

HH^0(prune sheaf oo)



HH^0(OO_X(1))

cotangentSheaf X;