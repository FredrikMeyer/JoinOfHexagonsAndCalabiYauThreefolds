restart

--w =  {6, 7, 11, 6, 7, 11, 9, 10, 1, 12, 7, 2}
T = ZZ/1009[x_0..x_11]--, Weights => w] -- ring for P^11
rms1 = toList apply(0..11, i -> random(T^3,T^3))
rms2 = toList apply(0..11, i -> random(T^3,T^3))
--- lag i stedet tilfeldig rang1-matriser
rms1 = toList apply(0..11, i -> (random(T^3,T^1) * random(T^1,T^3)))
rms2 = toList apply(0..11, i -> (random(T^3,T^1) * random(T^1,T^3)))
apply(oo, rank)
--
Z = 0* id_(ZZ^3)
blocks = apply(0..11, i-> (rms1#i | Z) || (Z | rms2#i))
--rank matrix toList apply(blocks, M -> flatten entries M) -- 12

A = sum toList apply(0..11, i-> x_i *blocksnew#i)
A = sum toList apply(0..5, i-> x_i *blocks#i)
IX =  (minors(2,A_{0..2}^{0..2}) + minors(2,A_{3..5}^{3..5}))






radical IX

mingens (IX + ideal(A_{3..5}^{3..5}))
dim ideal oo


dim IX
gens gb IX

minimalPresentation ideal mingens (ideal (A_{3..5}^{3..5}) + ideal det(A_{0..2}^{0..2}))

radical (ideal (A_{3..5}^{3..5}) + minors(2,A_{0..2}^{0..2}))


IX0 = ideal mingens minimalPresentation(IX + ideal(x_1-1))


sub(A, {x_0 => 1, x_1 => 0, x_2 => 0, x_3 => 0, x_4 => 0, x_5 => 0, x_6 => 0, x_7 => 0, x_8 => 0, x_9 => 0,x_10=>0, x_11 => 0})

ideal mingens ideal singularLocus IX0

gens gb IX




TT = ZZ/1009[x_0..x_7]
bb = apply(blocksnew, b -> sub(b,TT))
AA = sum toList apply(0..7, i -> x_i * bb#i)
IX2 = minors(2,AA)

X2 = Proj(TT/IX2)
mingens ideal singularLocus X2

cotangentSheaf X2
CT = prune oo
TS = tangentSheaf X2

euler X2
apply(0..3, i-> HH^i(CT))
apply(0..3, i-> HH^i(OO_X2))
apply(0..3, i -> HH^i(TS))



N = matrix toList apply(blocks, M -> flatten entries M)
N = matrix toList apply(rms1_{0..7}, M -> flatten entries M)

rotateMatrix = (M) -> (
    r := rank source M;
    c := rank target M;
    return matrix table(r,c, (i,j) -> M_(c-j-1,r-i-1));
    )

a = gens gb rotateMatrix N
l = rotateMatrix leadTerm a
a = rotateMatrix a

(entries a)#0
blocksnew = apply(entries a, ff -> matrix toList apply(0..5, i-> ff_{6*i..6*i+5}))


restart
R = QQ[x_0,x_1,x_2]
M = random(R^3,R^3)
M0 = random(R^3,R^1) * random(R^1,R^3)
M1 = random(R^3,R^1) * random(R^1,R^3)
M2 = random(R^3,R^1) * random(R^1,R^3)

minors(2,M + x_0*M0 + x_1*M1 + x_2*M2)

---

M1 = matrix toList apply(blocks, m -> flatten entries m)
-- V_0-basis
basisE = set apply({matrix{{1,0,0}},matrix{{0,1,0}},matrix{{0,0,1}}},transpose)
basisEE = apply(toList(basisE ** basisE), S -> S#0 * transpose S#1)
V0b = apply(basisEE, M -> ( M | Z) || (Z | Z))
-- V1-basis:
V1b = apply(basisEE, M -> ( Z | Z) || (Z | M))
M2 = matrix toList apply(V0b, m -> flatten entries m)
M2' = matrix toList apply(V1b, m -> flatten entries m)

U1 = gens ker (transpose M1 | transpose M2)
U1' = gens ker (transpose M1 | transpose M2')
first entries transpose U1
rank image  (transpose M1 | transpose M2)
ll0 = (entries (U1_0))_{0..11}
ll1 = (entries (U1_1))_{0..11}
ll2 = (entries (U1_2))_{0..11}
ll0' = (entries (U1'_0))_{0..11}
ll1' = (entries (U1'_1))_{0..11}
ll2' = (entries (U1'_2))_{0..11}
b0 = sum toList apply(0..11, i-> ll0#i * blocks#i)
b1 = sum toList apply(0..11, i-> ll1#i * blocks#i)
b2 = sum toList apply(0..11, i-> ll2#i * blocks#i)
b0' = sum toList apply(0..11, i-> ll0'#i * blocks#i)
b1' = sum toList apply(0..11, i-> ll1'#i * blocks#i)
b2' = sum toList apply(0..11, i-> ll2'#i * blocks#i)




R = QQ[x_0..x_2]
M = (x_0 * sub(b0,R) + x_1*sub(b1,R) + x_2*sub(b2,R))_{0..2}^{0..2}





blocksnew = apply(entries U1, ff -> matrix toList apply(0..5, i-> ff_{6*i..6*i+5}))



---- TODO
-- skriv H som sum H_0 + H_1 + HH hvor H_0=H \cap V_0 osv (V_0 er matrisene med rank 0 + rang *)
-- da ligger X helt inneholdt i HH som er seks-dim. 
-- vi burde da få X beskrevet som noe i P^5 (og det er *to* kandidater...!)