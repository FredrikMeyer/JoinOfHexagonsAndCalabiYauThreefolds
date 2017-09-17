restart
M = matrix{
 {-2, -2,  0,  0, -1},
 { -2,  0,  0,  0, -1},
 { 0, -2,  0,  0, -1},
 { 0,  0, -2, -2,  1},
 { 0,  0, -2,  0,  1},
 { 0,  0,  0, -2,  1},
 { 0,  0,  0,  2,  1},
 { 0,  0,  2,  0,  1},
 { 0,  0,  2,  2,  1},
 { 0,  2,  0,  0, -1},
 { 2,  0,  0,  0, -1},
 { 2,  2,  0,  0, -1}}

loadPackage "NormalToricVarieties"

P = convexHull transpose M
F = normalFan P
V = normalToricVariety  F

projEmb = (D) -> (
    X = variety D;
    L = OO D;
    m = rank HH^0(X,L);
    S = ring X;
    R = QQ[y_0..y_(m-1)];
    phi = map(S,R, basis(-first degrees L, S));
    kernel phi)

D = fold(plus,apply(0..(#rays V-1), i -> V_i))

cDiv V
wDiv V
r1 = fromCDivToPic V
r2 = fromCDivToWDiv V
r3 = fromPicToCl V
r4 = fromWDivToCl V

projEmb D

VV = makeSimplicial V
fan V
--apply(0..11, i->  rays ((cones(5,fan V))_0))
C = rays (cones(5,fan V))_0

-------------
C = posHull transpose matrix {{1, 0, 0, 0, -1}, {0, 1, 0, 0, -1}, {0, 0, -1, 1, 1}, {0, 0, 1, 0, 1}, {0, 0, 0, -1, 1}}

l = hilbertBasis dualCone C
transpose matrix apply(l, i -> flatten entries i)

loadPackage "defMethods"
rays C
I = toricIdeal oo

I


R = ring I
mins = apply(gens R, r ->  minimalPresentation sub(I, r => -1))
for i from 0 to (#mins-1) do (
    print dim ideal mingens ideal singularLocus mins#i
    )
