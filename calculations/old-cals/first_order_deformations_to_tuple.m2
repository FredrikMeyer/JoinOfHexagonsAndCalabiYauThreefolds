restart
loadPackage "VersalDeformations"
loadPackage "SimplicialComplexes"
-- Kuhnels 9 vertex
T = {(2, 5, 6, 7, 8), (1, 3, 5, 7, 9), (2, 4, 5, 8, 9), (1, 2, 5, 6, 8), (1, 5, 6, 8, 9),
    (1, 2, 3, 7, 8), (4, 6, 7, 8, 9), (2, 3, 6, 7, 9), (1, 2, 4, 5, 6), (1, 2, 4, 6, 7),
    (1, 3, 4, 5, 7), (2, 4, 6, 7, 9), (1, 2, 4, 5, 9), (3, 4, 5, 7, 8), (3, 5, 6, 7, 9),
    (1, 2, 3, 8, 9), (2, 3, 4, 8, 9), (2, 3, 4, 6, 9), (2, 3, 5, 7, 8), (2, 3, 5, 6, 7),
    (2, 3, 4, 5, 6), (2, 3, 4, 5, 8), (1, 3, 4, 5, 6), (3, 4, 6, 8, 9), (1, 2, 6, 7, 8),
    (1, 3, 6, 8, 9), (1, 2, 3, 7, 9), (1, 2, 5, 8, 9), (5, 6, 7, 8, 9), (1, 3, 5, 6, 9),
    (1, 4, 6, 7, 8), (4, 5, 7, 8, 9), (1, 2, 4, 7, 9), (1, 4, 5, 7, 9), (1, 3, 4, 7, 8), (1, 3, 4, 6, 8)}

S = ZZ/1009[x_1..x_9]
mons = apply(T, t -> product toList apply(t, i -> x_i))
K =  simplicialComplex mons
I = ideal K

F0 = gens I;
R0 = gens ker F0;
T1 = CT^1(0,F0);
(F,R) = firstOrderDeformations(F0,R0,T1);

f0 = sub(((F#0)_0)_0,S)
f1 = sub(sub((F#1)_0_0, t_1 => 1),S)

L = {f0,f1}
v0 = toSequence (exponents f0)#0
v1 = toSequence (exponents f1)#0
(v0,v1)

transpose sum F

-- enter without the t's
toTuple = L -> (
    f0 := L#0;
    f1 := L#1;
    v0 := toSequence (exponents f0)#0;
    v1 := toSequence (exponents f1)#0;
    (v0,v1)
)

toTuple L

listOfTuples = {}
for i from 0 to 20 do (
    i = 1
    f0 = sub(((F#0)_i)_0,S)
    f1 = sub(sub((F#1)_i_0, t_(i+1) => 1),S);
    L = {f0,f1};
    listOfTuples = listOfTuples | {toTuple L};
    )

(F#1)_0_0

