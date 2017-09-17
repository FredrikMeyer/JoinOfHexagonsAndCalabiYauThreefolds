restart
loadPackage "VersalDeformations"
kk = ZZ/32003
kk = QQ
w =  {1, 1, 4, 7, 7, 4, 1, 1, 4, 7, 7, 4, 1, 1}
w2 =  {7, 7, 8, 9, 8, 7, 1, 5, 14, 19, 14, 5, 5, 7}
w3 = {4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 1, 1}
T = kk[x_1..x_6,z_1..z_6,y_0,y_1,  Weights => w] --- w2
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
J = Ik + f Ik
R = kk[x_1..x_6, z_1..z_6, Weights=> w_{0..11}]
I = sub(ideal leadTerm J,R)

T1 = CT^1(0, gens I);
T2 = CT^2(0, gens I);

hilbertPolynomial I
hilbertPolynomial(I, Projective  => false)
     

(f1,r1,g1,c1) = versalDeformation(gens I, T1, T2, HighestOrder => 20, SmartLift => false);
T1 = T1_{6,13,14,21,34,41,42,49,62,69,70,77}
T1 = submatrix'(T1,{6,13,14,21,34,41,42,49,62,69,70,77}) 

ones = transpose matrix{{1,0,1,0,0,0}}
T11 = T1_{0..5} * ones
for i from 1 to 11 do (
    T11 = T11 | (T1_{(6*i)..(6*i+5)}*ones);
)
T11
ones = transpose matrix{{1,0,0,0,0,0}}
ones' = transpose matrix{{0,0,1,0,0,0}}


T11 = T1_{0..5} * ones | T1_{0..5} * ones'
for i from 1 to 11 do (
    T11 = T11 | (T1_{(6*i)..(6*i+5)}*ones);
        T11 = T11 | (T1_{(6*i)..(6*i+5)}*ones');
)
T11

(f1,r1,g1,c1) = versalDeformation(gens I, T11, T2, HighestOrder => 20);--, SanityCheck => false);
transpose sum f1
dim ideal sum g1 - dim R

transpose sum f1

IG = ideal mingens ideal g1
loadPackage "Binomials" -- prime ^
loadPackage "defMethods"
transpose mingens IG
IF = sub(transpose sum f1, toList apply(1..24, i -> t_i => -1))
IF = sub(transpose sum f1, {t_24 => 1, t_6 => 1, t_8 => 1, t_22 => 1, t_14 => 1, t_16 => 1, t_7 => 2, t_21 => 2, t_13 => 2, t_15 => 2, t_23 => 2, t_5 => 2, t_4 => 3, t_18 => 3, t_10 => 3, t_12 => 3,t_20 => 3, t_2 => 3, t_3 => 4, t_9 => 4, t_19 => 4, t_17 => 4, t_11 => 4, t_1 => 4})
sub(IG,  toList apply(1..24, i -> t_i => -1))

loadPackage "gfanInterface"

(ideal IF)_*
I_*
w = weightVector(I_*, (ideal IF)_*)
minp = ideal transpose mingens minimalPresentation (ideal IF + (x_4-1))
mingens ideal singularLocus minp
dim ideal oo

w
S = kk[x_1..x_6,z_1..z_6, Weights => w]
IX00  =  ideal sub(IF, S)

transpose mingens IX00
gens gb (IX00^2)
A = S/IX00
HH^0(sheaf prune Hom( prune ((module IX00) ** A), A))

transpose gens IX00

X00 = Proj(S/IX00)
HH^0(OO_X00(1))

IX00' = ideal transpose mingens sub(IX00, {x_1 => (x_1-3*x_3)/4, z_6 => z_6 + 24*x_2})


ideal transpose mingens ideal mingens  minimalPresentation(IX00 + (x_4-1))
mingens ideal singularLocus oo
radical ideal oo



apply(oo, dim)
for i from 1 to 6 do {
 minp =     ideal mingens  minimalPresentation(IX00 + (x_i-1));
 print dim ideal mingens ideal singularLocus minp;
    }

singlist = {}
for i from 1 to 6 do {
    sz = sub(IX00, z_i => 1);
    sx = sub(IX00, x_i => 1);
    singz = radical ideal mingens ideal singularLocus  minimalPresentation sz;
    singx = radical ideal mingens ideal singularLocus  minimalPresentation sx;
    singsx = (decompose singx);
    singsz = (decompose singz);
    sings = singsx | singsz;
    invz = sz.cache.minimalPresentationMap;
    invx = sx.cache.minimalPresentationMap;
    singlist = singlist | apply(singsx, I -> homogenize(preimage(invx, singx),x_i)) | apply(singsz, I -> homogenize(preimage(invz, singz),z_i));
    print i;
    }
2
