restart
kk = ZZ/101
R = QQ[x_0..x_15]--, Weights => {1,2,2,1,2,2,2,2}]
S = QQ[r,s,u,v,a,b,r',s',u',v',a',b']

f = map(S,R,{r*u*a, r*u*b,r*v*a,r*v*b, s*u*a,s*u*b,s*v*a,s*v*b,
       r'*u'*a', r'*u'*b',r'*v'*a',r'*v'*b', s'*u'*a',s'*u'*b',s'*v'*a',s'*v'*b'})
I = ker f
A = R/I
M = Proj(A)
B = A/sub(ideal( toList apply(1..6, i->random(1,R))),A)

CTM = time prune cotangentSheaf M;
     -- used 746.391 seconds (12 min)
(I,CTM) = value get "cotangentSheafP1join.txt"     ;
CTMX = sheaf ((CTM) ** B);

IMJ = sheaf image (A ** jacobian I)
IX = I + ideal toList apply(1..6, i->random(1,R));

--"cotangentsheafP1join.txt" << toString (I,module CTM) << endl << close

genera CTM
-- o16 = {4, 12, 30, 30, 25, 155, 397, 251}
eulers CTM
--i14 : eulers CTM
--
--o14 = {-3, 13, -29, 31, -24, 156, -396, 252}


eulers sheaf (I/I^2)
16*12
-192+168
-24+4*12
-- this proves that X_2 has Euler characeristic -48??
time eulers sheaf (I/I^2)
euler sheaf image (jacobian I ** A)

II = (sheaf (I/I^2));
imd = sheaf image (jacobian I ** A);
euler ((sheaf (I/I^2))(-2))
sum toList apply(0..6, i-> ((-1)^i) * binomial(6,i) * euler (II(-i)))
sum toList apply(0..6, i-> ((-1)^i) * binomial(6,i) * euler (imd(-i)))





---- cotangentSheaf M
---- compute cotangent sheaf of M.
---- Confirm Euler-char of X_1.
restart
kk = ZZ/101
R = QQ[x_0..x_17]
M1 = genericMatrix(R,3,3)
M2 = genericMatrix(R,x_9,3,3)
I = minors(2,M1) + minors(2,M2)
M = Proj(R/I)
A = R/I

eulers sheaf (I/I^2)
imd = sheaf image (jacobian I ** A);
eulers imd

CTM = time cotangentSheaf M;
