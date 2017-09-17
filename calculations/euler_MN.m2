restart

--- goal: compute the Euler char of X_3
kk = ZZ/1009
R = kk[x_0..x_8,y_0..y_7]
M1 = genericMatrix(R,3,3)

S = ZZ/1009[a,b,r,s,u,v]

R' = ZZ/1009[y_0..y_7]
f = map(S,R',{a*r*u,a*r*v,a*s*u,a*s*v,b*r*u, b*r*v, b*s*u, b*s*v})

I = sub(ker f,R) + minors(2,M1)

betti res I

eulers sheaf(I/I^2) --- -174

-204+174
--- conclusion euler-kar= -60

