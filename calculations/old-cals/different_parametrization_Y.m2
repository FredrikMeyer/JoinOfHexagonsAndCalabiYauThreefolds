restart
R = QQ[a,b,c,u,v,w,s,t]
S = QQ[x_1..x_6,z_1..z_6,y_0,y_1]

f = map(R,S,{c^4, a*c^3, a^2*b*c, a^2*b^2, a*b^2*c, b*c^3,
	w^4, u*w^3, u^2*v*w, u^2*v^2, u*v^2*w, v*w^3,a*b*c^2,u*v*w^2})

g = map(R,S,{s*c^4, s*a*c^3, s*a^2*b*c, s*a^2*b^2, s*a*b^2*c, s*b*c^3,
	t*w^4, t*u*w^3, t*u^2*v*w, t*u^2*v^2, t*u*v^2*w, t*v*w^3,s*a*b*c^2,t*u*v*w^2})

transpose mingens ker g

use R
intersect decompose radical saturate(saturate(saturate(ideal matrix g, ideal(s,t)),ideal(u,v,w)),ideal(a,b,c))

T = QQ[a,b,c,u,v,w,s,t,q_1..q_8]
II = sub(oo, T)
M = matrix{{a*b*u*v, c*u*v, a*b*w, c*w, a*b*s, c*s, u*v*t, w*t},{q_1,q_2,q_3,q_4,q_5,q_6,q_7,q_8}}
BI =first decompose minors(2,M)

BI




transpose mingens eliminate(eliminate(ker g, y_0), y_1)
sub(matrix g, c => 1)
