restart

{*
Computations on some examples in Klaus Altmann's article "The
versal Deformation of an isolated toric Gorenstein Singularity".
*}

R = QQ[t_1..t_6]
I = ideal(t_1+t_2-t_4-t_5,t_2+t_3-t_5-t_6,t_1^2+t_2^2-t_4^2-t_5^2, t_2^2+t_3^2-t_5^2-t_6^2);
S = QQ[s_1..s_15];
f = map(R,S,{t_1-t_2,t_1-t_3,t_1-t_4,t_1-t_5,t_1-t_6,t_2-t_3,t_2-t_4,t_2-t_5,t_2-t_6,t_3-t_4,t_3-t_5,t_3-t_6,t_4-t_5,t_4-t_6,t_5-t_6});
J = ker f;
minimalPresentation(S/J)

h = J.cache.minimalPresentationMapInv

incl =  (f * h )
minimalPresentation ideal mingens preimage(incl,I)
radical I== I

---
A = QQ[t_1..t_6,s_1..s_3,tt]
IA= sub(I,A)
sub(IA, {tt => t_1, s_1 => t_1-t_3, s_2 => t_4-t_2, s_3 => t_1-t_4, t_1 => tt, t_2 => tt-s_2-s_3, t_3 => tt-s_1, t_4 => tt-s_3, t_5 => tt-s_2, t_6 => tt-s_1-s_3})
mingens oo
