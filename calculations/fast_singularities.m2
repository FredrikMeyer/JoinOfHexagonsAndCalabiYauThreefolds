restart
{*
Script for fast comoputation of singular locus. Based on assumption
that equations are easier in local charts.
    *}
R = QQ[x_0..x_8,x_9]
M = genericMatrix(R,3,3)
I = minors(2,M)

fastSingularities = I -> (
    R := ring I;
    n := numgens R;
    gensR := gens R;
    singlist := {};
    for i from 0 to (n-1) do {
	affineChart := I + ideal(gensR_i - 1);
	sing        := radical ideal mingens ideal
	               singularLocus minimalPresentation affineChart;
	inv         := affineChart.cache.minimalPresentationMap;
	singlist = singlist | {(homogenize(preimage(inv,sing),gensR_i))};
	};
    saturate intersect(singlist)
    )

time fastSingularities I
time radical ideal singularLocus I

S = QQ[x_1..x_6,z_1..z_6,y]
M1 = matrix{{y,x_1,x_2},{x_4,y,x_3},{x_5,x_6,y}}
M2 = matrix{{y,z_1,z_2},{z_4,y,z_3},{z_5,z_6,y}}
J = minors(2,M1) + minors(2,M2)

time fastSingularities J
time  ideal singularLocus J
