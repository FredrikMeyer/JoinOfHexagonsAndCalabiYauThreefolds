restart
loadPackage "VersalDeformations"
kk = ZZ/3001
R = kk[x_1..x_12]

generateX2 = () -> (
    K = random(R^18,R^12);
    a = transpose gens gb K; -- samme bilde
    b = entries a;
    b = apply(0..11, i-> apply(b#i, z -> z*x_(i+1)));
    bb = sum toList b;
    bb1 = bb_{0..8};
    bb2 = bb_{9..17};
    M1 = matrix toList apply(0..2, i-> toList apply(0..2, j-> bb1#(3*i+j)));
    M2 = matrix toList apply(0..2, i-> toList apply(0..2, j-> bb2#(3*i+j)));
    I1 = minors(2, M1);
    I2 = minors(2, M2);
    I1+I2
    )

IX = generateX2()
IX = oo




--apply(1..12, i-> rank source gens gb minimalPresentation (IX + (x_i-1)))


T1 = time CT^1(0, gens IX);
source T1
--- fikk 39, etter 3 timer, karakteristikk 5
--- fikk 39 på nytt, karakteristikk 7, etter 3.3 timer
--- 39 nå også ,i kar 11, etter 7 timer