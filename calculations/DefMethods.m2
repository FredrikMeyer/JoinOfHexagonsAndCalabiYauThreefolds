-- *- coding: utf-8 -*-

newPackage(
        "DefMethods",
        Version => "1.0", 
        Date => "August 22, 2012",
        Authors => {{Name => "Fredrik Meyer", 
                  Email => "fredrme@math.uio.no", 
                  HomePage => "http://folk.uio.no/fredrme/"}},
        Headline =>   " - some methods that should have been standard",
        DebuggingMode => false,
	AuxiliaryFiles => false,
	Headline => "Some usueful methods"
        )

export {"oneMatrix", "findT2", "contains",
     "makeMonomial", "changeMatrix", "makeCone", "matrix2Sage", "toBinomial",
     "toricIdeal","endTerm", "makeBinomial", "polarization" ,  "distraction", "isMonomial",
     "isRadical", "projEmb","pp", "fastSingularities"}


-- returns the end term of f,
-- i.e. the term of least degree
endTerm = method()
endTerm(RingElement) := (f) -> (
    g := f;
    while (f != 0) do (
	g  = f;
	f  = f - leadTerm(f);
	);
    g
    )

--- transform a ring polynomial element
--- into a binomial
makeBinomial = method()
makeBinomial(RingElement) := (f) -> (
    while (#(terms f) > 2) do (
	f = f - endTerm f;
	);
    f
    )

makeBinomial(Ideal) := (I) -> (
    L := flatten entries gens gb I;
    ideal apply(L, f -> makeBinomial f)
    )

    
isMonomial = method()
isMonomial(Ideal) := (I) -> (
    F := mingens I;
    val := true;
    for i from 0 to (numcols F -1) do (
	if #terms(F_i_0) > 1 then (
	    val = false;
	    break;
	    );
	);
    val
    )

-- checks if a monomialideal is radical
isRadical = method()
isRadical(Ideal) := (I) -> (
    F := mingens I;
    val := true;
    for i from 0 to (numcols F - 1) do (
	if any((exponents(F_i_0))_0, i-> (i > 1)) then (
	    val = false;
	    break;
	);
    );
    val
    )

-- computes the polarization of a
-- monomial ideal
polarization = method()
polarization(Ideal) := I -> (
    n := numgens ring I;
    u := apply(numgens I, i -> first exponents I_i);
    Ilcm := max \ transpose u;
    PP := flatten apply(n, i -> apply(Ilcm#i, j -> pp_(i,j)));
    R := QQ(monoid[PP]);
    PP = gens R;
    p := apply(n, i -> sum((Ilcm)_{0..i-1}));
    monomialIdeal apply(u, e -> product apply(n, i ->product(toList(0..e#i-1), j -> PP#(p#i+j))))
    )

-- computes the distraction of a monomial ideal
distraction = method()
distraction(Ideal) := I -> (
    S := ring I;
    n := numgens S;
    X := gens S;
    J := polarization I;
    u := apply(numgens I, i -> first exponents I_i);
    Ilcm := max \ transpose u;
    W := flatten apply(n, i -> flatten apply(Ilcm#i, j -> X#i));
    section := map(S, ring J, apply(W, r -> r -
	    random(500)*X#(n-2) - random(500)*X#(n-1)));
    section ideal J
    )	

--- change the matrix M in position (m,n) to the
--- element r.
changeMatrix = method()
changeMatrix(Matrix,ZZ,ZZ,Thing) := (M,m,n,r) -> (
     Mm := mutableMatrix M;
     Mm_(m,n) = r;
     M = matrix Mm
     )
 
--- a matrix of ones
oneMatrix = method()
oneMatrix(ZZ, ZZ) := (m,n) -> (
    matrix toList apply(1..n, j -> toList apply(1..m, i-> 1))
     )

fastSingularities = method()
fastSingularities(Ideal) := I -> (
    R := ring I;
    n := numgens R;
    gensR := gens R;
    singlist := {};
    for i from 0 to (n-1) do {
	affineChart := I + ideal(gensR_i - 1);
	sing        := radical ideal mingens ideal singularLocus minimalPresentation affineChart;
	inv         := affineChart.cache.minimalPresentationMap;
	singlist = singlist | {(homogenize(preimage(inv,sing),gensR_i))};
	print(i);
	};
    saturate intersect(singlist)
    )

--- compute the t2 module
findT2 = method()
findT2(Ring, Ideal) := (P,I) -> (
     C := res(I, LengthLimit => 2);
     Rel0 := koszul(2, C.dd_1);
     A := P/I;
     HomR := Hom(image(C.dd_2)/image Rel0**A, A);
     Triv := image substitute(transpose C.dd_2, A);
     T2 := HomR/Triv;
     B := matrix basis(0, T2);
     T2Mat := (gens HomR) * B
     )

--checks if a list contains an element
-- equivalent to member(t,L), but better
-- name
contains = method() 
contains(BasicList,Thing) := (L,t) -> (
    member(t,L)
     )

--- given a 1xm matrix (m=numgens R),
--  make a monomial with exponents given by
-- M
makeMonomial = method()
makeMonomial(Matrix, PolynomialRing) := (M, R) -> (
     G := gens ideal gens R;
     p := 1_R;
     for z from 0 to (numColumns(G)-1) do (
	  if (M_z_0) < 0 then (
	       p = p / ((G_z_0)^(-M_z_0));
	       )
	  else (
	       p = p*((G_z_0)^(M_z_0));
	       )
	  );
     p
     )

-- given a prime binomial ideal
-- make a cone with the same
-- toric variety
makeCone = method()
makeCone(Ideal) := (I) -> (
     mGens := mingens I;
     M  := {};
     lll := {};
     for j from 0 to (numColumns(mGens)-1) do (
     	  lll = exponents (mGens)_j_0;
	  if (length lll == 1) then (
	       M = M | {lll#0};
	       )
	  else if (length lll == 2) then (
	       M = M | {(lll#0-lll#1)};
	       )
     	  );
     B := transpose matrix M;
     transpose LLL syz  matrix transpose B
     )

-- given a divisor on a toric variety, 
-- make an embedding into projective space 
projEmb = method()
projEmb(Thing) := (D) -> (
    X := variety D;
    L := OO D;
    m := rank HH^0(X,L);
    S := ring X;
    R := QQ[vars(0..m-1)];
    phi := map(S,R, basis(-first degrees L, S));
    kernel phi
    )
 
toBinomial = method()
toBinomial(List, Ring) := (b,R) -> (
    topp := 1_R; bottom := 1_R;
    scan(#b, i-> if b_i > 0 then topp = topp * R_i^(b_i)
	else if b_i < 0 then bottom = bottom * R_i^(-b_i));
    topp - bottom
    )


--- given a matrix generating a cone,
--- make the ideal 
toricIdeal = method();
toricIdeal(Matrix) := (A) -> (
    n := #(transpose entries A);    
    toricR := QQ[vars(0..n-1)]; --, Degrees => transpose entries A, MonomialSize => 16];
    B := transpose LLL syz matrix A;
    J := ideal apply(entries B, b-> toBinomial(b,toricR));
    scan(gens ring J, f-> J = saturate(J,f));
    J   
    )

-- given a macaulay2 matrix
-- output a string suitable for SAGE
-- {{..}} turns into [[..]]
matrix2Sage = method();
matrix2Sage(Matrix) := (M) -> (
     s := toString M;
     s = substring((7,length toString M), toString M);
     nystr := "";
     repll := {{"{","["},{"}","]"}};
     bool := true;
     for i from 0 to (length s - 1) do (
	  bool = true;
	  for j from 0 to (length repll - 1) do (	  
	       if s#i == repll#j#0 then (
		    nystr = nystr | repll#j#1;
		    bool = false;
		    break;
	       	    )
	       );
	  if bool then (
	       nystr = nystr | s#i;
	       )	  
     );
     nystr
     )

-----------    
beginDocumentation()
document {
    Key => DefMethods,
    Headline => "Some useful methods",
    PARA {
	"
	This package contains some useful methods from different fields of algebraic geometry.
	"
	}
    }

document {
    Key => {changeMatrix, (changeMatrix,Matrix,ZZ,ZZ,Thing)},
    Headline => "Change a matrix at a position"
    }

document {
    Key => {contains,(contains,BasicList,Thing)},
    Headline => "Contains",
    PARA {
	"This is synonymous with the built-in Macaulay2-function ", TO member, "."
	}
    }

document {
    Key => {distraction, (distraction,Ideal)},
    Headline => "Distraction of a monomial ideal"
    }

document {
    Key => {endTerm, (endTerm, RingElement)},
    Headline => "The trailing term"
    }

document{
    Key => {fastSingularities, (fastSingularities,Ideal)},
    Headline => "Compute singular locus",
    PARA {"For example, let us compute the singular locus of the cusp."},
    EXAMPLE {"R = QQ[x,y,z];",
	"I = ideal(x^2*z-y^3);",
	"fastSingularities I"
	}
    }

document{
    Key => {isMonomial, (isMonomial,Ideal)},
    Headline => "Checks if an ideal is monomial"
    }

document {
    Key => {isRadical, (isRadical,Ideal)},
    Headline => "Checks if a monomial ideal is radical"
    }

document{
    Key => {makeBinomial, (makeBinomial, Ideal),(makeBinomial,RingElement)},
    Headline => "Trim an ideal to a binomial ideal"
    }

document{
     Key => {makeCone, (makeCone,Ideal)},
     Headline => "Makes the cone of a toric ideal."
     }

document{
    Key => {makeMonomial, (makeMonomial,Matrix,PolynomialRing)},
    Headline => "From exponent vector to monomial"
    }

document{
     Key => {findT2, (findT2,Ring,Ideal)},
     Headline => "Find T^2.",
     TT "findT2(P,I)", " -- Returns a matrix representing a basis for T2."
     }
 
document {
    Key => {oneMatrix, (oneMatrix, ZZ,ZZ)},
    Headline => "A matrix of ones"
    } 
 
document{
     Key => {projEmb, (projEmb,Thing)},
     Headline => "Projective embedding"
     }
 
document {
    Key => {matrix2Sage, (matrix2Sage, Matrix)},
    Headline => "Convert to SAGE format"
    }

document {
    Key => {polarization, (polarization,Ideal)},
    Headline => "Polarization of a monomial ideal"
    }

 document {
     Key => pp,
     Headline => "A dummy variable"
     }

document {
    Key => {toBinomial, (toBinomial,List,Ring)},
    Headline => "From weightvector to binomial"
    }
 
 document {
     Key => {toricIdeal, (toricIdeal, Matrix)},
     Headline => "Toric ideal from semigroup cone"
     }