
-- Day 2 . 23.10.2017. Syzygies Reading group
-- Chapter 2 - Borel fixed monomial ideals
--


installPackage("GenericInitialIdeal")
installPackage("LexIdeals")

S=QQ[a,b,c,d]
--Note: since we didn't specify a monomial order, Macaulay2 will use the RevLex order with a>b>c>d.
--A lot of the following code is taken from Eisenbud-Grayson-Stillman-Sturmfels 'Computations in Algebraic Geometry with Macaulay2'
--Testing whether a monomial ideal is Borel-fixed (with respect to the obvious ordering of variables.)
isBorel monomialIdeal(a^2,a*b,b^2)
isBorel monomialIdeal(a^2,b^2)
isBorel monomialIdeal(b^3)
--'borel' gives the smallest Borel-fixed ideal containing a given monomial ideal.
borel monomialIdeal(b*c)
borel monomialIdeal(a,c^3)

--We can compute generic initial ideals with the 'GenericInitialIdeal' package:
needsPackage("GenericInitialIdeal")
I = ideal(b^2,c*d)
ginI = gin I

--The following code defines "by hand" a probabilistic algorithm to compute generic initial ideals (works with probability 1).
--(This is superfluous now since we have a package doing the work for us.)
mygin = method();
mygin Ideal := I -> (
    S := ring I;
    StoS := map(S, S, random(S^{0}, S^{numgens S:-1}));
    monomialIdeal StoS I);
mygin MonomialIdeal := I -> mygin ideal I;
mygin I

--Hilbert Series and (Castelnuovo-Mumford-)regularity stay the same when taking  passing through initial ideals, but the Betti numbers go up.
--(Remark: in all computations below, Macaulay2 computes the Hilbert polynomials, Betti numbers,... for the module S/I, not for the module I)
--K-Polynomials:
poincare I
poincare ginI
--Betti numbers and regularity:
res coker gens I
resI = res I
resginI = res ginI
betti resI
betti resginI
regularity I
regularity ginI

--With 'lexgin', we can compute the generic initial ideal with respect to the lex-order:
lexginI = lexgin I
poincare lexginI
--Here the regularity goes up (this can already be seen from the fact that there is a degree 4 generator).
betti res lexginI
regularity lexginI

--Example 2.10 in Miller-Sturmfels. 
genericForms = (p,q) ->  ideal(random(p,S), random(q,S));
gin genericForms(2,2)
gin genericForms(2,3)
lexgin genericForms(2,2)
lexgin genericForms(2,3)
--You can play around here by filling in various numbers for p and q and computing gin, lexgin, Betti tables,...

--Example of a non-generic initial ideal that's also Borel-fixed:
J=ideal(a^2,a*b+b^2,a*c)
ginJ=monomialIdeal gin J
inJ=monomialIdeal J
isBorel inJ and isBorel ginJ


-- Module Grobner bases


-- Example 2.19 - Construction of the minimal free resolution for in(M) where M is the module of syzygies for the monomial
-- ideal I.

restart
R=QQ[x1,x2,x3,x4,MonomialOrder=>Lex]

I=ideal(x1*x2*x4^4,x1*x2*x3*x4^2,x1*x3^6,x1*x2*x3^2,x2^6,x1*x2^2,x1^2)

matrixOfsyz = syz gens I -- there are 12 syzygies
leadingterms= leadTerm gens gb  image matrixOfsyz -- take leading terms of gb for syzygies
leadingterms_{1}
modulebasis=apply(rank target leadingterms, i->e_(i+1)) -- The e_i represent a module basis for R^7

S=R[modulebasis] -- take ring over R and the new variables e_i that represent a module basis
modulebasis=matrix{apply(rank target leadingterms,i->e_(i+1))}
leadTermModulebasis=modulebasis*leadingterms -- write leading terms of the grobner basis of syzygies in terms
-- of the module basis e_i
leadingVectors=apply(rank source modulebasis,i->mingens ideal(e_(i+1)*contract(e_(i+1),leadTermModulebasis)))
-- collect all the terms that have the leading term in the same position e_i.
use S
-- we construct a resolution of the leadingTerms of the module of syzygies by using Koszul complexes
 koszulComplexes=apply(#leadingVectors-1,j->koszul(leadingVectors_j))-- compute the Koszul complexes of the elements
 -- in leadingVectors
resLeadSyzygiesI=fold((i,j) -> i++j, koszulComplexes) -- add the koszul complexes of the leadingVectors, this gives the resolution
-- of the module of initial terms of the syzygies, note that there are 6 initial terms
resI = res coker gens I
resI
resLeadSyzygiesI
resI.dd
resLeadSyzygiesI.dd
-- compare the ranks of the modules in resI and resLeadSyzygiesI, also compare their matrices.

installPackage("MonomialIdealResolutions.m2")
input "MonomialIdealResolutions.m2"


J=monomialIdeal(x1*x2*x4^4,x1*x2*x3*x4^2,x1*x3^6,x1*x2*x3^2,x2^6,x1*x2^2,x1^2)
EKResolution J
o83.dd
