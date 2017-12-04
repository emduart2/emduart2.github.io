--

needsPackage "ChainComplexExtras"

-- Exercise problems for Day 3. Chapter 4 Miller-Sturmfels.

a2c a2b ab2 b2c bc2 ac2



-- Example 4.8
-- Go over the construction of the cellular resolution using the octahedron.
-- Is the minimal resolution of this ideal supported by a cellular resolution?
-- If so, which one?

R=QQ[x1,x2,x3,x4]

I = ideal(x1*x2, x1*x3, x1*x4, x2*x3, x2*x4, x3*x4)
resol = res I
resol.dd
syz gens I

-- For the ideal in example 4.8, compare its Taylor resolution with
-- the Koszul complex of thhe generators of I. Is the Koszul complex exact in this case?
-- what's the difference between these complexes.
I = monomialIdeal(x1*x2, x1*x3, x1*x4, x2*x3, x2*x4, x3*x4)
tay = taylorResolution I
kos = koszul gens I

kos.dd
tay.dd
-- For what kind of monomial ideals are these two complexes the same?

--=====================================================================-------------------
-- 4.13 Describe the hull resolution of the ideal  
-- <x1, x2, x3, x4>^m , m=1,2,3 and compare it with the Eliahou–Kervaire resolution.
--=====================================================================-------------------
-- Recall the code to construct an Eliahou-Kervaire resolution
-- for the initial of the syzygy module

restart
R=QQ[x1,x2,x3,x4,MonomialOrder=>Lex]

I=ideal(x1*x2*x4^4,x1*x2*x3*x4^2,x1*x3^6,x1*x2*x3^2,x2^6,x1*x2^2,x1^2) --- this is the ideal from Day 2.

I= (ideal(x1,x2,x3,x4))^3

matrixOfsyz = syz gens I -- there are 12 syzygies
leadingterms=leadTerm matrixOfsyz
--leadingterms= leadTerm gens gb  image matrixOfsyz -- take leading terms of gb for syzygies
leadingterms_{1}
modulebasis=apply(rank target leadingterms, i->e_(i+1)) -- The e_i represent a module basis for R^7

S=R[modulebasis] -- take ring over R and the new variables e_i that represent a module basis
modulebasis=matrix{apply(rank target leadingterms,i->e_(i+1))}
leadTermModulebasis=modulebasis*leadingterms -- write leading terms of the grobner basis of syzygies in terms
-- of the module basis e_i
leadingVectors=apply(rank source modulebasis,i->mingens ideal(e_(i+1)*contract(e_(i+1),leadTermModulebasis)))
-- WARNING: the leadingVectors list has first entry zero so we drop it.
leadingVectors=drop(leadingVectors,{0,0})
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
--=====================================================================-------------------


--=====================================================================-------------------
-- Computations with polyhedra in M2.
-- Use package polyhedra to compute Hull resolutions.
restart
needsPackage "Polyhedra"
V = matrix {{0,2,-2,0},{-1,1,1,1}}
P = convexHull V
vertices P
fVector P
-- Example 4.21 with t=2 and n=4
V= matrix{{2,1,1,1},{1,2,1,1},{1,1,2,1},{1,1,1,2}}
P1= convexHull V
fVector P1
L=faces(0,P1) -- access the faces of codimension 0
faces(1,P1)
faces(2,P1)
faces(3,P1)
--=====================================================================-------------------

-- 4.16 Describe an algorithm for computing the hull complex of a given monomial 
-- ideal. Analyze the running time of your algorithm.

-- First examples of Alexander duality via permutohedra and tree ideals 
restart

R=QQ[x,y,z]
I=monomialIdeal(x*y^2*z^3, x*y^3*z^2, x^2*y*z^3, x^2*y^3*z, x^3*y*z^2, x^3*y^2*z)
res I
tay= taylorResolution I
Ialexdual = dual I
res Ialexdual
res I

