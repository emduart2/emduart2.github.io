
--- Day 1 -- Syzygies Reading group
--- Some M2 computations with Hilbert functions, Hilbert series,
--- Chain complexes free resolutions, Betti tables and simplicial homology.
--- Koszul complex, Koszul homology,


-- 
restart


installPackage("SimplicialComplexes")
needsPackage("SimplicialComplex")
help "SimplicialComplexes"
R=QQ[x1,x2,x3,x4, Degrees=>{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}]
I=monomialIdeal(x2*x3*x4,x1*x3*x4,x1*x2*x4)
degree(x2*x3*x4)

simplicialComplex 
Lfacets={x1*x2*x3,x1*x4,x2*x4,x3*x4}
decompose I

hilbertFunction(3,I)
hilbertSeries(I) -- N-graded Hilbert Series of I
hilbertSeries(R)
help "poincare"
poincare I -- K-polynomial of I
poincare R -- K-polynomial of R


help "simplicialComplex"
CDelta = simplicialComplex I
Cdelta.dd

homCDelta= HH(CDelta)
prune hhomCDelta
simplicialComplex Lfacets
n=pdim ideal(x2*x3*x4,x1*x3*x4,x1*x2*x4)
view "pdim"
apply(,i->tally degree source o68.dd_i)

minres=res coker gens I

minres.dd
betti minres

help "hilbertFunction"

HH(o14,-1)
prune(o15)
help "prune"

decompose I

octahedron


-- Exercise 1 Ch 1 - Miller-Surtmfels

-- Let n=6 and Delta be the boundary of the octahedron.
--(a) Determine I_Delta and I_Delta^*
--(b) Compute their Hilbert series
--(c) Compute their minimal free resolutions
--(d) Interpret the Betti numbers in (c) in terms of simplicial homology





help "koszul"
KK=koszul gens ideal(x1,x2,x3,x4)
KK.dd
koszulHomology= HH(KK)
prune koszulHomology

KK=koszul gens I
KK.dd
koszulHomology=HH(KK)
prune koszulHomology

R=QQ[x1,x2,x3,x4, Degrees=>{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}]
I=ideal(x2*x3*x4,x1*x3*x4,x1*x2*x4)
minres=res coker gens I
n=pdim coker gens I
apply(n,i->tally degrees source minres.dd_i)
betti minres

help "pdim"
help "eagonNorthcott"

degrees source minres.dd_1
