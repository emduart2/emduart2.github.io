
--- Day 1 -- Syzygies Reading group
--- Some M2 computations with Hilbert functions, Hilbert series,
--- Chain complexes free resolutions, Betti tables and simplicial homology.
--- Koszul complex, Koszul homology,


-- 
restart


installPackage("SimplicialComplexes")
needsPackage("SimplicialComplex")
help "SimplicialComplexes"

R=QQ[x1,x2,x3,x4]
I=monomialIdeal(x2*x3*x4,x1*x3*x4,x1*x2*x4)
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
homCDelta= HH(CDelta)
prune hhomCDelta
simplicialComplex Lfacets


minres=res coker gens I
minres.dd
betti minres

help "hilbertFunction"

HH(o14,-1)
prune(o15)
help "prune"

decompose I



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

help "eagonNorthcott"

