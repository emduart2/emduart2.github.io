
--- Day 1 -- Syzygies Reading group 
--- Chapter 1 - 'Combinatorial Commutative algebra' Miller and Sturmfels
--- Some M2 computations with Hilbert functions, Hilbert series,
--- Chain complexes free resolutions, Betti tables and simplicial homology.
--- Koszul complex, Koszul homology,


-- 
restart


installPackage("SimplicialComplexes")
needsPackage("SimplicialComplexes")
help "SimplicialComplexes"

-- Singly graded examples
-- Example 1 on four vertices {1,2,3,4}, variables are x1,x2,x3,x4
-- 'Delta' denotes the simplicial complex.
-- Delta={{},all subsets of size two, all vertices, empty face}
-- This is a N-graded example. Also refered as singly graded

R=QQ[x1,x2,x3,x4]
I=monomialIdeal(x2*x3*x4,x1*x3*x4,x1*x2*x4) -- Stanly Reisner ideal associated to the non faces of Delta
hilbertFunction(3,I)  -- Hilbert function of ideal I in degree 3

hilbertSeries(I) -- N-graded Hilbert Series of I
help "poincare" -- this type of command shows the documentation for the function 'poincare'
poincare(I)     -- K-polynomial of I
hilbertSeries(R) -- N-graded Hilbert Series of R
poincare(R) -- K-polynomial of R
degree(x2*x3*x4) -- check the degree of a monomial

help "SimplicialComplexes" -- Documentation for the SimplicialComplexes package
help "simplicialComplex"
CDelta= simplicialComplex I -- Takes monomial ideal corresponding to nonfaces and returns a list of monomials
                             -- corresponding to the faces of Delta
Lfacets={x1*x2*x3,x1*x4,x2*x4,x3*x4}
simplicialComplex Lfacets
decompose I -- outputs prime decomposition of the ideal I
boundary(CDelta) -- topological boundary of Delta
homCDelta= prune HH CDelta -- computes simplicial homology of the complex associated to Delta

DeltaStar=dual(CDelta)
ideal DeltaStar
decompose I
intersect(ideal(x1,x2,x4),ideal(x1,x3,x4),ideal(x2,x3,x4)) -- prime decomposition of the Alexander dual 
-- following the book.

minres=res coker gens I -- computes minimal free resolution of I
minres.dd   -- we take a peek to thhe matrices appearing in the minimal free resolution
n=pdim coker gens ideal(x2*x3*x4,x1*x3*x4,x1*x2*x4) -- computex the projective dimension i.e maximal length of
-- a minimal free resolution
help "pdim"
apply(n+1,i->tally degrees source minres.dd_i)

betti minres -- Betti table corresponding to the minimal free resolution of I


-- Multigraded examples


restart
installPackage"SimplicialComplexes"
R=QQ[x1,x2,x3,x4,x5,Degrees=>{{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}}] -- multigraded ring

I=monomialIdeal(x1*x5,x4*x5,x2*x3*x5,x1*x2*x3*x4)

minres=res coker gens I  -- minimal resolution

n=pdim coker gens I      -- projective dimension

apply(n+1,i->tally degrees source minres.dd_i) -- tally the degrees of the source of the matrices in the complex
-- These are thhe multigraded Betti numbers

betti minres


help "koszul"
KK=koszul gens ideal(x1,x2,x3,x4) -- Koszul complex of the sequence {x1,x2,x3,x4}
KK.dd  -- matrices of the koszul complex
koszulHomology= HH(KK)
prune koszulHomology -- The  homology of the koszul complex is trivial, therefore the koszul complex
                     -- is a free resolution of ideal(x1,x2,x3,x4) 

KK=koszul gens I      -- we can compute the koszul complex of any sequence of polynomials, this
KK.dd                 -- complex is not necesarily exact as we can see in this example
koszulHomology=HH(KK)
prune koszulHomology

-- compare minres vs. KK
minres
KK

-- We can also compute Betti numbers using the definition of Tor
--
degrees Tor_0(module R/I,module R/ideal(x1,x2,x3,x4,x5))
degrees Tor_1(module R/I,module R/ideal(x1,x2,x3,x4,x5))
degrees Tor_2(module R/I,module R/ideal(x1,x2,x3,x4,x5))
degrees Tor_3(module R/I,module R/ideal(x1,x2,x3,x4,x5))


-- Exercise 1 Ch 1 - Miller-Surtmfels

-- Let n=6 and Delta be the boundary of the octahedron.
--(a) Determine I_Delta and I_Delta^*
--(b) Compute their Hilbert series
--(c) Compute their minimal free resolutions
--(d) Interpret the Betti numbers in (c) in terms of simplicial homology
