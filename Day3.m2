
restart
R=QQ[x,y]


--- Mahsa's example for monomial ideals in two variables

I=monomialIdeal(x^5*y,x^3*y^2,x^2*y^3,y^4)

R = QQ[x,y];
I = ideal (x^5*y, x^3*y^2,x^2*y^3,y^4)
--decompose is not relevant to what we want! this gives the Minimal associated Primes! :) 
decompose I
--Again not what we want! :)
primaryDecomposition I
primaryDecomposition(I, Strategy => Monomial)
-- Using Prop 3.2 we obtain an irreducible decomposition of the ideal
intersect(ideal(y),ideal(x^5,y^2),ideal(x^3,y^3),ideal(x^2,y^4))


-- Section 3.2 Examples in 3 variables
-- We define the polarization function
restart
S=QQ[x,y,z]
polarization = (I) ->(
    n = numgens ring I;
    u = apply(numgens I, i->first exponents I_i);
    Ilcmm = max \ transpose u;
    Z = flatten apply(n, i-> apply(Ilcmm#i,j->w_{i,j}));
    R = QQ(monoid[Z]);
    Z= gens R;
    p=apply(n,i->sum((Ilcmm)_{0..i-1}));
    monomialIdeal apply(u,e->product apply(n,i->product(toList(0..e#i-1),j->Z#(p#i+j))))
)

I=ideal(x^4,y^4,z^4,x^3*y^2*z,x*y^3*z^2,x^2*y*z^3)
primaryDecomposition(I,Strategy=>Monomial)
use R --  The ring R is the ring S~ in which we take the Stanley-Reisner ideal
J=polarization(I)
describe R
res J
codim J
use S
res I
codim I

--------------------------------------------

--Checking whether an ideal is Strongly Generic: 

IsStronglyGeneric = I ->(
    L = flatten entries gens I;
    exps=flatten apply(numgens I, i->exponents I_i);
    S1 = apply(rank source vars ring I,j->apply(#exps,i-> (exps_i)_j));
    S2=apply(#S1,i->delete(0,S1_i));
    S3 = apply (#S2,i->unique(S2_i));
    flatten S2 == flatten S3
)  

--test 2 examples of the book section 3.2 and 3.3:
Q = QQ[x,y,z];
J = ideal (x^4,y^4,z^4,x^3*y^2*z,x*y^3*z^2,x^2*y*z^3);
IsStronglyGeneric(J)
I' = ideal (x^2*z,x*y*z,y^2*z,x^3*y^5,x^4*y^4,x^5*y^3);
IsStronglyGeneric(I')

--------------------------------------------

-- Chapter 4
-- Example of a cellular resolution.

restart
S=QQ[a,b,c,d]
I=ideal(a^2*b,a^2*d,a*b^2,a^2*c,a*d^2,b^2*d,b*d^2,b^2*c,a*c^2,c*d^2,b*c^2,c^2*d)
J=polarization(I)
res J
res o2


--------------------------------------------


primaryDecomposition(J, Strategy=>Monomial)

-- Example of an ideal that whose minimal free resolution is not supported on
-- a CW-complex
restart
R=QQ[x1,x2,x3,x4,x5,x6,xv]
I=ideal(x3*x4*x5*x6,x2*x4*x5*x6,x1*x4*x5*x6,x3*x4*x5*xv,x2*x4*x5*xv,x1*x3*x5*xv,x1*x2*x4*xv,x1*x4*x6*xv,
    x1*x5*x6*xv,x3*x4*x6*xv,x2*x5*x6*xv,x2*x3*x6*xv,x1*x2*x3*xv)
resolute=res coker gens I
resolute.dd

-- Some of the coefficients in the monomials of the matrices of this resolution have
-- a coefficient that is 2. This means that the resolution cannot be obtained as the ceullar homology
-- of a CW - complex.
