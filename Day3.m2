
restart
R=QQ[x,y]


--- Mahsa's example for monomial ideals in two variables

I=monomialIdeal(x^5*y,x^3*y^2,x^2*y^3,y^4)

R = QQ[x,y];
I = ideal (x^5*y, x^3*y^2,x^2*y^3,y^4)
--decompose is not relevant to what we want! this gives the Minimal associated Primes! :) 
decompose I
--Again not we want! :)
primaryDecomposition I
primaryDecomposition(I, Strategy => Monomial)
-- Using Prop 3.2 we obtain an irreducible decomposition ode thhe ideal
intersect(ideal(y),ideal(x^5,y^2),ideal(x^3,y^3),ideal(x^2,y^4))


--Checking Wether an I deal is Strongly Generic: 
L = flatten entries gens I
S = subsets (L , 2) 
for i from 0 to #S -1 do 
    (for j from 0 to #L-1 do 
	(if exponents(lcm (S#i) > exponents L#j then )
	) 
    
    ;



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
use R
J=polarization(I)
describe R
res J
codim J
use S
res I
codim I
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

-- Some of the coefficients


--Checking Wether an I deal is Strongly Generic: 
L = flatten entries gens I
S = subsets (L , 2) 
for i from 0 to #S -1 do 
    (for j from 0 to #L-1 do 
	(if exponents(lcm (S#i) > exponents L#j then )
	) 
    
    ;
   




