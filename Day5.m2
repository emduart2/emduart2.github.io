

-- Day 5. Syzygies reading group. 
-- Sections 2B, 2C 'Geometry of syzygies' by Eisenbud.

-- Finding the Hilbert series of points in projective spaces.
-- We can find it by either using the ideal of the points,
-- or the points themselves.

-- Examples of a module that has the given hilbert function after proposition 2.7

restart
S=QQ[x0,x1,x2]

I=ideal(x0^2,x1^2,x2^2,x0*x1*x2)

betti res I

restart
S=QQ[x,y,z,w]
I=ideal(x,y,z) -- V(I)={[0;0;0;1]}

L={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}} -- list of points in PP^3

-- To evaluate the Hilbert function of points in degree d=2 using the points themselves,
-- compute a basis of the polynomial ring in degree 2 and evaluate this basis at each point
B  = super basis({2},S)

rows =  apply(#L,i->sub(B, matrix{L_i}))
M = fold((i,j)->i||j,rows)
rank M
-- The hilbert function of the points is the rank of the matrix M.

-- Come up with two sets of 7 points in PP^3 in general position, one whose points lie on a twisted cubic
-- and one whose points are not on a twisted cubic.


-- Prove Bezout's theorem for two plane curves of degrees a,b with states that the
-- two curves intersect on exactly a*b points by using a free resolution of the ideal of the points.
-- USe the free resolution to compute the Hilbert polynomial.


r=3
S=QQ[x_0..x_r]

--md denotes a vector whose entries are the degree d monomials of S.
m1=vars S
m2=monomials (m1**m1)
m3=monomials (m2**m1)
m4=monomials (m3**m1)

--define function evaluate which evaluates a polynomial at a point; UNFINISHED!!
evaluate = method()
evaluate (RingElement,Sequence) := RingElement => (f,s) ->(
    R=ring f;
    mu1=monomials R;
    if #mu1 != #s then( print "ERROR: the number of variables is different from the number of coordinates of the point"; break;)
Fmin)

-- choose 7 in general position in \PP^3
p1=(1,0,0,0)
p2=(0,1,0,0)
p3=(0,0,1,0)
p4=(0,0,0,1)
p5=(1,1,1,0)
p6=(1,1,0,1)
p7=(0,1,1,1)

--   where by evaluate(m1,p1) we mean 
--   sub(m1,x_0=>p1_0, x_1=>p1_1, x_2=>p1_2, x_3=>p1_3)

M1 = evaluate(m1,p1)  || evaluate(m1,p2) || evaluate(m1,p3) || evaluate(m1,p4) || evaluate(m1,p5) || evaluate(m1,p6) || evaluate(m1,p7)
M2 = evaluate(m2,p1)  || evaluate(m2,p2) || evaluate(m2,p3) || evaluate(m2,p4) || evaluate(m2,p5) || evaluate(m2,p6) || evaluate(m2,p7)
M3 = evaluate(m3,p1)  || evaluate(m3,p2) || evaluate(m3,p3) || evaluate(m3,p4) || evaluate(m3,p5) || evaluate(m3,p6) || evaluate(m3,p7)
M4 = evaluate(m4,p1)  || evaluate(m4,p2) || evaluate(m4,p3) || evaluate(m4,p4) || evaluate(m4,p5) || evaluate(m4,p6) || evaluate(m4,p7)

rank M1
rank M2
rank M3
rank M4

--prove that the twisted cubic has two linear syzygies