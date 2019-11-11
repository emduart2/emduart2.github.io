needsPackage "FourTiTwo"
S=QQ[p00,p01,p10,p11]
A= matrix{{1,1,0,0},{0,0,1,1},{1,0,1,0}}
toricMarkov(A,S)
toricMarkov(A)
viewHelp "MarkovBasis" 	  


restart
needsPackage "GraphicalModels"

R= markovRing(2,2,2) 
gens R

C1= {{{2},{3},{}}}
I1 = conditionalIndependenceIdeal(R,C1)
quadric1 = ideal ((p_(1,1,1)+p_(2,1,1))*(p_(1,2,2)+p_(2,2,2))-(p_(1,1,2)+p_(2,1,2))*(p_(1,2,1)+p_(2,2,1)))
quadric1 == I1

C2={{{1},{2},{3}}}
I2= conditionalIndependenceIdeal(R,C2)

C=C1|C2

IC= conditionalIndependenceIdeal(R,C)
I1+I2 == IC

comp =  decompose IC
netList comp

C3={{{2},{1,3},{}}}
I3=conditionalIndependenceIdeal(R,C3)

comp#2 == I3
