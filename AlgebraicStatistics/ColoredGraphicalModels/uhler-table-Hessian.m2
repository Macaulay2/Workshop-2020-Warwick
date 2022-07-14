-- GRAPH 1
restart
-- initialize all functions in the beginning of the file, then run this code
n=3
R=QQ[l_1..l_n]
K=matrix{{l_1,l_2,0,l_2},{l_2,l_1,l_3,0},{0,l_3,l_1,l_2},{l_2,0,l_2,l_1}}
f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI --2
M=cokernel jacobian ideal f
C=resolution M
B=betti C
H=jacobian ideal flatten entries jacobian ideal f
rank H --3


-- GRAPH 7
restart
-- initialize all functions in the beginning of the file, then run this code
n=5
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_3},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_3,0,l_5,l_1}}
f=det K
gradI= ideal flatten entries jacobian ideal f
M=cokernel jacobian ideal f
C=resolution M
B=betti C
codim gradI --2
H=jacobian ideal flatten entries jacobian ideal f
rank H --5


-- GRAPH 9
restart
-- initialize all functions in the beginning of the file, then run this code
n=6
R=QQ[l_1..l_n]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
f=det K
gradI= ideal flatten entries jacobian ideal f
M=cokernel jacobian ideal f
C=resolution M
B=betti C
codim gradI --3
H=jacobian ideal flatten entries jacobian ideal f
rank H--6

-- GRAPH 16
restart
-- initialize all functions in the beginning of the file, then run this code
n=7
R=QQ[l_1..l_n]
K=matrix{{l_1,l_4,0,l_7},{l_4,l_1,l_5,0},{0,l_5,l_2,l_6},{l_7,0,l_6,l_3}}
f=det K
gradI= ideal flatten entries jacobian ideal f
M=cokernel jacobian ideal f
C=resolution M
B=betti C
codim gradI --codim 3
H=jacobian ideal flatten entries jacobian ideal f
rank H --7


-- GRAPH 11
restart
n=5
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_4,0,l_5,l_2}}
f=det K
gradI= ideal flatten entries jacobian ideal f
M=cokernel jacobian ideal f
C=resolution M
B=betti C
codim gradI --2
H=jacobian ideal flatten entries jacobian ideal f
rank H --5

-- GRAPH 13
restart
n=6
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_6,0,l_5,l_2}}
f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI --3
M=cokernel jacobian ideal f
C=resolution M
B=betti C
H=jacobian ideal flatten entries jacobian ideal f
rank H --6

-- GRAPH 17
restart
n=7
R=QQ[l_1..l_n]
K=matrix{{l_1,l_4,0,l_7},{l_4,l_2,l_5,0},{0,l_5,l_1,l_6},{l_7,0,l_6,l_3}}

f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI --3
M=cokernel jacobian ideal f
C=resolution M
B=betti C
H=jacobian ideal flatten entries jacobian ideal f
rank H -- 7

-- GRAPH 18
restart
n=8
R=QQ[l_1..l_n]
K=matrix{{l_1,l_5,0,l_8},{l_5,l_2,l_6,0},{0,l_6,l_3,l_7},{l_8,0,l_7,l_4}}
f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI --3
M=cokernel jacobian ideal f
C=resolution M
B=betti C
H=jacobian ideal flatten entries jacobian ideal f
rank H --8 

