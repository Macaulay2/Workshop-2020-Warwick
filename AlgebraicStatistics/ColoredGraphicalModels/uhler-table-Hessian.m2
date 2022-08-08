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

resol=resolution(gradI)
betti(resol) --no linear syzygies
resol.dd_1
resol.dd_2
syz(gens trim gradI)

-- GRAPH 7
restart
-- initialize all functions in the beginning of the file, then run this code
n=5
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_3},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_3,0,l_5,l_1}}
f=det K
gradI= ideal flatten entries jacobian ideal f
support gradI
M=cokernel jacobian ideal f
C=resolution M
B=betti C
codim gradI --2
H=jacobian ideal flatten entries jacobian ideal f
rank H --5

resol=resolution(gradI)
betti(resol) --l=2 < maximal linear rank=n-1=4
resol.dd_1
resol.dd_2
rank oo -- maximal linear rank
syz(gens trim gradI)


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
C.dd_0
C.dd_1
C.dd_2
syz M
B=betti C
codim gradI --3
H=jacobian ideal flatten entries jacobian ideal f
rank H--6

resol=resolution(gradI)
resol_0
resol.dd_0
resol.dd_1
netList gradI_*
netList flatten entries resol.dd_1
resol.dd_2
--linear syzygies correspond to the first 3 columns
--the following 4 columns are quadratic syzygies and the last 3 are cubic syzygies

rank resol.dd_2

betti(resol)
--this can be read from the fact that in column 2 we have 3,4,3
--the graded resolution looks like this:
-- ...-> S(-4)^3+S(-5)^4+S(-6)^3 -> S(-3)^6 -> I -> 0
-- so l=3 in Simis' paper
-- l=3<n-1=5 not maximal linear rank

--it can also be seen directly here: 
syz (gens gradI)

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

resol=resolution(gradI)
betti(resol) --l=2 < maximal linear rank=n-1=6
resol.dd_1
mingens gradI
resol.dd_2
rank oo -- maximal linear rank
syz(gens trim gradI)


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


resol=resolution(gradI)
betti(resol) --l=2 < maximal linear rank=n-1=4
resol.dd_1
resol.dd_2
rank oo -- maximal linear rank
syz(gens trim gradI)


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

resol=resolution(gradI)
betti(resol) --l=1 < maximal linear rank=n-1=5
resol.dd_1
resol.dd_2
rank oo -- maximal linear rank
syz(gens trim gradI)


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

resol=resolution(gradI)
betti(resol) --l=3 < maximal linear rank=n-1=6
resol.dd_1
resol.dd_2
rank oo -- maximal linear rank
syz(gens trim gradI)

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

resol=resolution(gradI)
betti(resol) --l=3 < maximal linear rank=n-1=7
resol.dd_1
resol.dd_2
rank oo -- maximal linear rank
syz(gens trim gradI)


--test whether for one of their cloned matrices we get
--maximal linear rank
restart
n=5
R=QQ[l_1..l_n]
K=matrix{{l_1,l_2,l_3},{l_2,l_4,l_5},{l_3,l_5,l_4}}
-- 3-cycle with 2 equal vertices 
-- we expect no MLE existence for n=1 because 1 out of the 2 gens of 
-- IG1 doesn't seem to change signs

--Simis' notation: m=3 (size of matrix)
--The linear rank of the gradient ideal of f is m+1 choose 2 -2
m=3
binomial(m+1,2)-2==n-1

gradI= ideal flatten entries jacobian ideal det K
mingens gradI
loadPackage "Depth"
isRegularSequence mingens gradI
codim gradI
depth(gradI,R)
dim gradI
--not even Cohen-Macaulay

resol=resolution(gradI)
betti(resol) --l=4 = maximal linear rank
resol.dd_2
rank oo -- maximal linear rank but no linear presentation

--assume all 3 vertices equal

n=4
R=QQ[l_1..l_n]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_1}}
-- 3-cycle with 3 equal vertices 
-- p>0 for

m=3
binomial(m+1,2)-2==n-1 -- of course it's false: same m as before and different n
numgens gradI --but again n=mu(gradI)

gradI= ideal flatten entries jacobian ideal det K
codim gradI
dim gradI

resol=resolution(gradI)
betti(resol) --l=0 < maximal linear rank=n-1=3
resol.dd_2
rank oo 

--apparently the maximial linear rank condition only holds for
--two cloned entries in the main diagonal

--test whether for one of their strategic sparse matrices we get
--maximal linear rank

n=13
R=QQ[l_1..l_n]
m=6
r=4
m-r+2
m-2
K=matrix{{l_1,l_2,l_3,l_4,l_5,l_6},{l_2,l_7,l_8,l_9,l_10,0},
    {l_3,l_8,l_12,l_13,0,0},{l_4,l_9,l_13,0,0,0},
    {l_5,l_10,0,0,0,0},{l_6,0,0,0,0,0}}

gradI= trim ideal flatten entries jacobian ideal det K
codim gradI
dim gradI

resol=resolution(gradI)
betti(resol) --l=2 = maximal linear rank=mu(gradI)-1
resol.dd_2
rank oo --maximal linear rank and linearly presented



--QUESTION: n seems to be always the minimal number of 
--generators of the ideal in the colored 4-cycle 
--(and in the 3-cycle example too) ...does it make any sense?
-- it's clearly not always the case, as here,
support gradI
--although if you embed gradI in a ring with minimal support, then it would hold
