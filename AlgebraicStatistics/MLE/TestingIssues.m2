------------------
---#Issue 117----- 
------------------

restart
debug needsPackage "GraphicalModels"
R=markovRing(2,3,2)
R.markovVariables
R.markovVariables#(1,1,1)
T=ring(R.markovVariables#(1,1,1))
describe T
keys T
T.markovVariables

v=1
A=R
v=v-1
p:= i -> A.markovVariables#i
d=(2,3,2)
e=drop(d, {v,v})
--S = markovRing (e)
S=markovRing(e,Coefficients=>coefficientRing(A))
dv= d#v
F:= toList apply(((#e):1) .. e, i -> (
	       sum(apply(toList(1..dv), j -> (
			      newi := join(take(i,v), {j}, take(i,v-#d+1));
			      p newi)))))
map(A,S,F)

F
F_0
promote(F_0,R)

------Output from branch AlgebraicStatisticsRoser
i14 : F_0

o14 = p      + p
       1,1,1    2,1,1

o14 : R

------Output from branch AlgebraicStatistics

i14 : F_0

o14 = p      + p
       1,1,1    2,1,1

o14 : QQ[p     , p     , p     , p     , p     , p     , p     , p     , p     , p     , p     , p     ]
          1,1,1   1,1,2   1,2,1   1,2,2   1,3,1   1,3,2   2,1,1   2,1,2   2,2,1   2,2,2   2,3,1   2,3,2


-- The fact that the elements in F are not recognised as elements in R makes promote(F_i,R) to fail

A.markovVariables
p(1,1,1)
A.markovVariables#(1,1,1)

restart
needsPackage "GraphicalModels"
F = hiddenMap(1,markovRing(2,3,2));

------------------
---#Issue 116----- 
------------------

restart
needsPackage "GraphicalModels"
d=(2,3,4,5);
R = markovRing d
marginMap(1,R)
inverseMarginMap(1,R)

------------------
---#Issue 112----- 
------------------
restart
needsPackage "GraphicalModels"
G = digraph {{1,{}}, {2,{}}}
R = markovRing (2,2)
discreteVanishingIdeal (R,G)

G=graph{{1,2},{1,3},{2,3}}
discreteVanishingIdeal (R,G)
it doesn't say that the graph is incompatible. Instead, it gives
stdio:8:1:(3): error: Number of vertices of graph does not match size of ring
class G
D=digraph{{1,2}}
class G===Digraph
class D===Digraph
discreteVanishingIdeal(R,D)

------------------
---#Issue 121----- 
------------------
restart
needsPackage "GraphicalModels"
G=digraph{{1,2},{1,3},{2,3}}
R=gaussianRing G

-- Example from gaussianParametrization documentation

G = mixedGraph(digraph {{b,{c,d}},{c,{d}}},bigraph {{a,d}})
R = gaussianRing G
M = gaussianParametrization(R,SimpleTreks=>true)

-- Elina's example
restart
needsPackage "GraphicalModelsMLE"
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R = gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
scoreEquations(R,U)
solverMLE(G,U)	    


----------------------------------------------
----------------------------------------------
----------------------------------------------
restart
installPackage "GraphicalModels"

