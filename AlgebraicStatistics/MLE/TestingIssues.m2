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

------------------
---#Issue 115----- 
------------------
restart
needsPackage "GraphicalModelsMLE"

G= mixedGraph(digraph {{b,{c,d}},{c,{d}}},bigraph {{a,d}})
R=gaussianRing G
T=trekIdeal(R,G)
--stdio:4:3:(3): error: expected argument 2 to be a visible list
R.gaussianRingData
R.?gaussianRingData
R.?graphType
R.graphType
Stmts=trekSeparation G
vv=sort vertices G
SM=covarianceMatrix R
Stmts=={} --The issue was () instead of {}

R.graph#graph#Digraph
R.graph#graph#Digraph===digraph{} 
R.graph#graph#Bigraph===bigraph{}

restart
needsPackage "GraphicalModelsMLE"

G= digraph {{b,{c,d}},{c,{d}}}
R=gaussianRing G
T=trekIdeal(R,G)

restart
needsPackage "GraphicalModelsMLE"

G= graph {{b,{c,d}},{c,{d}}}
R=gaussianRing G
T=trekIdeal(R,G)

g=mixedGraph(G)
R2=gaussianRing g
T=trekIdeal(R2,g)


----------------------------------------------
--------------Issue #135----------------------
----------------------------------------------


restart
installPackage "GraphicalModels"

restart
needsPackage "GraphicalModels"

--examples for localMarkov
D = digraph {{1,{2,3}}, {2,{4,5}}, {3,{5,6}}, {4,{7}}, {5,{7}},{6,{7}},{7,{}}}
netList pack (3, localMarkov D) 
      
--examples for markovRing      
d=(2,3,4,5);
R = markovRing d;
numgens R
R_0, R_1, R_119 --here are some of the variables in the ring

restart
installPackage "GraphicalModelsMLE"

restart
installPackage "StatGraphs"
installPackage "GraphicalModels"
--/home/roserhp/.Macaulay2/local/share/Macaulay2/StatGraphs.m2:49:17:(2):[26]: 
--error: assignment to protected built-in global variable 'topologicalSort'

restart
uninstallPackage "Graphs"
uninstallPackage "StatGraphs"
uninstallPackage "GraphicalModels"
uninstallPackage "GraphicalModelsMLE"
uninstallPackage "EigenSolver"

restart
installPackage "Graphs"
check Graphs

installPackage "StatGraphs"
check StatGraphs
-- No existing tests

installPackage "GraphicalModels"
check GraphicalModels

-- running test 9 of package GraphicalModels on line 2990 in file ./GraphicalModels.m2
--    rerun with: check_9 "GraphicalModels"
--making test results
 ulimit -c unlimited; ulimit -t 700; ulimit -m 850000; ulimit -s 8192; ulimit -n 512;  cd /tmp/M2-4957-0/192-rundir/; GC_MAXIMUM_HEAP_SIZE=400M "/usr/bin/M2-binary" --int --no-randomize --no-readline --silent --stop --print-width 77 -e 'needsPackage("GraphicalModels", Reload => true, FileName => "/home/roserhp/Workshop-2020-Warwick/AlgebraicStatistics/MLE/GraphicalModels.m2")' <"/tmp/M2-4957-0/191.m2" >>"/tmp/M2-4957-0/191.tmp" 2>&1
/tmp/M2-4957-0/191.tmp:0:1: (output file) error: Macaulay2 exited with status code 1
/tmp/M2-4957-0/191.m2:0:1: (input file)
M2: *** Error 1

restart
needsPackage "GraphicalModels"
G = digraph {{a,{b,c}}, {b,{c,d}}, {c,{}}, {d,{}}}
R = gaussianRing G
assert(sort gens R === sort {s_(a,a), s_(a,b), s_(a,c), s_(a,d), s_(b,b), s_(b,c), s_(b,d), s_(c,c), s_(c,d), s_(d,d)})
--This fails because gaussianRing digraph has more variables now!
toExternalString gens R
assert(sort gens R === sort {l_(a,c),l_(a,b),l_(b,c),l_(b,d),p_(a,a),p_(b,b),p_(c,c),p_(d,d),s_(a,a),s_(a,b),
     s_(a,c),s_(a,d),s_(b,b),s_(b,c),s_(b,d),s_(c,c),s_(c,d),s_(d,d)})
gens R
directedEdgesMatrix R
bidirectedEdgesMatrix R

installPackage "EigenSolver"
installPackage "GraphicalModelsMLE"
check GraphicalModelsMLE

--TEST 6///
G=graph{{1,2},{2,3},{3,4},{1,4}}
R=gaussianRing(G)
U=random(ZZ^4,ZZ^4)
J=scoreEquations(R,U)
assert(dim J===0)
assert(degree J===5)
///   

--TEST 8///
L={matrix{{1,0},{0,1}},matrix{{-2,0},{0,1}},matrix{{sqrt(-1),0},{0,sqrt (-1)}},matrix{{0.0001*sqrt(-1),0},{0,0.0000001*sqrt (-1)}},matrix{{0,0},{0,0}}}
Y = checkPD(L)
Y2=checkPSD(L)
B = {matrix{{1, 0}, {0, 1}}}
B2 = {matrix{{1, 0}, {0, 1}},matrix{{0,0},{0,0}}}
assert(Y===B)	
assert(Y2==B2)
///





restart
needsPackage "GraphicalModels"
g = mixedGraph(digraph {{a,{b}},{b,{c}}},bigraph {{a,c}, {b,c}})
R = gaussianRing g
H = identifyParameters R
       
g#graph#Graph       
g#graph#Bigraph
G=graph collateVertices g
g#graph#Graph===graph{}


----------------------------------------------
--------------Issue #135----------------------
----------------------------------------------
restart
uninstallPackage "Graphs"
uninstallPackage "StatGraphs"
uninstallPackage "GraphicalModels"
uninstallPackage "GraphicalModelsMLE"
uninstallPackage "EigenSolver"

installPackage "Graphs"
check Graphs

installPackage "ReactionNetworks"
check ReactionNetworks
N = reactionNetwork "A <--> B"
S = stoichiometricMatrix N
E = stoichiometricConeKer N

N = oneSiteModificationA()
S = stoichiometricMatrix N
E = stoichiometricConeKer N

/usr/libexec/Macaulay2/bin/4ti2int64: error while loading shared libraries: libglpk.so.40: cannot open shared object file: No such file or directory
stdio:14:5:(3): error: error occurred while executing external program 4ti2: rays

 ulimit -c unlimited; ulimit -t 700; ulimit -m 850000; ulimit -s 8192; ulimit -n 512;  cd /tmp/M2-102-0/186-rundir/; GC_MAXIMUM_HEAP_SIZE=400M "/usr/bin/M2-binary" --int --no-randomize --no-readline --silent --stop --print-width 77 -e 'needsPackage("ReactionNetworks", Reload => true, FileName => "/usr/share/Macaulay2/ReactionNetworks.m2")' <"/tmp/M2-102-0/185.m2" >>"/tmp/M2-102-0/185.tmp" 2>&1
/tmp/M2-102-0/185.tmp:0:1: (output file) error: Macaulay2 exited with status code 1
stdio:3:9:(3):[1]: error: error occurred while executing external program 4ti2: rays
/tmp/M2-102-0/185.m2:0:1: (input file)
M2: *** Error 1

--the only error comes from an external program

installPackage "Visualize"
check Visualize

--examples ok
--no tests

installPackage "BinomialEdgeIdeals"


G={{1,2},{2,3},{3,1}}
d = disconnectors(G)
d = disconnectors(G,EffectiveOnly=>true)
stdio:21:5:(3): error: error occurred while executing external program 4ti2: markov

S={1}
isEffective(G,S)
isDisconnector(G,S)
disconnectors(G,EffectiveOnly=>true)
--errors come from an external program

check BinomialEdgeIdeals
--tests ok

installPackage "Chordal"
check Chordal
--ok

installPackage "SimplicialPosets"
check SimplicialPosets
--ok

installPackage "PhylogeneticTrees"
T = leafTree(4, {{0,1}})
phyloToric42(T, CFNmodel)
stdio:36:1:(3): error: error occurred while executing external program 4ti2: markov

check PhylogeneticTrees
--many errors, but all seem to come from external program 4ti2

installPackage "Matroids"
check Matroids
--ok

installPackage "NautyGraphs"
check NautyGraphs
--ok

installPackage "Posets" 
hibiRing booleanLattice 2
hibiRing chain 4
hibiRing(divisorPoset 6, Strategy => "4ti2")

P = poset {{1,2}, {2,4}, {3,4}, {3,5}}
pPartitionRing(divisorPoset 6, Strategy => "4ti2")

check Posets

needsPackage "FourTiTwo"

----------------------------
--- #Issue 155
---------------------------
restart
needsPackage "StatGraphs"
U = graph{{1,2},{2,3},{1,3}}
D = digraph{{1,4},{3,7},{8,9}}
B = bigraph{{4,5},{5,6},{7,9}}
G = mixedGraph(U,D,B)
partitionLMG G

