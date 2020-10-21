restart
needsPackage "StatGraphs"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

g=mixedGraph(G,D,B)
partitionLMG(g)

edG=edges G
edD=edges D
edB=edges B

EG={}
for i in edG do EG=join(EG,{toList i})
EB={}
for i in edB do EB=join(EB,{toList i})
EG
EB

G2=graph{join(EG,edD,EB)}
isSimple G2


L=join(edG,edD,edB)
L2=unique(L)
#L==L2

setG=set vertices G
setB=set vertices B
setint=setB*setG 
emptyset=set {}
setint===emptyset

listG= toList setG
g=mixedGraph(G,D,B)
v=vertices g


restart
debug needsPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}
g=mixedGraph(G,D,B)
R=gaussianRing g
undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

--test for gaussianMatrices
Stmts = {{{1,2},{3},{4}}, {{1},{3},{}}}
gaussianMatrices(R,Stmts)

--test for conditionalIndependenceIdeal
Stmts = {{{1,2},{3},{4}}, {{1},{3},{}}}
conditionalIndependenceIdeal(R,Stmts)

--test for gaussianParametrization
gaussianParametrization R
trekIdeal(R,g)

--tests for gaussianRing ZZ
S=gaussianRing 4
Stmts = {{{1,2},{3},{4}}, {{1},{3},{}}}
conditionalIndependenceIdeal(R,Stmts)
gaussianMatrices(S,Stmts)
covarianceMatrix S


keys R
R.graphType
keys R.graphType
R.data
keys R.data
R.data#sVar
R.data#compU
R.data#nn

R.graph
keys R.graph

R.gaussianVariables
keys R.gaussianVariables
R.gaussianVariables#(k_(2,2)) --doesn't work
R.gaussianVariables#k22
R.gaussianVariables.k_(2,2)

pos := (h, x) -> position(h, i->i===x)



R.mixedGraph
keys R.mixedGraph
R.mixedGraph#graph 
R.mixedGraph#graph#Graph
R.mixedGraph#graph#Digraph
R.mixedGraph#graph#Bigraph
keys R.gaussianRingData --key not found in hash table
R.gaussianRingData --key not found in hash table





partitionLMG MixedGraph := g -> (
   --check it's a simple graph
   --if isMixedGraphSimple(g)==false then print ("The input should be a simple mixedGraph.");
   --retrive graph, bigraph and digraph
   G= g#graph#Graph
   B= g#graph#Bigraph
   D= g#graph#Digraph
   --check D is a DAG: IS THIS ENOUGH???
   --if isCyclic D==false then error ("The directed part of the mixedGraph contains cycles.");
   --naive partition (vertices only adjacent to directed edges are not considered) 
   U=vertices G
   W=vertices B
   --check there are no common vertices for undirected and bidirected edges
   --if not (set U * set W===set {}) then 
   --error("Vertices cannot be adjacent to bidirected and undirected edges simultaneously.");
   --check there are no directed edges from set U to set W
   --for e in edges D do (if (member(e_0,set U) and member(e_1,set W)) 
   --then error("Directed edges cannot go from vertices adjacent to a bidirected edge to vertices adjacent to an undirected edge"));
   --remaining vertices (only adjacent to directed edges)
   V=set vertices g-set W-set U
   --if V===set{} then return (U,W);
   --------------------------------------------------
   --------------------------------------------------
   for v in toList V do (
      asc=forefathers(D,v);
      if asc*set W===set{} 
      then U=append(U,v) 
      else (if asc*set U==set{} then W=append(W,v) else return "Sorry");
   );
  

   asc=forefathers(D,5)
   asc*set W


--   );
  -- U,W
   -------------------------------------------------
   )

-----------------
--TEST EXAMPLES
-----------------

--EXAMPLE 1: LMG (check if 1,2,3 is a directed cycle or not)
restart
needsPackage "StatGraphs"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)

--({1, 2, 8, 3, 6, 7}, {5, 4}) OK but not ordered. Is it necessary?

--EXAMPLE 2: LMG with all vertices adjacent to undir o bidir
restart 
needsPackage "StatGraphs"
G=graph{{1,2}}
D=digraph{{2,3},{4,5}}
B=bigraph{{3,4},{5,6}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)

--({1, 2}, {3, 4, 5, 6}) OK



--EXAMPLE 3: Directed edges cannot go from vertices adjacent to a bidirected edge to vertices adjacent to an undirected edges
restart 
needsPackage "StatGraphs"
G=graph{{1,2}}
D=digraph{{1,3},{4,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)


--EXAMPLE 4:Vertices cannot be adjacent to bidirected and undirected edges simultaneously.
restart
needsPackage "GraphicalModels"
G=graph{{1,2},{9,10}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4},{6,10}}

g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)

--EXAMPLE 5: LMG
restart
needsPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4},{6,10}}

g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)
--({1, 2, 3}, {5, 4, 6, 10, 8, 7}) OK, not ordered

--EXAMPLE 6: not simple
restart
needsPackage "GraphicalModels"
G=graph{{1,1}}
D=digraph{{1,3},{2,1}}
B=bigraph{{4,1}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)

--EXAMPLE 7: Non empty intersection of anscestors and both set U and set W
-- not well-ordered (actually there exists another partition that makes it right without changing vertices: {1,2},{3,4,5,6,7})
restart
needsPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{2,4},{2,5},{3,5},{3,6}}
B=bigraph{{3,7}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
(U,W)=partitionLMG(g)
sort U
sort W
join(sort U,sort W)
R=gaussianRing g
K=undirectedEdgesMatrix R
L=directedEdgesMatrix R 
P=bidirectedEdgesMatrix R
covarianceMatrix R

--EXAMPLE 8: Non empty intersection of anscestors and both set U and set W
--well ordered
restart
debug needsPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{2,3},{2,5},{4,5},{4,6}}
B=bigraph{{4,7}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
(U,W)=partitionLMG(g)
sort U
sort W
join(sort U,sort W)
R=gaussianRing g
K=undirectedEdgesMatrix R
L=directedEdgesMatrix R 
P=bidirectedEdgesMatrix R
directSum(K,P)
covarianceMatrix R
R.gaussianVariables

    
--place remaining vertices in either U or W depending on their ancestors
   -- Case 1: all remaining vertices have higher number than vertices in U
   if max U < min toList V then return (U,join(W,toList V));
   -- Case 2: 
   for v in toList V do (
      asc:=forefathers(D,v);
      if asc*set W===set{} 
      then U=append(U,v) 
      else W=append(W,v);
   );
         

-----------------------------------
-----------------------------------
-----------------------------------
-- Olga's tests in issue #95
restart
needsPackage "GraphicalModels"
B=bigraph{{5,6},{6,7}} 
mixedGraph B
R=gaussianRing mixedGraph B
describe R
R.graph#graph#Graph
graph{}
R.graph#graph#Graph===graph{}
R.graph#graph#Bigraph===bigraph{}
R.graph#graph#Digraph===digraph{}
R2=gaussianRing B
describe R2
R2.gaussianRingData

G=graph{{1,2},{1,3},{2,3}}
gaussianRing mixedGraph G
gaussianRing G




-- Olga's tests in issue #100
restart
debug needsPackage "GraphicalModels"
D=digraph{{1,6},{4,7}}
R4=gaussianRing D
R4.gaussianRingData
R4.gaussianVariables
keys R4.gaussianVariables
covarianceMatrix R4
undirectedEdgesMatrix R4
directedEdgesMatrix R4
bidirectedEdgesMatrix R4
time gaussianVanishingIdeal R4
time gaussianVanishingIdeal(R4, oldVersion=>true)
conditionalIndependenceIdeal R4
Stmts = {{{1,6},{7},{4}}}
gaussianMatrices(R4,Stmts)
conditionalIndependenceIdeal(R4,Stmts)

debug needsPackage "GraphicalModelsMLE"
verticesInRing R4

R=R4

G = R.graph#graph#Digraph
vv = sort vertices G
n = #vv
v = (topSort G)#map
v = hashTable apply(keys v, i->v#i=>i)
v = apply(n,i->v#(i+1))
P = toList apply(v, i -> toList parents(G,i))
s = R.gaussianRingData#sVar
L = select(gens R, v -> first baseName v==s)
nx =#L -- Careful! old gaussianRing Digraph only had s variables
ny = max(P/(p->#p))
x := local x
y := local y
S = (coefficientRing R)[x_0 .. x_(nx-1),y_0 .. y_(ny-1)]
newvars = apply(ny, i -> y_i)
H = hashTable apply(nx,i->L#i=>x_i)
pos = (h, xx) -> position(h, i->i===xx)
sp = (i,j) -> if pos(vv,i) > pos(vv,j) then H#(s_(j,i)) else H#(s_(i,j));
I = trim ideal(0_S)
for i from 1 to n-1 do (
    J = ideal apply(i, j -> sp(v#j,v#i) - sum apply(#P#i, k ->y_k * sp(v#j,P#i#k)));
    I = eliminate(newvars, I + J);
)
F = map(R,S,apply(nx,i->x_i=>R.gaussianVariables#(L_i))|apply(newvars,i->i=>0))
F(I)

sp(v#1,v#1)     
s
s_(1,1)
L
H
H#(s_(1,1))
L_0
H#(L_0)
R.gaussianVariables#(L_0)

gens R
gens S
apply(nx,i->x_i=>L_i)
apply(newvars,i->i=>0)



R=gaussianRing mixedGraph D
R.gaussianVariables
keys R.gaussianVariables
covarianceMatrix R4


restart
needsPackage "GraphicalModels"
D=digraph{{1,6},{4,7}}
R4=gaussianRing D
time gaussianVanishingIdeal R4


R=gaussianRing mixedGraph D
time gaussianVanishingIdeal R

