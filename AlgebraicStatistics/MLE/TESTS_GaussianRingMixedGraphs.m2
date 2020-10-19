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

-----------------
--TEST EXAMPLES
-----------------

restart
needsPackage "StatGraphs"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)

vertices G
vertices B
vertices D

vertexSet G
vertexSet B
vertexSet D

for ed in edges D do (if member(ed_0,vertices B) and member(ed_1,vertices G) then print("Directed edges cannot go from vertices adjacent to a bidirected edge to vertices adjacent to an undirected edge"))


restart 
needsPackage "StatGraphs"
G=graph{{1,2}}
D=digraph{{2,3},{4,5}}
B=bigraph{{3,4},{5,6}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)






restart 
needsPackage "StatGraphs"
G=graph{{1,2}}
D=digraph{{1,3},{4,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)



restart
needsPackage "GraphicalModels"
G=graph{{1,2},{9,10}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4},{6,10}}

g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)

restart
needsPackage "GraphicalModels"
G=graph{{1,2},{9,10}}
D=digraph{{1,3},{2,1},{6,7},{7,8},{6,8}}
B=bigraph{{5,4},{6,10}}

g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)


restart
needsPackage "GraphicalModels"
G=graph{}
D=digraph{{1,3},{2,1}}
B=bigraph{{3,1}}

g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)


restart
needsPackage "GraphicalModels"
G=graph{{1,1}}
D=digraph{{1,3},{2,1}}
B=bigraph{{4,1}}
g=mixedGraph(G,D,B)
isMixedGraphSimple(g)
partitionLMG(g)
