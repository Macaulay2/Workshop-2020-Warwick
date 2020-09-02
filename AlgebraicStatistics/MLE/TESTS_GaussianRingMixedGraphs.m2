restart
loadPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

g=mixedGraph(G,D,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

keys R
R.mixedGraph
keys R.mixedGraph
R.mixedGraph#graph 
R.mixedGraph#graph#Graph
R.mixedGraph#graph#Digraph
R.mixedGraph#graph#Bigraph
keys R.gaussianRingData --key not found in hash table
R.gaussianRingData --key not found in hash table

restart
loadPackage "GraphicalModels"
G=graph{{1,2},{9,10}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

g=mixedGraph(G,D,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R
