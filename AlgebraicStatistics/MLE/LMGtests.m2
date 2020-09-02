restart
installPackage "Graphs"
installPackage "StatGraphs"
installPackage "GraphicalModels"

G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

-- MixedGraph with all components
g=mixedGraph(G,D,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without B
g=mixedGraph(G,D)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without D
g=mixedGraph(G,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without G
g=mixedGraph(D,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without D,B
g=mixedGraph(G)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without G,B
g=mixedGraph(D)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without G,D
g=mixedGraph(B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R