--Example for Bigraph
restart
loadPackage "StatGraphs"
B = bigraph{{1,5},{6}}
--returns a bigraph with bidirected edge between 1 and 5 and vertex 6 with no edges incidence to it. 


--Example for MixedGraph
--MixedGraph is a hashTable with keys Graph,Bigraph and Digraph
--vertices of a MixedGraph is union of all the vertices of Graph,Bigraph,Digraph
restart
loadPackage "StatGraphs"
G = graph{{1,2},{2,3},{3,4},{1,4},{1,5}}
D = digraph{{2,1},{3,1},{7,8}}
B = bigraph{{1,5}}
M = mixedGraph(G,D,B)
graph M
bigraph M
digraph M
vertices M


--Example for topologicalSort
--topologicalSort outputs a list of vertices in a topologically sorted order of a DAG.
--This method is used in topSort for creation of map.
restart
loadPackage "StatGraphs"
G = digraph{{2,1},{3,1}}
topologicalSort G


--Example for topSort 
--topSort outputs a HashTable with keys  digraph,map and newDigraph, where digraph is the original digraph,
--map is the relation between old ordering and the new ordering of vertices and 
--newDigraph is the Digraph with topologically sorted vertices.
restart
loadPackage "StatGraphs"
G = digraph{{2,1},{3,1}}
H = topSort G
H#digraph
H#map
H#newDigraph



--Example for collateVertices
--If M = mixedGraph(G,D,B) with vertices of G be V1, vertices of D be V2 and vertices of B be  V3.
--Then collateVertices(M) outputs a mixedGrpah with same edges as before but with V1 \cup V2 \cup V3 as the vertices of G,D and B.
restart
loadPackage "StatGraphs"
G = graph{{1,2},{2,3},{3,4},{1,4},{1,5}}
D = digraph{{2,1},{3,1},{7,8}}
B = bigraph{{1,5}}
M = mixedGraph(G,D,B)
collateVertices M


