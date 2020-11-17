--EXAMPLES WITH MULTIPLE EDGES THAT HAVE ZERO-DIMENSIONAL SCORE EQUATIONS IDEAL

--Graph with dir and bidir edges
--Type of multiple-edge: directed-bidirected
restart
loadPackage "GraphicalModelsMLE";
G = mixedGraph(digraph {{1,3},{2,4},{3,4}},bigraph {{3,4}});
S =  matrix {{7/20, 13/50, -3/50, -19/100}, {13/50, 73/100, -7/100, -9/100},{-3/50, -7/100, 2/5, 3/50}, {-19/100, -9/100, 3/50, 59/100}};
J = scoreEquations(gaussianRing G,S,SampleData=>false)
dim J,degree J
solverMLE(G,S,SampleData=>false)

--Graph with undir, dir and bidir edges
--Type of multiple-edge: directed-bidirected
restart
loadPackage "GraphicalModelsMLE";
G = mixedGraph(digraph {{1,3},{2,4},{3,4}},bigraph {{3,4}},graph{{1,2}});
R = gaussianRing G
S =  matrix {{7/20, 13/50, -3/50, -19/100}, {13/50, 73/100, -7/100, -9/100},{-3/50, -7/100, 2/5, 3/50}, {-19/100, -9/100, 3/50, 59/100}};
J = scoreEquations(gaussianRing G,S,SampleData=>false)
dim J,degree J
solverMLE(G,S,SampleData=>false)

--EXAMPLES WITH MULTIPLE EDGES THAT HAVE ONE-DIMENSIONAL SCORE EQUATIONS IDEAL

--Graph with undir, dir and bidir edges
--Type of multiple-edge: directed-undirected
restart
loadPackage "GraphicalModelsMLE";
G = mixedGraph(digraph {{1,3},{1,2},{2,4},{3,4}},bigraph {{3,4}},graph{{1,2}});
partitionLMG G
R = gaussianRing G
undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
S =  matrix {{7/20, 13/50, -3/50, -19/100}, {13/50, 73/100, -7/100, -9/100},{-3/50, -7/100, 2/5, 3/50}, {-19/100, -9/100, 3/50, 59/100}};
J=scoreEquations(R,S,SampleData=>false)
dim J, degree J  -- (1,4)


--Graph with undir and dir
--Type of multiple-edge: directed-undirected
restart
loadPackage "GraphicalModelsMLE";
G = mixedGraph(digraph {{1,3},{1,2},{2,4},{3,4}},graph{{1,2}});
partitionLMG G
R = gaussianRing G
undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
S =  matrix {{7/20, 13/50, -3/50, -19/100}, {13/50, 73/100, -7/100, -9/100},{-3/50, -7/100, 2/5, 3/50}, {-19/100, -9/100, 3/50, 59/100}};
J=scoreEquations(R,S,SampleData=>false)
dim J, degree J --(1,2)


--Graph with undir and dir edges
--Type of multiple-edge: directed-undirected
restart
loadPackage "GraphicalModelsMLE";
G = mixedGraph(digraph{{1,2}},graph{{1,2}});
partitionLMG G
R = gaussianRing G
undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
X=random(QQ^2,QQ^2)
S=X*transpose(X)
J=scoreEquations(R,S,SampleData=>false)
dim J, degree J --(1,2)

--Graph with dir and bidir edges
--Type of multiple-edge: directed-bidirected
restart
loadPackage "GraphicalModelsMLE";
G = mixedGraph(digraph{{1,2}},bigraph{{1,2}});
partitionLMG G
R = gaussianRing G
undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
X=random(QQ^2,QQ^2)
S=X*transpose(X)
J=scoreEquations(R,S,SampleData=>false)
dim J, degree J --(1,2)

