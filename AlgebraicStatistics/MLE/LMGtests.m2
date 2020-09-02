restart
installPackage "Graphs"
installPackage "StatGraphs"
installPackage "GraphicalModels"

loadPackage "Graphs"
loadPackage "StatGraphs"
loadPackage "GraphicalModels"

G=graph{{1,2},{1,3},{2,3}}
D=digraph{{1,6},{4,7}}
B=bigraph{{5,6},{6,7}}

W:=vertexSet B;
U:=collateVertices g - W;


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

-- random function with options
r = method(Options => {Test => true});
r RR := opts -> x -> (
    
    if opts.Test then x^2 else x+1
    );
r(0.5)
r(0.5, Test => true)
r(0.5, Test => false)

r = method(Options => {Test => true});
r RR := opts -> x -> (
    y:=x;
    if opts.Test then y=x^2;
    return y
    );
r(0.5)
r(0.5, Test => false)


