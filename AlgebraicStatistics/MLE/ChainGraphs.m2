restart;

load "GraphicalModelsMLE.m2"

--***************************--
--INTERNAL ROUTINES
--***************************--

------------------------------------------------------
-- Confirm that digraph satisfies the conditions of CG 
-- for digraphs
------------------------------------------------------


confirmDir=(D, vertD, vertG,  compG, compCountG,adjD)->(
    -- check that D is acyclic
    if isCyclic D then error("digraph must be acyclic");
    
        
    -- check that the vertices of the digraph are the same as in the graph	
    if not isSubset(vertD, vertG) then error("the vertices of the digraph should be the same 
	as in the graph"); 
	
    -- check that D does not contain any edges between vertices from the same component	
    block:={};
    ind:={};
    for i from 0 to compCountG - 1 do
    (
	-- check which vertices of D are in the i-th chain component
	block_i:=set compG#i * set  vertD;
	-- extract their indices in the adjacency matrix adjD
	ind_i:= apply(block_i,j->position(vertD, k->k===j));
	-- check if the corresponding submatrix is a zero matrix
	if submatrix(adjD,ind_i,ind_i)=!=0 then error ("all edges 
	    inside a chain component should be undirected");
	);
    if #block =!= compCountG then error ("The list of chain components L and the digraph D are incompatible");
    if #block =!= # ind then error ("Internal error: mistake in computation of index set");
    );

-------------------------------------------------------
-- Confirm that there are no directed cycles
-------------------------------------------------------
confirmNoDirCycles=(adjD, compCountG)->(
    -- form the adjacency matrix of the graph of chain components 
       adjCCG:= mutableMatrix{{compCountG:0},{ compCountG:0}};
	   for i from 0 to  compCountG - 1 do
	   (
	       for j from 0 to  compCountG - 1 do
	       (if submatrix(adjD,ind_i,ind_j)=!=0 then adjCCG_(i,j)=1 else adjCCG_(i,j)=0
	       );
     -- check that the graph of chain components is acyclic
       if isCyclic (digraph adjCCG) then error("the chain graph should not contain directed cycles");    
    );

--**************************--
--  METHODS 	      	   	  --
--**************************--

MLEsolverCG= method();
MLEsolverCG(Graph, Digraph, Matrix):=(G,D,U)->(
    -------------------------------------------------------------------------
    -- Pre-process input
    -------------------------------------------------------------------------    
    -- confirm input has the correct type
     if class G =!= Graph then error("the first argument  must be an undirected graph
	    (M2 object of type Graph)");
    
     if class D =!= Digraph then error("the second argument  must be a directed graph
	    (M2 object of type Digraph)");
	    
     if class U =!= Matrix then error("the third  argument  must be a directed matrix");		    
    -- extract information about the graph G
    vertG:=vertexSet G;
    varCount:=# vertG;
    edgG:=edges G;
    compG:= connectedComponents G;
    compCountG:= # compG;
    
    -- extract information about the digraph D
    vertD:=vertexSet D;
    underlyingEdgesD:==edges (underlyingGraph D);
    adjD := adjacencyMatrix D;

    -- confirm D satisfies the conditions
    confirmDir(D, vertD, vertG,  compG, compCountG,adjD);
    
    -- confirm that there are no directed cycles
    confirmNoDirCycles=(adjD, compCountG);
    
    -- confirm U is compatible with the chain graph
    if numgens target U =!= varCount then error("the rows/columns of the data matrix should be in 
	bijection to the variables in the chain ");
    
    -------------------------------------------------------------------------
    -- Compute the MLE using (5.52) from Lauritzen, "Graphical Models", 1996 
    -------------------------------------------------------------------------
    -- initialize all-zero K and M of correct dimension
    K:= mutableMatrix{{varCount:0},{varCount:0}};
    -- combine the directed and undirected edges
    combCG:=addEdges'(G,underlyingEdgesD);
    
    -- iterate over chain components    
    for C in compG do (
	-- determine the pa(c)
	paC:=set {};
	for a in c do paC=paC+parents(D,a);
    	-- generate the undirected graph c_union_paC
	-- 1. Get the vertex set
	vertC:=c + paC;
	-- 2. Add undirected edges between the parents
	adjC:= mutableMatrix{{#vertC:0},{#vertC:0}};   	   
    	for i in paC do(
	       for j in paC  do adjCCG_(i,j)=1;
	       );
	-- 3. Add any edge that was originally in the chain graph
	
	-- Estimate K_C for the undirected graph  c_union_paC
	
	-- Take (K_C)^Gamma
	KCG:=
	-- Compute ssd_paC
	
	-- Take  (ssd_paC^{-1})^Gamma
	SSDInvG:=
	
	-- Compute M (define n!)
	M:=KCG - n*SSDInvG;
	-- Compute K
	K=K+M;
	);
    return K;
    );


end


MLEsolverCG = method();
MLEsolverCG(List, Digraph, Matrix):= (L,D,U) -> (
    -- Add Caveat to documentation: vertices that are incident to only directed edges 
    -- should be entered as Graphs with 1 vertex.
    
    -- add a check for connectivity components?
    ------------------------------------------------------
    -- construct the vertex set of the CG
    ------------------------------------------------------
    -- combine the vertices but remember the chain component
    vert:={}
    for l in L do vert=append(vert, vertexSet l);

    -- count the total number of vertices
    n:=# flatten vert;
    
    --all undirected edges 
    edgUndir:={};
    for l in L do edgUndir=append(edgUndir, edges  l);
    
    -- all directed edges
    edgDir:= edges D;
    ------------------------------------------------------
    -- confirm input is compatible
    ------------------------------------------------------
    -- all elements of the list are undirected graphs
    for l in L do 
       (if class l =!= Graph then error("each chain component must be an undirected graph
	    (M2 object of type Graph)"));
    
    -- all vertices are unique
    if # unique(flatten vert)=!=n then error("vertices should
	be unique");
	
    
    -- check that D is acyclic
    if isCyclic D then error("digraph must be acyclic");
    
 
    -- check that U has the correct dimension nxn
    if numgens target U =!= n then error("the rows/columns of the data matrix should be in 
	bijection to the variables in the chain; the variables should be unique ");
	
    -- check that the vertices of the digraph are the same as in the chain components	
    if not isSubset(vertexSet D, flatten vert) then error("every vertex of the digraph should appear
	in exactly one chain component"); 
	
    -- check that D does not contain any edges between vertices from the same component	
    adjD := adjacencyMatrix D;
    block:={};
    ind:={};
    for i from 0 to #L-1 do
    (
	-- check which vertices of D are in the i-th chain component
	block_i:=set vert#i * set  vertexSet D;
	-- extract their indices in the adjacency matrix adjD
	ind_i:= apply(block_i,j->position(vertexSet D, k->k===j));
	-- check if the corresponding submatrix is a zero matrix
	if submatrix(adjD,ind_i,ind_i)=!=0 then error ("all edges 
	    inside a chain component should be undirected");
	);
    if #block =!= #L then error ("The list of chain components L and the digraph D are incompatible");
    if #block =!= # ind then error ("Internal error: mistake in computation of index set");
    -- check that all edges between components are directed
    
    -- (this is done earlier in the check for uniqueness of vertices.
    -- because if there is an undirected edge between two chain
    -- components then it will be given in a third chain component,
    -- so the vertices are not unique)
    -----------------------------------------------------------------
    -- form a graph of chain components
    -----------------------------------------------------------------
    -- form the adjacency matrix of the graph of chain components 
       adjCCG:= mutableMatrix{{#block:0.9},{#block:0.9}};
	   for i from 0 to # block - 1 do
	   (
	       for j from 0 to # block - 1 do
	       (if submatrix(adjD,ind_i,ind_j)=!=0 then adjCCG_(i,j)=1 else adjCCG_(i,j)=0
	       );
       if ring adjCCG=!= ZZ then error ("Internal error: mistake in adjCCG computation"); 
       
    -- isCyclic CCG?      
	        	    
    ----------------------------------------------------
    -- Compute the MLE using (5.52) from Lauritzen, "Graphical Models", 1996 
    ----------------------------------------------------
    --iterate over chain components	
    for i from 0 to #L-1 do
    (
    	-- determine the parents of a chain component
	pa_i:= positions(flatten entries (adjacencyMatrix G2)_{i}, j -> j != 0);
	--merge the vertices of L#i and pa_i
	V_i:= flatten (take (vert,{i,pa_i}));
	-- generate the corresponding undirected graph
	
	-- vertices in pa_i
	paV_i:= flatten (take (vert,pa_i));
	-- edges between vertices in pa_i
	G_i:=graph(paV_i, subsets(paV_i, 2), EntryMode => "edges");
	
	--undirected edges between a pair $(\alpha,\beta)$ if there is an edge, directed or undirected, between them in 
	--the chain graph $G$.
	edgUndir
	 edgDir
	-- estimate K^*for that undirected graph
	-- estimate the S for pa_i
	-- check the difference between #Gamma and n
	-- assign K=to the formula
	-- iterate for all i
    );

);

parents = method()
parents (Digraph, Thing) := Set => (G, v) -> (
    i := position(vertexSet G, u -> u === v);
    if i === null then error "v is not a vertex of G.";
    set (vertexSet G)_(positions(flatten entries (adjacencyMatrix G)_{i}, j -> j != 0))
    )

children = method()
children (Digraph, Thing) := Set => (G, v) -> (
    i := position(vertexSet G, u -> u === v);
    if i === null then error "v is not a vertex of G.";
    set (vertexSet G)_(positions(first entries (adjacencyMatrix G)^{i}, j -> j != 0))
    )