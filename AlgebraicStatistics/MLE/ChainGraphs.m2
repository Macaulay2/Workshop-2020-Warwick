restart;

load "GraphicalModelsMLE.m2"


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
    A := adjacencyMatrix D;
    for i from 0 to #L-1 do
    (
	-- check which vertices of D are in the i-th chain component
	block_i:=set vert#i * set  vertexSet D;
	-- extract their indices in the adjacency matrix A
	ind_i:= apply(block_i,j->position(vertexSet D, k->k===j));
	-- check if the corresponding submatrix is a zero matrix
	if submatrix(A,ind_i,ind_i)=!=0 then error ("all edges 
	    inside a chain component should be undirected");
	);
    -- check that all edges between components are directed
    -- (this is done earlier in the check for uniqueness of vertices.
    -- because if there is an undirected edge between two chain
    -- components then it will be given in a third chain component,
    -- so the vertices are not unique)
    
    -- form a graph of chain components	
    ----------------------------------------------------
    -- Compute the MLE using (5.52) from Lauritzen, "Graphical Models", 1996 
    ----------------------------------------------------
    -- sort the vertices in D
    orderedComp := sort vertexSet D;
    --iterate over chain components	
    for i from 0 to #L-1 do
    (
    	-- determine the parents of a chain component
	pa_i:=parents(D,orderedComp#i);
	-- identify the corresponding chain components
	orderedpa_i:=apply(pa_i,j->position(orderedCom,k->k===j));
	--merge the vertices of L#i and pa_i
	V_i:= join(take (vert,{i,orderedpa_i}))
	-- generate the corresponding undirected graph
	-- estimate K^*for that undirected graph
	-- estimate the S for pa_i
	-- check the difference between #Gamma and n
	-- assign K=to the formula
	-- iterate for all i
    );

);
