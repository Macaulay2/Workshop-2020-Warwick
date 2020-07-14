restart;

load "GraphicalModelsMLE.m2"


MLEsolverCG = method();
MLEsolverCG(List, Digraph, Matrix):= (L,D,U) -> (
    ------------------------------------------------------
    -- preprocess input
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
	
    -- check that the digraph represents the chain	
    if #L =!= #vertexSet D then error("the number of vertices in the digraph should 
	correspond to the number of chain components"); 
   	
    ----------------------------------------------------
    -- Compute the MLE using (5.52) from Lauritzen, "Graphical Models", 1996 
    ----------------------------------------------------
    -- sort the vertices in D
    orderedComp = sort vertexSet D;
    --iterate over chain components	
    for i from 0 to #L-1 do
    (
    	-- determine the parents of a chain component
	pa_i:=parents(D,orderedComp#i);
	--merge the vertices of L#i and pa_i
	--!! needs to be updated to use orderedComp
	V_i:= join(take (vert,{i,pa_i}))
	-- generate the corresponding undirected graph
	-- estimate K^*for that undirected graph
	-- estimate the S for pa_i
	-- check the difference between #Gamma and n
	-- assign K=to the formula
	-- iterate for all i
    );
    -----------------------------------------------
    -- Add a self-check function for users to confirm the correspondence between element of L and
    -- vertices of D
    -----------------------------------------------

);
