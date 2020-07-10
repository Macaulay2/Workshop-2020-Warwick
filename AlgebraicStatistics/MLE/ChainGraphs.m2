restart;

load "GraphicalModelsMLE.m2"


MLEsolverCG = method();
MLEsolverCG(List, Digraph, Matrix):= (L,D,U) -> (
    vert:={}
    -- combine the vertices from different chain components
    -- but keep as separate list for the later use
    for l in L do vert=append(vert, vertexSet l)
    -- count the total number of vertices
    n:=# flatten vert;
    ---------------------------------------------------
    -- confirm input is compatible
    ---------------------------------------------------
    -- check that D is acyclic
    if isCyclic D then error("digraph must be acyclic");
    
    --!!!!!!!!!!!TBD!!!!!!!------------------------------
    -- add a test for checking that D contains only numbered variables
    -------------------------------------------------
    -- check that U has the correct dimension nxn
    if numgens target U =!= n then error("the rows/columns of the data matrix should be in 
	bijection to the variables in the chain; the variables should be unique ");
    -- check that the digraph represents the chain	
    if #L =!= #vertexSet D then error("the number of vertices in the digraph should 
	correspond to the number of chain components"); 
    ---!!!!!!!TBD!!!!!!!-------------------------------	
    -- possibly add a test for confirming uniqueness of variables (count as set and as a sequence
    -- for example)
    -------------------------------------------------
    --iterate over chain components	
    for i from 0 to #L-1 do
    (
    	-- check that each chain component is an undirected graph
	if class L#i =!= Graph then error("each chain component must be an undirected graph
	    (M2 object of type Graph)");
	pa_i:=parents(D,i);
	--merge the vertices of L#i and pa_i
	V_i:= join(take (vert,{i,pa_i}) )
	-- generate the corresponding undirected graph
	-- estimate K^*for that undirected graph
	-- estimate the S for pa_i
	-- check the difference between #Gamma and n
	-- assign K=to the formula
	-- iterate for all i
    );


);
