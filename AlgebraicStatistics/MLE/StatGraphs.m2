newPackage(
        "StatGraphs",
        Version => "0.1", 
        Date => "3 August 2020",
        Authors => {{Name=> "Carlos Amendola", 
	   Email=> "carlos.amendola@tum.de",
	   HomePage=>"http://www.carlos-amendola.com/"},
       
	  {Name => "Luis David Garcia Puente",
	   Email => "lgarcia@shsu.edu",
	   HomePage => "http://www.shsu.edu/~ldg005"},
       
          {Name=> "Roser Homs Pons", 
	   Email=> "roserhp@gmail.com",
	   HomePage=>"https://personal-homepages.mis.mpg.de/homspons/index.html"},
       
          {Name=> "Olga Kuznetsova", 
	   Email=> "kuznetsova.olga@gmail.com",
	   HomePage=>"https://okuznetsova.com"},
       
          {Name=> "Harshit J Motwani", 
	   Email=> "harshitmotwani2015@gmail.com"}},
        Headline => "Graphs specific for algebraic statistics",
        DebuggingMode => true,
	PackageExports => {"Graphs"}
        )

export {
    "graphComponents",
    "graphFunctions",
     "topologicalSort",
    "topSort",
    "SortedDigraph",
    "Bigraph",
    "bigraph",
    "LabeledGraph",
    "labeledGraph",
    "MixedGraph",
    "mixedGraph",
    "newDigraph",
    "collateVertices",
    "partitionLMG",
    "isLoopless",
    "hasDirCycles",
    "hasMultipleEdges"
    }

if Graphs.Options.Version < "0.3.2" then error "StatGraphs requires Graphs version 0.3.2 or later"

topologicalSort = method(TypicalValue =>List)
topologicalSort Digraph := List => D -> topologicalSort(D, "")
topologicalSort (Digraph, String) := List => (D,s) -> (
    if instance(D, Graph) or isCyclic D then error "Topological sorting is only defined for acyclic directed graphs.";
    s = toLower s;
    processor := if s == "random" then random
        else if s == "min" then sort
        else if s == "max" then rsort
        else if s == "degree" then L -> last \ sort transpose {apply(L, v -> degree(D, v)), L}
        else identity;
    S := processor sources D;
    L := {};
    v := null;
    while S != {} do (
        v = S_0;
        L = append(L, v);
        S = processor join(drop(S, 1), select(toList children (D, v), c -> isSubset(parents(D, c), L)));
        );
    L
    )




SortedDigraph = new Type of HashTable;

-- Keys:
--      digraph: the original digraph
--      NewDigraph: the digraph with vertices labeled as integers obtained from sorting
--      map: the map giving the sorted order

topSort = method(TypicalValue =>HashTable)
topSort Digraph := SortedDigraph => D -> (
L := topologicalSort D;
g := graph D;
new SortedDigraph from {
digraph => D,
newDigraph => digraph hashTable apply(#L, i -> i + 1 => apply(toList g#(L_i), j -> position(L, k -> k == j) + 1)),
map => hashTable apply(#L, i -> L_i => i + 1)
}
)


Bigraph = new Type of Graph

bigraph = method(TypicalValue =>Bigraph, Options => {Singletons => null, EntryMode => "auto"})
bigraph HashTable := opts -> g -> new Bigraph from graph(g, opts)
bigraph List := opts -> L -> new Bigraph from graph(L, opts)
bigraph (List, List):= opts -> (V,L) -> new Bigraph from graph(V,L, opts)
bigraph (List, Matrix) :=  opts -> (V,A) -> new Bigraph from graph(V,A, opts)
bigraph Matrix := opts -> A -> new Bigraph from graph(A, opts)


graphData = "graphData"
labels = "labels"

LabeledGraph = new Type of HashTable

labeledGraph = method(TypicalValue =>LabeledGraph)
labeledGraph (Digraph,List) := (g,L) -> (
    C := new MutableHashTable;
    C#cache = new CacheTable from {};
    lg := new MutableHashTable;
    lg#graphData = g;
    label := new MutableHashTable;
    if instance(g,Graph) then (
        sg := simpleGraph g;
        scan(L, i -> 
            if (sg#graph#(i#0#0))#?(i#0#1) then label#(i#0) = i#1
            else if (sg#graph#(i#0#1))#?(i#0#0) then label#({i#0#1,i#0#0}) = i#1
            else error (toString(i#0)|" is not an edge of the graph");
            );
        )
    else (
        scan(L, i -> 
            if (g#graph#(i#0#0))#?(i#0#1) then label#(i#0) = i#1
            else error (toString(i#0)|" is not an edge of the graph");
            );
        );
    lg#labels = new HashTable from label;
    C#graph = lg;
    new LabeledGraph from C
    )

net LabeledGraph := g -> horizontalJoin flatten (
     net class g,
    "{",
    stack (horizontalJoin \ sort apply(pairs (g#graph),(k,v) -> (net k, " => ", net v))),
    "}"
    )

toString LabeledGraph := g -> concatenate(
    "new ", toString class g#graph,
    if parent g#graph =!= Nothing then (" of ", toString parent g),
    " from {",
    if #g#graph > 0 then demark(", ", apply(pairs g#graph, (k,v) -> toString k | " => " | toString v)) else "",
    "}"
    )

graph LabeledGraph := opts -> g -> g#graph  --used to transform the LabeledGraph into a hashtable




MixedGraph = new Type of HashTable

mixedGraph = method(TypicalValue =>MixedGraph)
mixedGraph (Graph, Digraph, Bigraph) := (g,d,b) -> (
    C := new MutableHashTable;
    C#cache = new CacheTable from {};
    h := new MutableHashTable;
    h#Graph = g;
    h#Digraph = d;
    h#Bigraph = b;
    C#graph = new HashTable from h;
    new MixedGraph from C)
mixedGraph Graph := g -> mixedGraph(g,digraph {}, bigraph {})
mixedGraph Digraph := d -> mixedGraph(graph {},d, bigraph {})
mixedGraph Bigraph := b -> mixedGraph(graph {},digraph {}, b)
mixedGraph (Digraph, Bigraph) := (d,b) -> mixedGraph(graph {},d,b)
mixedGraph (Bigraph, Digraph) := (b,d) -> mixedGraph(graph {},d,b)
mixedGraph (Graph, Digraph) := (g,d) -> mixedGraph(g,d,bigraph {})
mixedGraph (Digraph, Graph) := (d,g) -> mixedGraph(g,d,bigraph {})
mixedGraph (Graph, Bigraph) := (g,b) -> mixedGraph(g,digraph {},b)
mixedGraph (Bigraph, Graph) := (b,g) -> mixedGraph(g,digraph {},b)
mixedGraph (Digraph, Bigraph, Graph) := (d,b,g) -> mixedGraph(g,d,b)
mixedGraph (Bigraph, Graph, Digraph) := (b,g,d) -> mixedGraph(g,d,b)
mixedGraph (Graph, Bigraph, Digraph) := (g,b,d) -> mixedGraph(g,d,b)
mixedGraph (Bigraph, Digraph, Graph) := (b,d,g) -> mixedGraph(g,d,b)
mixedGraph (Digraph, Graph, Bigraph) := (d,g,b) -> mixedGraph(g,d,b)

net MixedGraph := g -> horizontalJoin flatten (
     net class g,
    "{",
    stack (horizontalJoin \ sort apply(pairs (g#graph),(k,v) -> (net k, " => ", net v))),
    "}"
    )

toString MixedGraph := g -> concatenate(
    "new ", toString class g#graph,
    if parent g#graph =!= Nothing then (" of ", toString parent g),
    " from {",
    if #g#graph > 0 then demark(", ", apply(pairs g#graph, (k,v) -> toString k | " => " | toString v)) else "",
    "}"
    )

graph MixedGraph := opts -> g -> g#graph     	     --used to transform the MixedGraph into a hashtable
digraph MixedGraph := opts -> g -> g#graph#Digraph
bigraph MixedGraph := opts -> g -> g#graph#Bigraph
vertices MixedGraph := G -> toList sum(apply(keys(G#graph),i->set keys(graph (G#graph)#i)))
vertexSet MixedGraph := G -> vertices G

descendents (MixedGraph, Thing) := (G,v) -> descendents(digraph G, v)
nondescendents (MixedGraph, Thing) := (G,v) -> nondescendents(digraph G, v)
parents (MixedGraph, Thing) := (G,v) -> parents(digraph G, v)
foreFathers (MixedGraph, Thing) := (G,v) -> foreFathers(digraph G, v)
children (MixedGraph, Thing) := (G,v) -> children(digraph G, v)
neighbors (MixedGraph, Thing) := (G,v) -> neighbors(G#graph#Graph, v)
nonneighbors (MixedGraph, Thing) := (G,v) -> nonneighbors(G#graph#Graph, v)


collateVertices = method()
collateVertices MixedGraph := g -> (
    v := vertices g;
    hh := new MutableHashTable;
    G := graph g;
    -- Graph
    x := graph G#Graph;
    scan(v,j->if x#?j then hh#j=x#j else hh#j={});
    gg := graph(new HashTable from hh);
    -- Digraph
    x = graph G#Digraph;
    scan(v,j->if x#?j then hh#j=x#j else hh#j={});
    dd := digraph(new HashTable from hh);
    -- Bigraph
    x = graph G#Bigraph;
    scan(v,j->if x#?j then hh#j=x#j else hh#j={});
    bb := bigraph(new HashTable from hh);
    mixedGraph(gg,dd,bb))

indexLabelGraph MixedGraph := MixedGraph => G -> (
    V := vertices G;
    h := hashTable apply(#V, i -> V_i => i);
    U := G#graph#Graph;
    B := G#graph#Bigraph;
    D := G#graph#Digraph;
    
    inputG:=new MutableHashTable;
    inputG#graphComponents={U,B,D};
    inputG#graphFunctions={graph,bigraph,digraph};
    inputG=new HashTable from inputG;
    
    G=mixedGraph toSequence(
      for i to #inputG#graphComponents-1 list (
      E := apply(toList \ edges inputG#graphComponents_i, e -> {h#(e_0), h#(e_1)});
      inputG#graphFunctions_i(flatten E, E,EntryMode => "edges")
      )
     )
    );


-- Makes a partition U\cup W of the vertices V of a loopless mixed graph (inputed as a mixedGraph) 
-- such that U contains all the vertices adjacent to undirected edges, 
-- W contains all the vertices adjacent to bidirected edges 
-- and there are no directed edges from W to U
-- and all vertices in U have lower value than those in W.
partitionLMG = method()
partitionLMG MixedGraph := g -> (
   --check it's a simple graph
   if isLoopless(g)==false then print ("Warning: the expected input is a loopless mixed graph.");
   --retrieve graph, bigraph and digraph
   G:= g#graph#Graph;
   B:= g#graph#Bigraph;
   D:= g#graph#Digraph;
   --check there are no directed cycles
   if hasDirCycles g  then print ("Warning: the expected input is a loopless mixed graph without directed cycles.");
   --naive partition (vertices only adjacent to directed edges are not considered) 
   U:=vertices G;
   W:=vertices B;
   --check there are no common vertices for undirected and bidirected edges
   if not (set U * set W===set {}) then 
   error("Vertices cannot be adjacent to bidirected and undirected edges simultaneously.");
   --check that vertices in U (undir) have lower value than vertices in W (bidir):
   if max U > min W then error ("Vertex ordering issues: vertices adjacent to undirected edges should have lower value than vertices adjacent to bidirected edges");
   --check there are no directed edges from set U to set W
   for e in edges D do (if (member(e_0,set W) and member(e_1,set U)) 
   then error("Directed edges cannot go from vertices adjacent to a bidirected edge to vertices adjacent to an undirected edge"));
   --check directed edges always go from lower value to higher value
   for e in edges D do (
       if e_0 > e_1 then error ("Vertex ordering issues: directed edges must go from low to high");
   );
   --check whether there are remaining vertices (only adjacent to directed edges)
   V:=set vertices g-set W-set U;
   if V===set{} then return (U,W);
   --place remaining vertices in either U or W depending on their value
   for v in toList V do (if v < max U then U=append(U,v) else W=append(W,v););
   U,W
)

--Check whether a graph is loopless in each type of edges
isLoopless = method()
isLoopless MixedGraph := Boolean => g -> (
   --retrieve graph, digraph and bigraph
   G:= g#graph#Graph;
   B:= g#graph#Bigraph;
   D:= g#graph#Digraph;
   --retrieve underlying undirected edges 	
   --edG:=edges G;
   --edD:=edges underlyingGraph D;
   --edB:=edges B;
   --build list of underlying edges with and without repetitions
   --ed:=join(edG,edD,edB);
   --e:=unique(ed);
   --check there are no loops and no repetitions
   isSimple(G) and isSimple(underlyingGraph D) and isSimple(B) --and #ed==#e
   )

isLoopless Graph := Boolean => g -> (
   isSimple(g)
   )

isLoopless Bigraph := Boolean => g -> (
   isSimple(g)
   )

isLoopless Digraph := Boolean => g -> (
   isSimple(underlyingGraph g)
   )

-- Checked whether a MixedGraph contains multiple edges
hasMultipleEdges = method()
hasMultipleEdges MixedGraph := Boolean => g -> (
   --retrieve graph, digraph and bigraph
   G:= g#graph#Graph;
   B:= g#graph#Bigraph;
   D:= g#graph#Digraph;
   --retrieve underlying undirected edges 	
   edG:=edges G;
   edD:=edges underlyingGraph D;
   edB:=edges B;
   --build list of underlying edges with and without repetitions
   d:=join(edG,edD,edB);
   e:=unique(d);
   --check there are no repetitions
   #d==#e
   )

-- Check whether a MixedGraph does not contain any directed cycles
hasDirCycles=method()
hasDirCycles MixedGraph := Boolean => g -> (
    G:=indexLabelGraph g;
    U:= graph(sort vertices G#graph#Graph,edges G#graph#Graph);
    B:= bigraph(sort vertices G#graph#Bigraph,edges G#graph#Bigraph);
    D:= digraph(vertices G,edges G#graph#Digraph);
    compU:=connectedComponents U;
    compB:=connectedComponents B;
    vertOnlyDir:=vertices D - set vertices U - set vertices B;
    allComp:=flatten  {connectedComponents U,connectedComponents B, pack(vertOnlyDir,1)};
    n:=# compU + # compB + #vertOnlyDir;
    if n==1 then not edges D=={} else(
    adjMG:=mutableMatrix(ZZ,n,n);
    -- form the adjacency matrix of the graph of chain components 
    for i from 0 to  n - 1 do
	   (
	       for j from 0 to  n - 1 do
	       (if not submatrix(adjacencyMatrix D, toList(set(allComp_i)*set(vertices D)), toList(set(allComp_j)*set(vertices D)))==0 then adjMG_(i,j)=1 else adjMG_(i,j)=0
	       ));

    adjMG=matrix adjMG;
    isCyclic (digraph adjMG))  
    );

--******************************************--
-- DOCUMENTATION     	       	    	    -- 
--******************************************--

beginDocumentation()

doc ///
    Key
        StatGraphs
    Headline
        a package for graphs that are used primarily in statistical models.
    Description        
        Text
            This package contains the types of graphs that are used in algebraic 
      	    statistics: @TO Bigraph@,  @TO LabeledGraph@,  @TO MixedGraph@ and  @TO SortedDigraph@. 
	    
	    A Bigraph is a simple graph with bidirected edges. 
	    A MixedGraph is a graph with undirected, directed and bidirected edges. 
	    
        Example	   
	    G = bigraph {{3,4},{1,2},{2,4}}
	    	    
        Example	   
	    G = mixedGraph(graph{{1,2}},digraph {{1,3},{2,3}},bigraph {{3,4}})

	    
   	
    Caveat
       StatGraphs requires  @TO Graphs@ version 0.3.2 or later.
       
    SeeAlso
       Graphs
    	   
///

--------------------------------
-- Documentation
--------------------------------

--------------------------------------------
-- Documentation Bigraph
--------------------------------------------

doc ///
  Key 
     Bigraph

  Headline
     a simple graph with  bidirected edges. 
  Description
     Text  
         Bigraph is a simple graph that has  bidirected edges.
	 To create a Bigraph, use @TO bigraph@. 
	
  SeeAlso
     bigraph

///

--------------------------------------------
-- Documentation bigraph
--------------------------------------------

doc ///
  Key
     bigraph
     (bigraph, HashTable)
     (bigraph, List)
     (bigraph, List, List)
     (bigraph, List, Matrix)
     (bigraph, Matrix)

    
  Headline
     this function creates a Bigraph 
  Usage
     G= bigraph(H) 
     G= bigraph(L) 
     G= bigraph(V,L)
     G= bigraph(V,A)
     G= bigraph(A)   
 
  Inputs
     H:HashTable
       hashtable of edges	 	 	 
     L:List
       list of edges
     V:List 
       list of vertices
     A:Matrix
       adjacency matrix   
    
  Outputs
     :Bigraph
    
  Description
    Text
        This is a constructor of a simple graph of type Bigraph.  One can use the same input types
	as in @TO graph@. 
    Example
        G = bigraph {{3,4},{1,2},{2,4}}

  SeeAlso
    Bigraph

///

--------------------------------------------
-- Documentation MixedGraph
--------------------------------------------

doc ///
  Key 
     MixedGraph

  Headline
     a graph that has undirected, directed and bidirected edges. 
  Description
     Text  
         MixedGraph is a graph that has undirected, directed and bidirected edges.
	 To create a MixedGraph, use @TO mixedGraph@. Each subgraph (undirected,
	 directed and bidirected) is a simple graph.
	
  SeeAlso
     mixedGraph
     collateVertices

///

--------------------------------------------
-- Documentation mixedGraph
--------------------------------------------

doc ///
  Key
     mixedGraph
     (mixedGraph, Graph, Digraph, Bigraph)
     (mixedGraph, Graph,Bigraph,Digraph) 
     (mixedGraph, Digraph,Graph,Bigraph)
     (mixedGraph, Digraph,Bigraph,Graph)
     (mixedGraph, Bigraph,Graph,Digraph)
     (mixedGraph, Bigraph,Digraph,Graph) 
     (mixedGraph, Graph, Digraph)
     (mixedGraph, Digraph,Graph)
     (mixedGraph, Digraph, Bigraph)
     (mixedGraph, Bigraph,Digraph)
     (mixedGraph, Graph, Bigraph)
     (mixedGraph, Bigraph, Graph)
     (mixedGraph, Graph)
     (mixedGraph, Digraph)
     (mixedGraph, Bigraph)
    
  Headline
     this function creates a MixedGraph from a combination of Graph, Digraph and Bigraph 
  Usage
     G= mixedGraph(U, D, B) 
     G= mixedGraph G 
     G= mixedGraph D
     G= mixedGraph B
     G= mixedGraph(G, D)
     G= mixedGraph(G, B)
     G= mixedGraph(D, B)
 
  Inputs
     U:Graph
       component that contains all undirected edges of the graph	 	 	 
     D:Digraph
       component that contains all directed edges of the graph
     B:Bigraph
       component that contains all bidirected edges of the graph
    
  Outputs
     :MixedGraph
    
  Description
    Text
        This is a constructor of graphs of type MixedGraph from a combination of Graph, Digraph and Bigraph.  One can also input any subset 
        and any permutation of the arguments.
    
        Note that this constructor does not check the input satisfies the properties of loopless mixed graphs from 
        Sadeghi and Lauritzen, 2020 <@HREF"https://arxiv.org/pdf/1109.5909.pdf"@>.
      
    Example
        G = mixedGraph(graph{{1,2}},digraph {{1,3},{2,3}},bigraph {{3,4}})

  SeeAlso
    MixedGraph

///

--------------------------------------------
-- Operations on MixedGraph
--------------------------------------------

--------------------------------------------
-- Documentation bigraph(MixedGraph)
--------------------------------------------

doc ///
  Key
     (bigraph, MixedGraph)
  Headline
     Extract the Bigraph component of a MixedGraph
  Usage
     G=bigraph G
 
  Inputs
     G:MixedGraph
    
  Outputs
     :Bigraph
    
  Description
    Text
        This method extracts the Bigraph component of a MixedGraph
      
    Example
        G= mixedGraph(graph{{a,b},{b,c}},digraph {{a,d},{c,e},{f,g}},bigraph {{d,e}})
        bigraph G

  SeeAlso
    MixedGraph
    (digraph, MixedGraph) 
    (graph, MixedGraph) 

///

--------------------------------------------
-- Documentation digraph(MixedGraph)
--------------------------------------------

doc ///
  Key
     (digraph, MixedGraph)
  Headline
     Extract the Digraph component of a MixedGraph
  Usage
     G=digraph G
 
  Inputs
     G:MixedGraph
    
  Outputs
     :Digraph
    
  Description
    Text
        This method extracts the Digraph component of a MixedGraph
      
    Example
        G= mixedGraph(graph{{a,b},{b,c}},digraph {{a,d},{c,e},{f,g}},bigraph {{d,e}})
        digraph G

  SeeAlso
    MixedGraph
    (bigraph, MixedGraph) 
    (graph, MixedGraph) 

///

--------------------------------------------
-- Documentation collateVertices 
--------------------------------------------

doc /// 
    Key
        collateVertices
        (collateVertices, MixedGraph) 
	
    Headline
        converts a MixedGraph so that each component subgraph has the same set of vertices.
    Usage
        collateVertices(G)
    Inputs
        G:MixedGraph
    Outputs
         :MixedGraph 
    Description 
        Text
	    Let G=mixedGraph(U,D,B) and denote the vertices of U by V1, 
	    the vertices of D by V2 and the vertices of B by  V3.
	    Then the method collateVertices(G) outputs a mixedGraph with same 
	    edges as before but with V1 \cup V2 \cup V3 as the vertices of U,D 
	    and B.
	    
        Example
	   U = graph{{1,2},{2,3},{3,4},{1,4},{1,5}}
	   D = digraph{{2,1},{3,1},{7,8}}
	   B = bigraph{{1,5}}
	   G = mixedGraph(U,D,B)
	   collateVertices G
   ///
--------------------------------------------
-- Documentation indexLabelGraph
--------------------------------------------

doc ///
  Key
     (indexLabelGraph, MixedGraph)
  Headline
     Relabels the vertices of a MixedGraph according to their indices, indexed from 0.
  Usage
     G=indexLabelGraph G
 
  Inputs
     G:MixedGraph
    
  Outputs
     :MixedGraph
    
  Description
    Text
        This method relabels the vertices of a MixedGraph according to their indices. The method indexes from 0 to the number of vertices minus one.
	This is an adaptation of @TO indexLabelGraph@.
      
    Example
        G= mixedGraph(graph{{a,b},{b,c}},digraph {{a,d},{c,e},{f,g}},bigraph {{d,e}})
	indexLabelGraph G

  SeeAlso
    MixedGraph
    indexLabelGraph

///

--------------------------------------------
-- Documentation isLoopless 
--------------------------------------------

doc /// 
    Key
        isLoopless 
        (isLoopless, MixedGraph) 
	(isLoopless, Graph) 
	(isLoopless, Bigraph) 
	(isLoopless, Digraph) 
    Headline
        checks whether a graph contains a loop
    Usage
        isLoopless(G)
    Inputs
        G:
	 graph of Type MixedGraph, Graph or Digraph
    Outputs
         :Boolean 
    Description 
        Text
	  This method checks whether a graph contains a loop. 
	  
	  If the input is a  @TO Graph@ or a  @TO Bigraph@, then this is equivalent to 
	   @TO isSimple@. 
	  
	  If the input is a  @TO Digraph@, then this is equivalent to checking whether 
	  the  @TO underlyingGraph@ @TO isSimple@.
	  
	  If the input is a  @TO MixedGraph@, then this checks whether the undirected,
	  directed and bidirected subgraphs separately contain loops.
        Example
	   U = graph{{1,2},{2,3},{3,4}}
	   D = digraph{{2,5}}
	   B = bigraph{{5,6}}
	   G = mixedGraph(U,D,B)
	   isLoopless G
	   
	Example
	   U = graph{{1,1}}
	   isLoopless U   
   ///

--------------------------------------------
-- Documentation noDirCycles 
--------------------------------------------

doc /// 
    Key
        hasDirCycles 
        (hasDirCycles, MixedGraph) 
	
    Headline
        checks whether a MixedGraph contains a directed cycle
    Usage
        hasDirCycles(G)
    Inputs
        G:MixedGraph
    Outputs
         :Boolean 
    Description 
        Text
	   This method checks whether a @TO MixedGraph@ a directed cycle. 
	   
	   A directed cycle is a cycle in the @TO Digraph@ constructed from G by
	   identifying all connected components on bidirected and undirected edges.
	   Such a connected component contain either bidirected edges only or 
	   undirected edges only.
	    
        Example
	   U = graph{{1,2},{2,3},{3,4},{1,4},{1,5}}
	   D = digraph{{2,1},{3,1},{7,8}}
	   B = bigraph{{1,5}}
	   G = mixedGraph(U,D,B)
	   hasDirCycles G
        Example
	   U = graph{{1,2},{3,4}}
	   D = digraph{{1,3},{4,2}}
	   G = mixedGraph(U,D)
	   hasDirCycles G 
	Example
	   U = graph{{1,2}}
	   B = bigraph{{3,4}}
	   D = digraph{{1,3},{4,2}}
	   G = mixedGraph(U,D,B)
	   hasDirCycles G      
   ///
--------------------------------------------
-- Operations on vertices of a  MixedGraph
--------------------------------------------  

--------------------------------------------
-- Documentation children
--------------------------------------------
doc ///
  Key
     (children,MixedGraph,Thing)
    
  Headline
    returns the children of a vertex of a MixedGraph
  Usage
    children (G,v)
 
  Inputs
    G: MixedGraph
    v: Thing
       a vertex of G
    
  Outputs
     :Set
    
  Description
    Text
        The children of v are the all the vertices u such that v,u is in the edge 
	set of the @TO MixedGraph@ G. So the children of a vertex v are exactly those 
	vertices on a MixedGraph graph that v points to.
    Example
        G = mixedGraph(graph{{3,1}},digraph {{1,2},{2,3}},bigraph {{3,4},{2,4}})
	children (G,1)
	children (G,2)
    	children (G,3)

  SeeAlso
    MixedGraph

///
 
--------------------------------------------
-- Documentation vertices and vertexSet
--------------------------------------------

doc ///
  Key
     (vertices, MixedGraph)
     (vertexSet, MixedGraph)
    
  Headline
     this function creates a union of all the vertices of Graph, Bigraph, Digraph components of a MixedGraph.
  Usage
     V=vertices G
     V=vertexSet G
 
  Inputs
     G:MixedGraph
    
  Outputs
     :List
    
  Description
    Text
        This function creates a union of all the vertices of Graph, Bigraph, Digraph components of a MixedGraph.
	This is an adaptation of vertices and vertexSet from @TO Graphs@.
      
    Example
        G = mixedGraph(graph{{3,1}},digraph {{1,2},{2,3}},bigraph {{3,4}})
	vertices G
	vertexSet G

  SeeAlso
    MixedGraph

///
--------------------------------------------
-- Topological sorting functions
--------------------------------------------

--------------------------------------------
-- Documentation topologicalSort
--------------------------------------------

doc /// 
    Key
        topologicalSort
        (topologicalSort, Digraph) 
	(topologicalSort, Digraph, String) 
    Headline
        outputs a list of vertices in a topologically sorted order of a DAG.
    Usage
        topologicalSort(D)
    Inputs
        D:Digraph
	S:String
    Outputs
         :List 
    Description 
        Text
	    This function outputs a list of vertices in a topologically sorted order of a directed acyclic graph (DAG).
        Example
	   D = digraph{{2,1},{3,1}}
	   topologicalSort D
   ///

--------------------------------------------
-- Documentation topSort
--------------------------------------------
doc /// 
    Key
        topSort
        (topSort, Digraph) 
    Headline
        outputs a list of vertices in a topologically sorted order of a DAG.
    Usage
        topologicalSort(D)
    Inputs
        D:Digraph
	  needs to be a directed acyclice graph DAG  
    Outputs
         :HashTable 
    Description 
        Text
	    This method outputs a HashTable with keys digraph, map and newDigraph, where digraph is the original digraph,
	    map is the relation between old ordering and the new ordering of vertices and newDigraph is the Digraph with 
	    topologically sorted vertices.

        Example
	   D = digraph{{2,1},{3,1}}
	   H = topSort D
	   H#digraph
	   H#map
	   H#newDigraph
   ///

--******************************************--
-- TESTS     	       	    	      	    --
--******************************************--

--------------------------------------
--------------------------------------
end--
--------------------------------------
--------------------------------------

