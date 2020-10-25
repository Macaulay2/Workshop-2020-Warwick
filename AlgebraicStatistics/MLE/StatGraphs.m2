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
    "isMixedGraphLoopless",
    "indexLabelMixedGraph",
    "noDirCycles"
    }


topologicalSort = method()
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

topSort = method()
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

bigraph = method(Options => {Singletons => null, EntryMode => "auto"})
bigraph HashTable := opts -> g -> new Bigraph from graph(g, opts)
bigraph List := opts -> L -> new Bigraph from graph(L, opts)
bigraph (List, List):= opts -> (V,L) -> new Bigraph from graph(V,L, opts)
bigraph (List, Matrix) :=  opts -> (V,A) -> new Bigraph from graph(V,A, opts)
bigraph Matrix := opts -> A -> new Bigraph from graph(A, opts)


graphData = "graphData"
labels = "labels"

LabeledGraph = new Type of HashTable

labeledGraph = method()
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

mixedGraph = method()
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


indexLabelMixedGraph = method()
indexLabelMixedGraph MixedGraph := MixedGraph => G -> (
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
   if isMixedGraphLoopless(g)==false then error ("The input should be a simple mixedGraph.");
   --retrieve graph, bigraph and digraph
   G:= g#graph#Graph;
   B:= g#graph#Bigraph;
   D:= g#graph#Digraph;
   --check there are no directed cycles
   if not noDirCycles g  then error ("A loopless mixed graph should not contain directed cycles.");
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

--Check whether a MixedGraph is loopless in each type of edges
isMixedGraphLoopless = method()
isMixedGraphLoopless MixedGraph := Boolean => g -> (
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

-- Check whether a Mixed Graph does not contain any directed loops
noDirCycles=method()
noDirCycles MixedGraph := Boolean => g -> (
    G:=indexLabelMixedGraph g;
    U:= graph(sort vertices G#graph#Graph,edges G#graph#Graph);
    B:= bigraph(sort vertices G#graph#Bigraph,edges G#graph#Bigraph);
    D:= digraph(vertices G,edges G#graph#Digraph);
    compU:=connectedComponents U;
    compB:=connectedComponents B;
    vertOnlyDir:=vertices D - set vertices U - set vertices B;
    allComp:=flatten  {connectedComponents U,connectedComponents B, pack(vertOnlyDir,1)};
    n:=# compU + # compB + #vertOnlyDir;
    if n==1 then edges D=={} else(
    adjMG:=mutableMatrix(ZZ,n,n);
    -- form the adjacency matrix of the graph of chain components 
    for i from 0 to  n - 1 do
	   (
	       for j from 0 to  n - 1 do
	       (if not submatrix(adjacencyMatrix D, toList(set(allComp_i)*set(vertices D)), toList(set(allComp_j)*set(vertices D)))==0 then adjMG_(i,j)=1 else adjMG_(i,j)=0
	       ));

    adjMG=matrix adjMG;
    not isCyclic (digraph adjMG))  
    );


