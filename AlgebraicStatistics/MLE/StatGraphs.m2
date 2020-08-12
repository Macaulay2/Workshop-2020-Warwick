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
     "topologicalSort",
    "topSort",
    "SortedDigraph",
    "Bigraph",
    "bigraph",
    "LabeledGraph",
    "labeledGraph",
    "MixedGraph",
    "mixedGraph",
    "collateVertices"
    }


--check whats the difference between topologicalSort and topSort
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

bigraph = method(Options => {Singletons => null})
bigraph HashTable := opts -> g -> new Bigraph from graph(g, opts)
bigraph List := opts -> g -> new Bigraph from graph(g, opts)


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

graph LabeledGraph := opts -> g -> g#graph




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
mixedGraph (Graph, Digraph) := (g,d) -> mixedGraph(g,d,bigraph {})
mixedGraph (Graph, Bigraph) := (g,b) -> mixedGraph(g,digraph {},b)

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

graph MixedGraph := opts -> g -> g#graph
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


beginDocumentation()

doc ///
  Key
   StatGraphs
  Headline
    Graphs specific for algebraic statistics
  Description
    Text
      This package contains the types of graphs that are used in algebraic 
      statistics: bidirected graphs and mixed graphs. Mixed graphs are loopless
      graphs with directed, bidirected and undirected edges.
      
      !!!TBD: talk more about the functions!!!
    Example
     
  SeeAlso
   Graphs
   GraphicalModels
   GraphicalModelsMLE 
///

--------------------------------------------
-- Documentation MixedGraph
--------------------------------------------

doc ///
  Key 
  MixedGraph

  Headline
  MixedGraph is a graph that has undirected, directed and bidirected edges. 
  
  Description
  Text  
  To create a MixedGraph, use mixedGraph (Graph, Digraph, Bigraph). MixedGraph is a type, 
  with ancestor classes HashTable < Thing. 
  
  Common ways to use MixedGraph are the following functions in the GraphicalModels package:
   - gaussianRing(MixedGraph)
   - trekIdeal(Ring,MixedGraph)
   - trekSeparation(MixedGraph)
   - collateVertices(MixedGraph)
  
  NOTE 1: While the type MixedGraph allows for any combination of undirexted, directed and bidirected edges, the functions
  using MixedGraph as input make the following assumptions (in line with loopless mixed graphs from 
      Sadeghi and Lauritzen, 2020 <@HREF"https://arxiv.org/pdf/1109.5909.pdf"@>):
  The vertices of $G$ can be partitioned into $V=U\cup W$ such that
  - if $i-j$ in $G$ then $i,j\in U$
  - if $i\leftrightarrow j$ in $G$ then $i,j\in W$ 
  - if $i\to j$ in $G$ then $i\in U$ and $j\in W$ 
 and
 the directed edges form a DAG, i.e., there are no directed cycles. Here, the nodes of a directed cycles can
 be
 - the vertices $V$ of  $G$  
 - undirected connected components 
 - bidirected connected components
 For example,   $\{1 \to 2,2\leftrightarrow 3,3 \to 1\}$ is a directed cycle.
 
 NOTE 2: Current implementations of trekIdeal and trekSeparation do not allow undirected edges.
 
 To extract information about MixedGraph use:
 - graph MixedGraph
 - digraph MixedGraph
 - bigraph MixedGraph
 - vertices MixedGraph
 
 - descendents (MixedGraph, vertex)
 - parents (MixedGraph, vertex)
 - foreFathers (MixedGraph, vertex) 
 - children (MixedGraph, vertex) 
 - neighbors (MixedGraph, vertex) 
 - nonneighbors (MixedGraph, vertex)    

  SeeAlso
  mixedGraph, gaussianRing, trekIdeal, trekSeparation, Graph, Digraph, Bigraph, collateVertices

///

--------------------------------------------
-- Documentation mixedGraph
--------------------------------------------

doc ///
  Key
  mixedGraph
  (mixedGraph, Graph, Digraph, Bigraph)
  (mixedGraph, Graph, Digraph)
  (mixedGraph, Digraph, Bigraph)
  (mixedGraph, Graph, Bigraph)
  (mixedGraph, Graph)
  (mixedGraph, Digraph)
  (mixedGraph, Bigraph)
    
  Headline
  This function creates a MixedGraph from a combination of Graph, Digraph and Bigraph. One can also input any subset 
  of the arguments.
   
  Usage
  G= mixedGraph(U, D, B) 
  G= mixedGraph G 
  G= mixedGraph D
  G= mixedGraph B
  G= mixedGraph(G, D)
  G= mixedGraph(G, B)
  G= mixedGraph(D, B)
 
  Inputs
  U: Graph
  D: Digraph
  B: Bigraph
    
  Outputs
  G: MixedGraph
      The graph with undirected edges given by U, directed edges given by D and bidirected edges given by B.
    
  Description
    Text
    This is a constructor of graphs of type MixedGraph. Note that this constructor does not check the input
    satisfies the properties of loopless mixed graphs from 
    Sadeghi and Lauritzen, 2020 <@HREF"https://arxiv.org/pdf/1109.5909.pdf"@>.
      
    Example
    G = mixedGraph(graph{{1,2}},digraph {{1,3},{2,3}},bigraph {{3,4}})

  SeeAlso
  MixedGraph

///

end--
doc ///
  Key

  Headline
   
  Usage
 
  Inputs
    
  Outputs
    
  Description
    Text
      
    Example

  SeeAlso

///


TEST ///
 
///



