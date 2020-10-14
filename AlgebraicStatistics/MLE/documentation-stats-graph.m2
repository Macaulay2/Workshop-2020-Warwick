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
