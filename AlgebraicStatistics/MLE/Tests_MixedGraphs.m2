restart
loadPackage "GraphicalModels"
G=graph{{1,3}}
D=digraph{{2,1},{3,2}}
B=bigraph{{3,4}}

LMG=mixedGraph(G,D,B)
R=gaussianRing(LMG)
describe R
undirectedEdgesMatrix(R)
directedEdgesMatrix(R)
bidirectedEdgesMatrix(R)
covarianceMatrix(R)

----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
restart
loadPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{5,3}}
B=bigraph{{5,4}}

LMG=mixedGraph(G,D,B)
R=gaussianRing(LMG)
describe R
undirectedEdgesMatrix(R)
directedEdgesMatrix(R)
bidirectedEdgesMatrix(R)
covarianceMatrix(R)

----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
restart
loadPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{1,3},{3,2}}
B=bigraph{{5,4}}

g=mixedGraph(G,D,B)
R=gaussianRing(LMG)
describe R
undirectedEdgesMatrix(R)
directedEdgesMatrix(R)
bidirectedEdgesMatrix(R)
covarianceMatrix(R)

connectedComponentsMG = method()
connectedComponentsMG MixedGraph := List => G -> (
    V := vertices G;
    while #V != 0 list (
        C := {first V};
        i := 0;
        while i!= #C do (
            N := toList neighborsMG2(G, C_i);
            C = unique(C | N);
            V = V - set C;
            i = i + 1;
            if #V == 0 then break;
            );
        C
        )
    )

V=vertices LMG

neighborsMG(LMG,1)

neighborsMG(LMG,2)

neighborsMG(LMG,3)

neighborsMG(LMG,4)

neighborsMG(LMG,5)

neighborsMG2(LMG,1)

neighborsMG2(LMG,2)

neighborsMG2(LMG,3)

neighborsMG(LMG,4)

neighborsMG(LMG,5)

neighborsMG(LMG,10)

connectedComponentsMG LMG





restart
loadPackage "GraphicalModels"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

g=mixedGraph(G,D,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

isSubset((set vertexSet G)*(set vertexSet B),set {})

L=connectedComponentsMG g

U={};
W={};
for l in L do (
i=0;
while(member(l_i,vertexSet G)==false and member(l_i,vertexSet B)==false and i==#l) do i=i+1;
    if i==#l then (if U=={} then W=flatten append(W,l) else U=flatten append(U,l); continue;);
    if member(l_i,vertexSet G)==true then U=flatten append(U,l) else W=flatten append(W,l);
    ) 

i=2

i==#l_i

U={1,2,3}
U=={}
U,W


-----------------------------------------------------------
-----------------------------------------------------------
connectedComponentsMG2 = method()
connectedComponentsMG2 MixedGraph := List => G -> (
    V := vertices G;
    while #V != 0 list (
        C := {first V};
        i := 0;
        while i!= #C do (
            N := toList neighborsMG2(G, C_i);
            C = unique(C | N);
            V = V - set C;
            i = i + 1;
            if #V == 0 then break;
            );
        C
        )
    )

-----------------------------------------------------------
-----------------------------------------------------------
connectedComponentsMG = method()
connectedComponentsMG MixedGraph := List => G -> (
    V := vertices G;
    while #V != 0 list (
        C := {first V};
        i := 0;
        while i!= #C do (
            N := toList neighborsMG(G, C_i);
            C = unique(C | N);
            V = V - set C;
            i = i + 1;
            if #V == 0 then break;
            );
        C
        )
    )
-----------------------------------------------------------
-----------------------------------------------------------
-- neighbors considering directed edges as undirected edges
neighborsMG = method()
neighborsMG(MixedGraph, Thing) := Set => (g,v) -> (
    i := position(vertices g, u -> u === v);
    if i === null then error "v is not a vertex of g.";
    G:= collateVertices g;
    S=set {};
    for h in keys G#graph do S=S+neighbors(underlyingGraph G#graph#h,v);
    S
 )
-----------------------------------------------------------
-----------------------------------------------------------

-----------------------------------------------------------
-----------------------------------------------------------
-- neighbors considering directed edges
neighborsMG2 = method()
neighborsMG2 (MixedGraph, Thing) := Set => (g,v) -> (
    i := position(vertices g, u -> u === v);
    if i === null then error "v is not a vertex of g.";
    G:= collateVertices g;
    V={};
    for h in keys G#graph do V=append(V,(vertexSet G#graph#h)_(positions(first entries (adjacencyMatrix G#graph#h)^{i}, j -> j != 0)));
    set flatten unique V
 )
-----------------------------------------------------------
-----------------------------------------------------------
pos := (h, x) -> position(h, i->i===x)
gaussianRingList := new MutableHashTable;


gaussianRingMG = method(Dispatch=>Thing, Options=>{Coefficients=>QQ, sVariableName=>"s", lVariableName=>"l", 
	  pVariableName=>"p", kVariableName=>"k"})
gaussianRingMG MixedGraph := Ring => opts -> (g) -> (
     -- convert mixedGraph to hash table
     gg:= graph g;
     -- necessary condition (not sufficient) to have partition V=U\union W
     if(isSubset((set vertexSet gg#Graph)*(set vertexSet gg#Bigraph),set {})==false) then error "undirected and bidirected edges in same connected component";
     -- make partition V=U\union W
     L:=connectedComponentsMG g;
     U:={};
     W:={};
     for l in L do (
	 i=0;
	 while(member(l_i,vertexSet gg#Graph)==false and member(l_i,vertexSet gg#Bigraph)==false and i==#l) do i=i+1;
    	 if i==#l then (if U=={} then W=flatten append(W,l) else U=flatten append(U,l); continue;);
    	 if member(l_i,vertexSet gg#Graph)==true then U=flatten append(U,l) else W=flatten append(W,l);
    	 );
     -- sort vertices only according to vertex number
     vv := sort vertices g;
     -- add all vertices to all graphs and convert them to hash tables
     G := graph collateVertices g;
     dd := graph G#Digraph;
     bb := graph G#Bigraph;
     uu := graph G#Graph;
     -- set ring variables
     s := toSymbol opts.sVariableName;
     l := toSymbol opts.lVariableName;
     p := toSymbol opts.pVariableName;
     k := toSymbol opts.kVariableName; --ADDED
     kk := opts.Coefficients;          
     if (not gaussianRingList#?(kk,s,k,l,p,vv)) then ( --k ADDED
	  --(kk,s,k,l,p,vv) uniquely identifies gaussianRing in case of MixedGraph input.
     sL := delete(null, flatten apply(vv, x-> apply(vv, y->if pos(vv,x)>pos(vv,y) then null else s_(x,y))));
     -- start adding from gaussianRing Graph
     kL := join(apply(U, i->k_(i,i)),delete(null, flatten apply(U, x-> apply(toList uu#x, y->if pos(vv,x)>pos(vv,y) then null else k_(x,y)))));
     -- m := #kL; --eliminate the k's 
     -- end adding
     lL := delete(null, flatten apply(vv, x-> apply(toList dd#x, y->l_(x,y))));	 
     pL := join(apply(W, i->p_(i,i)),delete(null, flatten apply(W, x-> apply(toList bb#x, y->if pos(vv,x)>pos(vv,y) then null else p_(x,y)))));
     m := #kL+#lL+#pL; --#kL ADDED
     R := kk(monoid [kL,lL,pL,sL,MonomialOrder => Eliminate m, MonomialSize=>16]); --kL ADDED
     -- create gaussianVariables hash table: (symbol s)_(i,j) => ring var with the same name, same for l, p.
     H := new MutableHashTable;
     nextvar := 0;
     for v in kL do (H#v = R_nextvar; nextvar = nextvar+1); --ADDED
     for v in lL do (H#v = R_nextvar; nextvar = nextvar+1);
     for v in pL do (H#v = R_nextvar; nextvar = nextvar+1);
     for v in sL do (H#v = R_nextvar; nextvar = nextvar+1);
     R.gaussianVariables = new HashTable from H;
     R#numberOfEliminationVariables = m;
     R.gaussianRingData = {#vv,s,k,l,p}; -- k ADDED
     R.mixedGraph = g;
     gaussianRingList#((kk,s,k,l,p,vv)) = R;); -- k ADDED
     gaussianRingList#((kk,s,k,l,p,vv)) -- k ADDED
     )

----------------------------------------------------------------------------------
----------------------------------------------------------------------------------


--undirected graph
LMG=mixedGraph(G)
R=gaussianRing(LMG)
describe R
undirectedEdgesMatrix(R)
--comparison with gaussianRing Graph
S=gaussianRing(G)
describe S
undirectedEdgesMatrix(S)
covarianceMatrix(S)

--directed graph
LMG=mixedGraph(D)
R=gaussianRing(LMG)
describe R
directedEdgesMatrix(R)
covarianceMatrix(R)
--comparison with gaussianRing Digraph
S=gaussianRing(D)
describe S
covarianceMatrix(S)
directedEdgesMatrix(S) -- ONLY WORKS FOR GAUSSIAN RINGS COMING FROM MIXED GRAPHS

--undirected and directed edges
LMG2=mixedGraph(D,B)
R2=gaussianRing(LMG2)
describe R2

--directed and bidirected
LMG=mixedGraph(D,B)
R=gaussianRing(LMG)
describe R
covarianceMatrix(R)
directedEdgesMatrix(R)
bidirectedEdgesMatrix(R)
----------------------------------------------------------------------------------
----------------------------------------------------------------------------------


-- ISSUES WITH GRAPH OF MIXEDGRAPH
C=graph collateVertices LMG
vertices C
vertices G
vertices D
vertices B

graph C  -- does not work as it is but can be substituted by C#graph#Graph
digraph C
bigraph C

C#graph

graph C#graph

C#graph#Graph
C#graph#Digraph
C#graph#Bigraph

keys C
keys C#graph#Graph

keys graph C
     dd := graph G#Digraph;
     bb := graph G#Bigraph;
     uu := graph G#Graph;

keys graph LMG    

vv=vertices LMG
uu=graph C#Graph

U
kL = join(apply(U, i->k_(i,i)),delete(null, flatten apply(U, x-> apply(toList uu#x, y->if pos(U,x)>pos(U,y) then null else k_(x,y)))));
kL   
join(apply(U, i->k_(i,i)),delete(null, flatten apply(U, x-> apply(toList (graph gg#Graph)#x, y->if pos(U,x)>pos(U,y) then null else k_(x,y)))))

(graph G#Graph)#3  
(graph gg#Graph)#3

U

g

collateVertices g

apply(U, i->k_(i,i))
delete(null,flatten apply(U, x-> apply(toList (graph G#Graph)#x, y->if pos(U,x)>pos(U,y) then null else k_(x,y))))
flatten apply(U, x-> apply(toList (graph gg#Graph)#x, y->if pos(U,x)>pos(U,y) then null else k_(x,y)))

