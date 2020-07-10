--************************************************************************
-- Relevant methods from GaussianModels.m2    --
--************************************************************************

gaussianRingData = local gaussianRingData
gaussianVariables = local gaussianVariables

----------------------------------------------
-- gaussianRing MixedGraph
----------------------------------------------

gaussianRing MixedGraph := Ring => opts -> (g) -> (
     G := graph collateVertices g;
     dd := graph G#Digraph;
     bb := graph G#Bigraph;
     uu := G#Graph;
     if #(edges uu) > 0 then error "mixedgraph must have no undirected part ";
     vv := sort vertices g;
     s := toSymbol opts.sVariableName;
     l := toSymbol opts.lVariableName;
     p := toSymbol opts.pVariableName;
     kk := opts.Coefficients;          
     if (not gaussianRingList#?(kk,s,l,p,vv)) then ( 
	  --(kk,s,l,p,vv) uniquely identifies gaussianRing in case of MixedGraph input.
     sL := delete(null, flatten apply(vv, x-> apply(vv, y->if pos(vv,x)>pos(vv,y) then null else s_(x,y))));
     lL := delete(null, flatten apply(vv, x-> apply(toList dd#x, y->l_(x,y))));	 
     pL := join(apply(vv, i->p_(i,i)),delete(null, flatten apply(vv, x-> apply(toList bb#x, y->if pos(vv,x)>pos(vv,y) then null else p_(x,y)))));
     m := #lL+#pL;
     R := kk(monoid [lL,pL,sL,MonomialOrder => Eliminate m, MonomialSize=>16]);
     -- create gaussianVariables hash table: (symbol s)_(i,j) => ring var with the same name, same for l, p.
     H := new MutableHashTable;
     nextvar := 0;
     for v in lL do (H#v = R_nextvar; nextvar = nextvar+1);
     for v in pL do (H#v = R_nextvar; nextvar = nextvar+1);
     for v in sL do (H#v = R_nextvar; nextvar = nextvar+1);
     R.gaussianVariables = new HashTable from H;
     R#numberOfEliminationVariables = m;
     R.gaussianRingData = {#vv,s,l,p};
     R.mixedGraph = g;
     gaussianRingList#((kk,s,l,p,vv)) = R;); 
     gaussianRingList#((kk,s,l,p,vv))
     )
 
 ------------------------------------------------------------------
-- undirectedEdgesMatrix Ring 
------------------------------------------------------------------

undirectedEdgesMatrix = method()
undirectedEdgesMatrix Ring := Matrix =>  R -> (
     if not (R.?graph and R.?gaussianRingData) then error "expected a ring created with gaussianRing of a Graph";
     g := R.graph;
     bb:= graph g;
     vv := sort vertices g;
     n := R.gaussianRingData#0; --number of vertices
     k := R.gaussianRingData#2; 
     H := R.gaussianVariables;
     PM := mutableMatrix(R,n,n);
     scan(vv,i->PM_(pos(vv,i),pos(vv,i))=H#(k_(i,i)));
     scan(vv,i->scan(toList bb#i, j->PM_(pos(vv,i),pos(vv,j))=if pos(vv,i)<pos(vv,j) then H#(k_(i,j)) else H#(k_(j,i))));
     matrix PM) 



------------------------------------------------------------------
-- directedEdgesMatrix Ring 
------------------------------------------------------------------

directedEdgesMatrix = method()
directedEdgesMatrix Ring := Matrix => R -> (
     if not (R.?mixedGraph and R.?gaussianRingData) then error "expected a ring created with gaussianRing of a MixedGraph";     
     g := R.mixedGraph;
     G := graph collateVertices g;
     dd := graph G#Digraph;
     vv := sort vertices g;
     n := R.gaussianRingData#0;
     l := R.gaussianRingData#2;
     H := R.gaussianVariables;
     LM := mutableMatrix(R,n,n);
     scan(vv,i->scan(toList dd#i, j->LM_(pos(vv,i),pos(vv,j))=H#(l_(i,j))));
     matrix LM) 
 
 ------------------------------------------------------------------
-- covarianceMatrix(Ring)
------------------------------------------------------------------

covarianceMatrix = method()
covarianceMatrix(Ring) := Matrix => (R) -> (
       if not R.?gaussianRingData then error "expected a ring created with gaussianRing";    
       if R.?graph then (  
     	    g:=R.graph;
	    vv := sort vertices g;
     	    n := R.gaussianRingData#0;
     	    s := R.gaussianRingData#1;
            H := R.gaussianVariables;
     	    SM := mutableMatrix(R,n,n);
     	    scan(vv,i->scan(vv, j->SM_(pos(vv,i),pos(vv,j))=if pos(vv,i)<pos(vv,j) then H#(s_(i,j)) else H#(s_(j,i))));
     	    matrix SM	    
	    ) 
       else if R.?mixedGraph then (  
     	    g = R.mixedGraph;
	    vv = sort vertices g;
     	    n = R.gaussianRingData#0;
     	    s = R.gaussianRingData#1;
            H = R.gaussianVariables;
     	    SM = mutableMatrix(R,n,n);
     	    scan(vv,i->scan(vv, j->SM_(pos(vv,i),pos(vv,j))=if pos(vv,i)<pos(vv,j) then H#(s_(i,j)) else H#(s_(j,i))));
     	    matrix SM	    
	    ) 
       else (
	    n =R.gaussianRingData; 
	    genericSymmetricMatrix(R,n)
	    )
  )