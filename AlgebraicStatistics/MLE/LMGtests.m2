restart
--installPackage "Graphs"
--installPackage "StatGraphs"
--installPackage "GraphicalModels"

loadPackage "Graphs"
loadPackage "StatGraphs"
loadPackage "GraphicalModels"
debug loadPackage "GraphicalModelsMLE"
--debug needsPackage "GraphicalModelsMLE"
--------------------------
-- Testing LMG function
--------------------------
-- Test 1: Input
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

-- Test 1: Confirm the function gives correct output
J=scoreEquationsFromCovarianceMatrix(R,U);
I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J===I)

-- Test 1: Confirm no saturation works
JnoSat=scoreEquationsFromCovarianceMatrix(R,U,Saturate=>false);
dim JnoSat
degree JnoSat
   
-- Test 2: Input (best to restart and reload packages first)
restart

loadPackage "Graphs"
loadPackage "StatGraphs"
loadPackage "GraphicalModels"
debug loadPackage "GraphicalModelsMLE"

U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
G=graph{{1,2},{2,3},{3,4},{1,4}}
G=mixedGraph G
R=gaussianRing(G)

-- Test 2: Run with the LMG function
J=scoreEquationsFromCovarianceMatrix(R, U)
assert(dim J===0)
assert(degree J===5)
test= ideal(58*k_(4,4)+47*k_(1,4)-21*k_(3,4)-8,
      27*k_(3,3)+14*k_(2,3)-42*k_(3,4)-16,
      10*k_(2,2)-13*k_(1,2)+7*k_(2,3)-8,
      115*k_(1,1)-26*k_(1,2)+94*k_(1,4)-16,
      49587264*k_(3,4)^2+159578097*k_(1,2)-401770875*k_(1,4)+86063425*k_(2,3
      )-241038399*k_(3,4)-9279488,
      289984*k_(2,3)*k_(3,4)-2996077*k_(1,2)+3236687*k_(1,4)-1267133*k_(2,3
      )+2001475*k_(3,4), 289984*k_(1,4)*k_(3,4)+572663*k_(1,2)-1267133*k_(1,
      4)+247207*k_(2,3)-713625*k_(3,4),
      289984*k_(1,2)*k_(3,4)-1267133*k_(1,2)+1588223*k_(1,4)-634637*k_(2,3)+
      786099*k_(3,4), 12469312*k_(2,3)^2+159578097*k_(1,2)-401770875*k_(1,4
      )+94182977*k_(2,3)-216679743*k_(3,4)-9279488,
      289984*k_(1,4)*k_(2,3)-1428105*k_(1,2)+2001475*k_(1,4)-713625*k_(2,3)+
      983079*k_(3,4), 289984*k_(1,2)*k_(2,3)+2001475*k_(1,2)-3960705*k_(1,4
      )+786099*k_(2,3)-2523789*k_(3,4),
      163260992*k_(1,4)^2+159578097*k_(1,2)-347253883*k_(1,4)+86063425*k_(2,
      3)-216679743*k_(3,4)-9279488,
      289984*k_(1,2)*k_(1,4)-713625*k_(1,2)+786099*k_(1,4)-302505*k_(2,3)+
      482391*k_(3,4), 58866752*k_(1,2)^2+144498929*k_(1,2)-401770875*k_(1,4
      )+86063425*k_(2,3)-216679743*k_(3,4)-9279488);
assert(J===test)

------------------------------------------------
-- Tests for gaussianRing, ...EdgesMatrix
------------------------------------------------
G=graph{{1,2},{1,3},{2,3}}
D=digraph{{1,6},{4,7}}
B=bigraph{{5,6},{6,7}}

G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}

-- MixedGraph with all components
g=mixedGraph(G,D,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without B
g=mixedGraph(G,D)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without D
g=mixedGraph(G,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without G
g=mixedGraph(D,B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without D,B
g=mixedGraph(G)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without G,B
g=mixedGraph(D)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- MixedGraph without G,D
g=mixedGraph(B)
R=gaussianRing g

undirectedEdgesMatrix R
directedEdgesMatrix R
bidirectedEdgesMatrix R
covarianceMatrix R

-- random function with options
r = method(Options => {Test => true});
r RR := opts -> x -> (
    
    if opts.Test then x^2 else x+1
    );
r(0.5)
r(0.5, Test => true)
r(0.5, Test => false)

r = method(Options => {Test => true});
r RR := opts -> x -> (
    y:=x;
    if opts.Test then y=x^2;
    return y
    );
r(0.5)
r(0.5, Test => false)