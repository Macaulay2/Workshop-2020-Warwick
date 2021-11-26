restart;

uninstallPackage("GraphicalModelsMLE")
installPackage("GraphicalModelsMLE")
check GraphicalModelsMLE
help GraphicalModelsMLE

----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------

--DOCUMENTATION REVIEW

--checkPD, checkPSD
restart
loadPackage "GraphicalModelsMLE"

help checkPD
help checkPSD


--Remove roundMatrix, use numeric
restart
loadPackage "GraphicalModelsMLE"
help scoreEquations
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph{{3,4}});
R = gaussianRing(G);
U = random(RR^4,RR^4)
J = scoreEquations(R,U)

V=sampleCovarianceMatrix(U)

n=5
numeric(n,V_(0,0))
numeric(n,V_(2,3))
numeric(n,V_(2,3))^QQ
lift(numeric(n,V_(2,3)),QQ)
lift(numeric(n,V_(0,0)),QQ)

VV=matrix apply(entries V,r->r/(v->lift(numeric(53,v),QQ)))
scoreEquations(R,VV,SampleData=>false,RealPrecision=>5)

--Output of solver MLE: sequence documentation, type of first entry(RR instead of old CC)

restart
uninstallPackage "GraphicalModelsMLE"
installPackage "GraphicalModelsMLE"
help solverMLE

G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph{{3,4}});
R = gaussianRing(G);
U = random(RR^4,RR^4)
(mx,MLE,ml)=solverMLE(G,U)
mx


--Issues with documentation of RealPrecision
restart
uninstallPackage "GraphicalModelsMLE"
installPackage "GraphicalModelsMLE"

help [solverMLE,RealPrecision]
help [scoreEquations,RealPrecision]

--Similar issues with documentation of optional inputs
restart
uninstallPackage "GraphicalModelsMLE"
installPackage "GraphicalModelsMLE"

help [scoreEquations,DoSaturate]
help [solverMLE,ChooseSolver]
