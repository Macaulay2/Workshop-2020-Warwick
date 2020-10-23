restart

installPackage "GraphicalModels"
debug loadPackage "GraphicalModels"

--Test pairMarkov
G = graph({{a,b},{b,c},{c,d},{d,e},{e,a}})
pairMarkov G

D = digraph {{1,{2,3}}, {2,{4,5}}, {3,{5,6}}, {4,{7}}, {5,{7}},{6,{7}},{7,{}}}
pairMarkov D

D1=digraph{{1,2},{2,1}}
pairMarkov D1

-- Test localMarkov
localMarkov G
localMarkov D
localMarkov D1

-- Test globalMarkov
globalMarkov G
globalMarkov D
globalMarkov D1

-- test markovRing
d=(2,3,4,5);
R = markovRing d
numgens R
R_0, R_1, R_119
class R

-- test gaussianRing
G=graph{{1,2},{1,3},{2,3}}
D=digraph{{1,6},{4,7}}
B=bigraph{{5,6},{6,7}}
M=mixedGraph(G,B)

R1=gaussianRing G
class R1

R2=gaussianRing D
class R2

R3=gaussianRing B
class R3

R4=gaussianRing M
class R4

-- test undirectedEdgesMatrix

G=graph{{1,2},{1,3},{2,3}}
D=digraph{{1,6},{4,7}}
B=bigraph{{5,6},{6,7}}

M1=mixedGraph(G,B)
M2=mixedGraph(G)
M3=mixedGraph(D,B)

R1=gaussianRing M1
R2=gaussianRing M2
R3=gaussianRing M3
R4=gaussianRing D
R5=gaussianRing G

undirectedEdgesMatrix R1
undirectedEdgesMatrix R2
undirectedEdgesMatrix R5
undirectedEdgesMatrix R3
undirectedEdgesMatrix R4

-- test directedEdgesMatrix

G=graph{{1,2},{1,3},{2,3}}
D=digraph{{1,6},{4,7}}
B=bigraph{{5,6},{6,7}}

M1=mixedGraph(D,B)
M2=mixedGraph(D)
M3=mixedGraph(G,B)

R1=gaussianRing M1
R2=gaussianRing M2
R3=gaussianRing M3
R4=gaussianRing D
R5=gaussianRing G

directedEdgesMatrix R1
directedEdgesMatrix R2
directedEdgesMatrix R5
directedEdgesMatrix R3
directedEdgesMatrix R4

-- test covarianceMatrix

R=gaussianRing 5
covarianceMatrix R

covarianceMatrix R2
covarianceMatrix R4
covarianceMatrix gaussianRing B
covarianceMatrix gaussianRing mixedGraph B
covarianceMatrix gaussianRing G

-- test gaussianMatrices
Stmts = {{{1,2},{3},{4}}, {{1},{3},{}}}
gaussianMatrices (R, Stmts)
gaussianMatrices (R2, Stmts) -- testing error message

Stmts = {{{1,6},{4},{7}}, {{6},{7},{}}}
gaussianMatrices (gaussianRing D, Stmts)
gaussianMatrices (gaussianRing mixedGraph D, Stmts)

Stmts = {{{1,2},{3},{}}, {{1},{2},{}}}
gaussianMatrices(gaussianRing G, Stmts)
gaussianMatrices(gaussianRing mixedGraph G, Stmts)

Stmts = {{{5,6},{7},{}}, {{5},{6},{}}}
gaussianMatrices(gaussianRing B, Stmts)
gaussianMatrices(gaussianRing mixedGraph B, Stmts)

-- test conditionalIndependenceIdeal
GM = graph {{1,2},{2,3},{3,4},{4,1}}
DM = digraph {{1,{}},{2,{1}},{3,{1}},{4,{2,3}}}
RM = markovRing (2,2,2,2)
conditionalIndependenceIdeal (RM, globalMarkov(GM))
conditionalIndependenceIdeal (RM, globalMarkov(DM))

Stmts = {{{1,2},{3},{}}, {{1},{2},{}}}
conditionalIndependenceIdeal(R5, Stmts)
conditionalIndependenceIdeal(gaussianRing mixedGraph G, Stmts)

Stmts = {{{1,6},{4},{7}}, {{6},{7},{}}}
conditionalIndependenceIdeal (gaussianRing D, Stmts)
conditionalIndependenceIdeal (gaussianRing mixedGraph D, Stmts)

-- test gaussianParametrization
gaussianParametrization R1
gaussianParametrization R2
gaussianParametrization R3
gaussianParametrization gaussianRing mixedGraph D
gaussianParametrization gaussianRing D

-- test gaussianVanishingIdeal Ring
gaussianVanishingIdeal QQ[x]
gaussianVanishingIdeal gaussianRing 3
gaussianVanishingIdeal gaussianRing D
gaussianVanishingIdeal gaussianRing G
gaussianParametrization R1
gaussianParametrization R2
gaussianParametrization R3
gaussianParametrization gaussianRing mixedGraph D

-- test discreteVanishingIdeal
discreteVanishingIdeal (QQ[x],D)
d=(2,3,4,5);
R = markovRing d
discreteVanishingIdeal (R,D)

restart
debug loadPackage "GraphicalModels"
D1 = digraph {{1,{}}, {2,{}}}
Rnew = markovRing (2,2)
discreteVanishingIdeal (Rnew,D1)
discreteVanishingIdeal (Rnew,G)
-- test trekSeparation

GTS = mixedGraph(digraph {{b,{c,d}},{c,{d}}},bigraph {{a,d}})
S = trekSeparation GTS

GTS = mixedGraph(digraph {{b,{c,d}},{c,{d}}})
S = trekSeparation GTS

GTS = mixedGraph(bigraph {{a,d}})
S = trekSeparation GTS
-- test trekIdeal
G= mixedGraph(digraph {{b,{c,d}},{c,{d}}},bigraph {{a,d}})
R=gaussianRing GTI
T=trekIdeal(RTI,GTI)

g = graph{{a,b},{b,c},{c,d},{a,d}}
R = gaussianRing g
T = trekIdeal(R,g)

g = mixedGraph(graph{{a,b},{b,c},{c,d},{a,d}})
R = gaussianRing g
T = trekIdeal(R,g)

g = digraph{{1,{4}},{2,{4}},{3,{4,5}},{4,{5}}}
R = gaussianRing g
T = trekIdeal(R,g)

g=mixedGraph(g)
R = gaussianRing g
T = trekIdeal(R,g)

-- test marginMap
d=(2,3,4,5);
R = markovRing d
marginMap(1,R)
marginMap(5,R)

-- test inverseMarginMap
d=(2,3,4,5);
R = markovRing d
inverseMarginMap(1,R)
inverseMarginMap(5,R)

-- test hiddenMap
hiddenMap(1,R)
hiddenMap(5,R)

-- test identifyParameters
G = mixedGraph(digraph {{a,{b}},{b,{c}}},bigraph {{a,c}, {b,c}})
R = gaussianRing G
H = identifyParameters R
