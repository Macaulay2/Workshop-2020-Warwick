restart
--installPackage "Graphs"
--installPackage "StatGraphs"
--installPackage "GraphicalModels"

loadPackage "Graphs"
loadPackage "StatGraphs"
loadPackage "GraphicalModels"
loadPackage "GraphicalModelsMLE"
debug needsPackage "GraphicalModelsMLE"
--loadPackage "GraphicalModelsMLEorig"
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
   
-- Test 2: Input
U=random(ZZ^4,ZZ^4)
G1=graph{{1,2},{2,3},{3,4},{1,4}}
G2=mixedGraph G1
R1=gaussianRing(G1)
R2=gaussianRing(G2)

-- Test2: Run with the code for undirected edges only
J=scoreEquationsFromCovarianceMatrixUndir(R1,U)
assert(dim J===0)
assert(degree J===5)

-- Test 2: Run with the LMG function
J_LMG=scoreEquationsFromCovarianceMatrix(R2,U);
J===J_LMG


-- Test 2: Debug code
L = directedEdgesMatrix R2
K = undirectedEdgesMatrix R2
P = bidirectedEdgesMatrix R2

d = numRows L
(F,lpR)=changeRing(d,R2)
FR = frac(lpR);
 
Kinv=inverse substitute(K, FR)
P=substitute(P,FR)
W= directSum(Kinv,P)

L= substitute (L,FR)
IdL = inverse (id_(FR^d)-L)
S = (transpose IdL) * W * IdL
--describe S
Sinv= inverse S

V = sampleCovarianceMatrix(U)
C1 = trace(Sinv * V)/2;
C1derivative = JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
LL = (substitute(jacobian( substitute( matrix {{det S}},lpR)),FR))*matrix{{(-1/(2*det(S)))}} - (C1derivative);
LL=flatten entries(LL)
denoms = apply(#LL, i -> lift(denominator(LL_i), lpR));
prod = product(denoms);
J=ideal apply(#LL, i -> lift(numerator(LL_i),lpR));
J = saturate(J, prod);
J

JUndir=scoreEquationsFromCovarianceMatrixUndir(R1,U)
J===JUndir
J
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


