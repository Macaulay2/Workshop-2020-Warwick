restart
loadPackage "GraphicalModelsMLE"
G=digraph{{1,2},{2,3}}
Gmixed=mixedGraph(G,bigraph{})
R=gaussianRing Gmixed
describe R
L=directedEdgesMatrix R
W=bidirectedEdgesMatrix R
covarianceMatrix R

--produce sample covariance matrix
d=numrows L
X={}
for i from 1 to d do (X=append(X,random(ZZ^1,ZZ^d)))
X

L = for i from 0 to d-1 list {X_i};
L
U1 = matrix(L)
U = matrix({{X_0},{X_1},{X_2}})

S=sampleCovarianceMatrix(X)   
S1 = sampleCovarianceMatrix(U)   -- checking the new written code for matrices


changeRing=(d,R)->(
    -- count the number of S variables. d is the number of
    -- rows in a relevant matrix (clarify!!)
    numSvars:=lift(d*(d+1)/2,ZZ);
    --lp ring is the ring without the s variables
    lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
    KK:=coefficientRing(R);
    lpR:=KK[lpRvarlist];
    -- here i is taken to numgens(R)-numSvars-1 because
    -- indexing starts from 0. But for subscripts it
    -- starts from 1.
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    return (F,lpR)
    );


removeSvar = (R) ->(
    numSvars := #set support covarianceMatrix R;
    lpRvarlist := gens R - set support covarianceMatrix R;
    KK := coefficientRing(R);
    lpR := KK[lpRvarlist];
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    return {F,lpR};
);
G=graph{{1,2},{2,3},{3,4},{1,4}}
R=gaussianRing(G)
 K:=undirectedEdgesMatrix R

F = removeSvar(R)

F_0

F_1



 

K = F_0(K)



























-- change rings
R2=coefficientRing(R)[select(gens R,v-> first baseName v!=symbol s)]
describe R2
R2target=gens R2
for i from 1 to numgens R-numgens R2 do R2target=append(R2target,0)
R2target
F=map(R2,R,R2target)
  
--compute MLE
L = F(L)
W = F(W) 
Sigma = (transpose inverse (id_(R2^d)-L)) * W * inverse (id_(R2^d)-L)
FR = frac(R2)
describe FR
K = inverse sub(Sigma,FR)

-- describe the partial derivatives of -n/2 log det Sigma - n/2 tr (S*K) w.r.t the variables of R2 in the field of fractions FR 
LL = (sub((jacobian(matrix{{det(Sigma)}})), FR))*matrix{{(-1/(2*det(Sigma)))}} - JacobianMatrixOfRationalFunction(trace(K * S)/2)
-- lift the equations from FR to R2
LL=flatten entries(LL)
denoms = apply(#LL, i -> lift(denominator(LL_i), R2))
prod = product(denoms)
J=ideal apply(#LL, i -> lift(numerator(LL_i),R2))
J = saturate(J, prod)

dim J, degree J
sols=zeroDimSolve(J)

--convert list of points to list of lists
T={}
for l in sols do (T = T|{coordinates(l)})	

-- substitute solutions on K and Sigma
MK = {};
for t in T do ( m = substitute(K,matrix{t}),MK = MK|{m});    
MSigma = {};
for t in T do ( ms = substitute(Sigma,matrix{t}),MSigma = MSigma|{ms});    
-- check one is the inverse of the other (true up to numerical error)
MK_0
inverse MSigma_0
-- check positive definiteness
eigenvalues MSigma_0




restart
needsPackage "GraphicalModelsMLE"
G = mixedGraph(digraph {{1,3},{2,4}},bigraph {{3,4}})
S = matrix{{2.93, -1.7, 0.76, -0.06},
{-1.7, 1.64, -0.78, 0.1},
{0.76, -0.78, 1.66, -0.78},
{-0.06, 0.1, -0.78, 0.81}}
solverMLE(G,S)



roundMatrix = method() -- only accepts real matrices
roundMatrix (ZZ, Matrix) := Matrix => (n, A) -> matrix apply(entries A, r -> r/(e -> (round(n,0.0+e))^QQ))

S
roundMatrix(53,S)

M = matrix{{sqrt(2),sqrt(3)}}

roundMatrix(51,M)
for i from 1 to 100 do 
(
    print i;
    print roundMatrix(i,M);
)

roundList = method(); 
roundList (ZZ, List) := List => (n, L) -> apply(L, r -> r/(e -> (round(n,0.0+e))^QQ))

roundList(5,S)

n = numrows S
X = for i to n-1 list S^{i};
X


X = for i to n-1 list flatten entries S^{i};
X

-- I think it should be 









--checks for installations

restart
uninstallPackage "Graphs"
uninstallPackage "StatGraphs"
uninstallPackage "EigenSolver"
uninstallPackage "GraphicalModels"
uninstallPackage "GraphicalModelsMLE"


restart
installPackage "Graphs"
installPackage "StatGraphs"
installPackage "EigenSolver"
installPackage "GraphicalModels"
installPackage "GraphicalModelsMLE"

help StatGraphs
