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
U = matrix({{X_0},{X_1},{X_2}})

S=sampleCovarianceMatrix(X)   
S1 = sampleCovarianceMatrix(U)   -- checking the new written code for matrices

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
