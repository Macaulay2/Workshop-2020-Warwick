restart
needsPackage "TateOnProducts"
kk=ZZ/101
n={2}
(S,E)=productOfProjectiveSpaces(n,CoefficientField=>kk)
p1234=ideal(S_0*S_1 - S_0*S_2, S_0*S_2 - S_1*S_2) -- projective frame
degree p1234
p5=ideal random(S^1,S^{2:-1})
I=intersect(p1234,p5) -- ideal of a set of 5 general points
degree I
-- Resolution of I : 0->2O(-4)->2O(-3)+O(-2)->I->0
exactSeq = res I
betti exactSeq

T1=tateResolution(C=module I,{-8},{8})
T2=tateResolution(B=S^{1:-2,2:-3},{-8},{8})
T3=tateResolution(A=S^{2:-4},{-8},{8})
netList{betti T3, betti T2, betti T1}
-- they are mixed in lower degrees, but for a sufficiently large degree d,
-- it simply decomposes as a short exact sequence of vector spaces at deg d.
presentation C

-- the previous error was due to the fact that
-- (res I) does not resolve the S-module I but of S_Z; in particular,
-- (res I)_0 = S^1. 

-- Yeongrak: I think it is already implemented; if we have a map phi:A->B 
-- then truncate(d, phi) computes an induced map on their truncations at d.
-- Could you guys try to find why 'truncate(5,psi)' does not give an answer?
phi=exactSeq.dd_2; 
psi=map(C,B,sub(matrix{{1,0,0},{0,1,0},{0,0,1}},S))
isHomogeneous psi

truncate(5,phi) -- truncate(5,phi) looks okay; maps between free modules
truncate(5,psi) -- error!


-- complete this code to build two induced maps
-- A_d -> B_d, B_d -> C_d for d sufficiently large.
-- Lawrence's answer:
inclusion = inducedMap (target psi, target inducedMap (C, truncate(5,C)))
d1 = psi * inducedMap (B, truncate (5,B)) // inclusion // inducedMap (C, truncate(5,C))
d2 = phi * inducedMap (A, truncate (5,A)) // inducedMap (B, truncate(5,B))

source d2 == truncate(5,A)
target d2 == truncate(5,B)
source d1 == target d2
target d1 == truncate(5,C)

-- which looks too complicated.
------

