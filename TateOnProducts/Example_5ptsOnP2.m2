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
betti (resI=res I)
-- Resolution of I : 0->2O(-4)->2O(-3)+O(-2)->I->0

T1=tateResolution(C=module I,{-8},{8})
T2=tateResolution(B=S^{2:-3,1:-2},{-8},{8})
T3=tateResolution(A=S^{2:-4},{-8},{8})
netList{betti T3, betti T2, betti T1}
-- they are mixed in lower degrees, but for a sufficiently large degree d,
-- it simply decomposes as a short exact sequence of vector spaces at deg d.

phi=resI.dd_2; psi=resI.dd_1;
matrix truncate(5,phi); -- truncate(5,phi) looks okay; maps between free modules
matrix truncate(5,psi)

-- complete this code to build two induced maps
-- A_d -> B_d, B_d -> C_d for d sufficiently large.
