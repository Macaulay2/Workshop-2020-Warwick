restart
--installPackage "Graphs"
--installPackage "StatGraphs"
--installPackage "GraphicalModels"

loadPackage "Graphs"
loadPackage "StatGraphs"
loadPackage "GraphicalModels"
loadPackage "GraphicalModelsMLE"
debug needsPackage "GraphicalModelsMLE"
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

-- Test 3 (Roser)
restart
 debug needsPackage("GraphicalModelsMLE")
 G=graph{{1,2},{2,5},{5,6},{2,4},{4,5},{3,4}}
 g=mixedGraph G
 R2=gaussianRing(g)
 U=matrix{{1, 2, 9, 6, 0, 0}, {2, 7, 7, 3, 2, 2}, {6, 3, 4, 1, 5, 5}, {5, 5, 8, 8, 7, 6}, {3, 2, 3, 8, 7, 5}, {8, 0, 5, 3, 8, 5}}
 J2=scoreEquationsFromCovarianceMatrix(R2,U)
 
 L = directedEdgesMatrix R2
 K= undirectedEdgesMatrix R2 
 P= bidirectedEdgesMatrix R2
 
 d = numRows L
 (F,lpR)=changeRing(d,R2)
 FR = frac(lpR)
 
  -- Kinv
    K=substitute(K, FR)
    Kinv=inverse K
    P=substitute(P,FR)
       
     --Omega
    if K==0 then W=P else (if P==0 then W=Kinv else W = directSum(Kinv,P))
    --W:= directSum(Kinv,P);
    
    -- move to FR, the fraction field of lpR
    L= substitute (L,FR)
    
    -- Sigma
    if L==0 then S=W else (
	IdL = inverse (id_(FR^d)-L);
    	S = (transpose IdL) * W * IdL
	);
    if S == Kinv then Sinv= K else Sinv = inverse S;
    
    -- Sample covariance matrix
    V = sampleCovarianceMatrix(U)
     
    -- Compute ideal J   
    C1 = trace(Sinv * V)/2;
    C1derivative = JacobianMatrixOfRationalFunction(C1);
    LL =JacobianMatrixOfRationalFunction (det S)*matrix{{(-1/(2*det(S)))}} - (C1derivative);
    LL=flatten entries(LL);
    denoms = apply(#LL, i -> lift(denominator(LL_i), lpR));
    prod = product(denoms);
    J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR))
    J = saturate(J, prod)
    
    -- Solution
    -- Saturated ideal
    JSat= ideal(452*k_(5,6)+627, 28733*k_(4,5)+639, 1703*k_(3,4)+72, 28733*k_(2,5)-1781,
      28733*k_(2,4)+309, 524*k_(1,2)-51, 452*k_(6,6)-915, 3961131380*k_(5,5)-4311459839,
      12575600843*k_(4,4)-1906007136, 3406*k_(3,3)-771, 557075404*k_(2,2)-148839607,
      524*k_(1,1)-111)
  
  JNoSat= ideal(-209*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+209*k_(3,3)*k_(4,4)*k_(5,5)*k
      _(6,6)*k_(1,2)^2+209*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+209*k_(1,1)*k_(3,3)*k_(4,4)*
      k_(6,6)*k_(2,5)^2+209*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-209*k_(5,5)*k_(6,6)*k_(1,2
      )^2*k_(3,4)^2-209*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-418*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k
      _(2,5)*k_(4,5)+209*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-209*k_(3,3)*k_(6,6)*k_(1,2)^2*
      k_(4,5)^2+209*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-209*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,
      6)^2-209*k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-209*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2+209*k
      _(1,2)^2*k_(3,4)^2*k_(5,6)^2+36*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)-36*k_(3,3)*k_(5,5)*
      k_(6,6)*k_(2,4)^2-36*k_(3,3)*k_(4,4)*k_(6,6)*k_(2,5)^2-36*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2
      +36*k_(6,6)*k_(2,5)^2*k_(3,4)^2+72*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k_(4,5)-36*k_(2,2)*k_(3,
      3)*k_(6,6)*k_(4,5)^2-36*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2+36*k_(3,3)*k_(2,4)^2*k_(5,6)^2+36
      *k_(2,2)*k_(3,4)^2*k_(5,6)^2,
      -185*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+185*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k
      _(1,2)^2+185*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+185*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*
      k_(2,5)^2+185*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-185*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,
      4)^2-185*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-370*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k
      _(4,5)+185*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-185*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^
      2+185*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-185*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2-185
      *k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-185*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2+185*k_(1,2)^2*
      k_(3,4)^2*k_(5,6)^2+36*k_(1,1)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)-36*k_(1,1)*k_(5,5)*k_(6,6)*k
      _(3,4)^2-36*k_(1,1)*k_(3,3)*k_(6,6)*k_(4,5)^2-36*k_(1,1)*k_(3,3)*k_(4,4)*k_(5,6)^2+36*k_(1,
      1)*k_(3,4)^2*k_(5,6)^2, -14*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+14*k_(3,3)*k_(4
      ,4)*k_(5,5)*k_(6,6)*k_(1,2)^2+14*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+14*k_(1,1)*k_(3,
      3)*k_(4,4)*k_(6,6)*k_(2,5)^2+14*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-14*k_(5,5)*k_(6,6
      )*k_(1,2)^2*k_(3,4)^2-14*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-28*k_(1,1)*k_(3,3)*k_(6,6)*k_(
      2,4)*k_(2,5)*k_(4,5)+14*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-14*k_(3,3)*k_(6,6)*k_(1,2
      )^2*k_(4,5)^2+14*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-14*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(
      5,6)^2-14*k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-14*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2+14*k_(
      1,2)^2*k_(3,4)^2*k_(5,6)^2+3*k_(1,1)*k_(2,2)*k_(4,4)*k_(5,5)*k_(6,6)-3*k_(4,4)*k_(5,5)*k_(6
      ,6)*k_(1,2)^2-3*k_(1,1)*k_(5,5)*k_(6,6)*k_(2,4)^2-3*k_(1,1)*k_(4,4)*k_(6,6)*k_(2,5)^2+6*k_(
      1,1)*k_(6,6)*k_(2,4)*k_(2,5)*k_(4,5)-3*k_(1,1)*k_(2,2)*k_(6,6)*k_(4,5)^2+3*k_(6,6)*k_(1,2)^
      2*k_(4,5)^2-3*k_(1,1)*k_(2,2)*k_(4,4)*k_(5,6)^2+3*k_(4,4)*k_(1,2)^2*k_(5,6)^2+3*k_(1,1)*k_(
      2,4)^2*k_(5,6)^2, -257*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+257*k_(3,3)*k_(4,4)*
      k_(5,5)*k_(6,6)*k_(1,2)^2+257*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+257*k_(1,1)*k_(3,3
      )*k_(4,4)*k_(6,6)*k_(2,5)^2+257*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-257*k_(5,5)*k_(6,
      6)*k_(1,2)^2*k_(3,4)^2-257*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-514*k_(1,1)*k_(3,3)*k_(6,6)*
      k_(2,4)*k_(2,5)*k_(4,5)+257*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-257*k_(3,3)*k_(6,6)*k
      _(1,2)^2*k_(4,5)^2+257*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-257*k_(3,3)*k_(4,4)*k_(1,2
      )^2*k_(5,6)^2-257*k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-257*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6
      )^2+257*k_(1,2)^2*k_(3,4)^2*k_(5,6)^2+36*k_(1,1)*k_(2,2)*k_(3,3)*k_(5,5)*k_(6,6)-36*k_(3,3
      )*k_(5,5)*k_(6,6)*k_(1,2)^2-36*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,5)^2-36*k_(1,1)*k_(2,2)*k_(3,3
      )*k_(5,6)^2+36*k_(3,3)*k_(1,2)^2*k_(5,6)^2,
      -305*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+305*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k
      _(1,2)^2+305*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+305*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*
      k_(2,5)^2+305*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-305*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,
      4)^2-305*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-610*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k
      _(4,5)+305*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-305*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^
      2+305*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-305*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2-305
      *k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-305*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2+305*k_(1,2)^2*
      k_(3,4)^2*k_(5,6)^2+36*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(6,6)-36*k_(3,3)*k_(4,4)*k_(6,6)*k
      _(1,2)^2-36*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)^2-36*k_(1,1)*k_(2,2)*k_(6,6)*k_(3,4)^2+36*k_(6,
      6)*k_(1,2)^2*k_(3,4)^2, -161*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+161*k_(3,3)*k
      _(4,4)*k_(5,5)*k_(6,6)*k_(1,2)^2+161*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+161*k_(1,1)*
      k_(3,3)*k_(4,4)*k_(6,6)*k_(2,5)^2+161*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-161*k_(5,5
      )*k_(6,6)*k_(1,2)^2*k_(3,4)^2-161*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-322*k_(1,1)*k_(3,3)*k
      _(6,6)*k_(2,4)*k_(2,5)*k_(4,5)+161*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-161*k_(3,3)*k
      _(6,6)*k_(1,2)^2*k_(4,5)^2+161*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-161*k_(3,3)*k_(4,4
      )*k_(1,2)^2*k_(5,6)^2-161*k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-161*k_(1,1)*k_(2,2)*k_(3,4)^2
      *k_(5,6)^2+161*k_(1,2)^2*k_(3,4)^2*k_(5,6)^2+36*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)-36*
      k_(3,3)*k_(4,4)*k_(5,5)*k_(1,2)^2-36*k_(1,1)*k_(3,3)*k_(5,5)*k_(2,4)^2-36*k_(1,1)*k_(3,3)*k
      _(4,4)*k_(2,5)^2-36*k_(1,1)*k_(2,2)*k_(5,5)*k_(3,4)^2+36*k_(5,5)*k_(1,2)^2*k_(3,4)^2+36*k_(
      1,1)*k_(2,5)^2*k_(3,4)^2+72*k_(1,1)*k_(3,3)*k_(2,4)*k_(2,5)*k_(4,5)-36*k_(1,1)*k_(2,2)*k_(3
      ,3)*k_(4,5)^2+36*k_(3,3)*k_(1,2)^2*k_(4,5)^2,
      85*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)-85*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k_(1,
      2)^2-85*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2-85*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*k_(2,5
      )^2-85*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2+85*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,4)^2+85*
      k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2+170*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k_(4,5)-85*
      k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2+85*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^2-85*k_(1,1
      )*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2+85*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2+85*k_(1,1)*k_(3,
      3)*k_(2,4)^2*k_(5,6)^2+85*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2-85*k_(1,2)^2*k_(3,4)^2*k_(5,6
      )^2-36*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k_(1,2)+36*k_(5,5)*k_(6,6)*k_(1,2)*k_(3,4)^2+36*k_(3
      ,3)*k_(6,6)*k_(1,2)*k_(4,5)^2+36*k_(3,3)*k_(4,4)*k_(1,2)*k_(5,6)^2-36*k_(1,2)*k_(3,4)^2*k_(
      5,6)^2, -k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k
      _(1,2)^2+k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*k_(2,5)^
      2+k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,4)^2-k_(1,1)*k_(
      6,6)*k_(2,5)^2*k_(3,4)^2-2*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k_(4,5)+k_(1,1)*k_(2,2)*
      k_(3,3)*k_(6,6)*k_(4,5)^2-k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^2+k_(1,1)*k_(2,2)*k_(3,3)*k_(4,
      4)*k_(5,6)^2-k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2-k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-k_(1,1
      )*k_(2,2)*k_(3,4)^2*k_(5,6)^2+k_(1,2)^2*k_(3,4)^2*k_(5,6)^2-36*k_(1,1)*k_(3,3)*k_(5,5)*k_(6
      ,6)*k_(2,4)+36*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,5)*k_(4,5)+36*k_(1,1)*k_(3,3)*k_(2,4)*k_(5,6)^2
      , 83*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)-83*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k_(
      1,2)^2-83*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2-83*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*k_(2
      ,5)^2-83*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2+83*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,4)^2+
      83*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2+166*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k_(4,5)-
      83*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2+83*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^2-83*k_(1
      ,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2+83*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2+83*k_(1,1)*k_(
      3,3)*k_(2,4)^2*k_(5,6)^2+83*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2-83*k_(1,2)^2*k_(3,4)^2*k_(5
      ,6)^2-36*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*k_(2,5)+36*k_(1,1)*k_(6,6)*k_(2,5)*k_(3,4)^2+36*k
      _(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(4,5),
      -4*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+4*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k_(1,2
      )^2+4*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+4*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*k_(2,5)^2
      +4*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-4*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,4)^2-4*k_(1,1
      )*k_(6,6)*k_(2,5)^2*k_(3,4)^2-8*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k_(4,5)+4*k_(1,1)*k
      _(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-4*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^2+4*k_(1,1)*k_(2,2)*k
      _(3,3)*k_(4,4)*k_(5,6)^2-4*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2-4*k_(1,1)*k_(3,3)*k_(2,4)^2*
      k_(5,6)^2-4*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2+4*k_(1,2)^2*k_(3,4)^2*k_(5,6)^2-3*k_(1,1)*k
      _(2,2)*k_(5,5)*k_(6,6)*k_(3,4)+3*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,4)+3*k_(1,1)*k_(6,6)*k_(2,5
      )^2*k_(3,4)+3*k_(1,1)*k_(2,2)*k_(3,4)*k_(5,6)^2-3*k_(1,2)^2*k_(3,4)*k_(5,6)^2,
      -41*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+41*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k_(1
      ,2)^2+41*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+41*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*k_(2,
      5)^2+41*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-41*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,4)^2-41
      *k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-82*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k_(4,5)+41*
      k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-41*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^2+41*k_(1,1
      )*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-41*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2-41*k_(1,1)*k_(3,
      3)*k_(2,4)^2*k_(5,6)^2-41*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2+41*k_(1,2)^2*k_(3,4)^2*k_(5,6
      )^2+36*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)-36*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)+
      36*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5),
      -209*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)+209*k_(3,3)*k_(4,4)*k_(5,5)*k_(6,6)*k
      _(1,2)^2+209*k_(1,1)*k_(3,3)*k_(5,5)*k_(6,6)*k_(2,4)^2+209*k_(1,1)*k_(3,3)*k_(4,4)*k_(6,6)*
      k_(2,5)^2+209*k_(1,1)*k_(2,2)*k_(5,5)*k_(6,6)*k_(3,4)^2-209*k_(5,5)*k_(6,6)*k_(1,2)^2*k_(3,
      4)^2-209*k_(1,1)*k_(6,6)*k_(2,5)^2*k_(3,4)^2-418*k_(1,1)*k_(3,3)*k_(6,6)*k_(2,4)*k_(2,5)*k
      _(4,5)+209*k_(1,1)*k_(2,2)*k_(3,3)*k_(6,6)*k_(4,5)^2-209*k_(3,3)*k_(6,6)*k_(1,2)^2*k_(4,5)^
      2+209*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)^2-209*k_(3,3)*k_(4,4)*k_(1,2)^2*k_(5,6)^2-209
      *k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)^2-209*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)^2+209*k_(1,2)^2*
      k_(3,4)^2*k_(5,6)^2-36*k_(1,1)*k_(2,2)*k_(3,3)*k_(4,4)*k_(5,6)+36*k_(3,3)*k_(4,4)*k_(1,2)^2
      *k_(5,6)+36*k_(1,1)*k_(3,3)*k_(2,4)^2*k_(5,6)+36*k_(1,1)*k_(2,2)*k_(3,4)^2*k_(5,6)-36*k_(1,
      2)^2*k_(3,4)^2*k_(5,6))
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

restart

myFunction = method(Options =>  options saturate ++ {ABC => true} );
options myFunction
myFunction (Ideal,Ideal):= opts -> (I,J)->(
    if opts.ABC then (saturate (I,J, opts)) 
    else return I;
    ); 
R = ZZ/32003[a..d];
I = ideal(a^3-b, a^4-c);
Ih = homogenize(I,d);
saturate(Ih,d)
J = ideal d
myFunction(Ih,J)
myFunction(Ih,J, ABC=> true)
myFunction(Ih,J, ABC=> false)
myFunction(Ih,d, ABC=> true, DegreeLimit=>1)

restart
myFunction1 = method(Options =>  options saturate );

myFunction1 (Ideal,Ideal):= opts -> (I,J)->(
   saturate (I,J, opts)
    ); 

R = ZZ/32003[a..d];
I = ideal(a^3-b, a^4-c);
Ih = homogenize(I,d);
saturate(Ih,d)
J = ideal d
myFunction1 (Ih,J)
myFunction1 (Ih,J)===saturate(Ih,d)

myFunction1(Ih,J, DegreeLimit=>1)

restart
myFunction2 = method(Options =>  {doSaturate => true, saturateOptions => options saturate});
myFunction2 (Ideal,Ideal):= opts -> (I,J)->(
    if opts.doSaturate then (
	g:=opts.saturateOptions  >>opts-> args ->(args, opts);
	saturate (g(I,J)))
    else return I
    ); 

R = ZZ/32003[a..d];
I = ideal(a^3-b, a^4-c);
Ih = homogenize(I,d);
saturate(Ih,d)
J = ideal d
myFunction2 (Ih,J)
myFunction2(Ih,J)===saturate(Ih,d)
myFunction2(Ih,J, doSaturate=> true)
myFunction2(Ih,J, doSaturate=> true)===saturate(Ih,d)
myFunction2(Ih,J, doSaturate=> false)
myFunction2(Ih,J, doSaturate=> true, saturateOptions => {DegreeLimit=>1})
myFunction2(Ih,J, doSaturate=> true, saturateOptions => {DegreeLimit=>1})===saturate(Ih,d, DegreeLimit=>1)
myFunction2(Ih,J, doSaturate=> true, saturateOptions => {DegreeLimit=>1, MinimalGenerators => false})===saturate(Ih,d, DegreeLimit=>1,MinimalGenerators => false )

g=options saturate >>opts-> args ->(args, opts)

g={DegreeLimit=>1,MinimalGenerators => false}>>opts-> args ->(args, opts)
g x
g (I,J)
saturate (g (I,J))
saturate (I,J)===saturate (g (I,J))
g (I,J,DegreeLimit=>1)
g (I,J,DegreeLimit=>1,MinimalGenerators => false)
g (I,J,{DegreeLimit=>1,MinimalGenerators => false})

g(I,J, {DegreeLimit=>1})
saturate (I,J,DegreeLimit=>1)===saturate (g (I,J,DegreeLimit=>1))

g(Ih,J)


myFunction2 = method(Options =>  {doSaturate => true, saturateOptions => options saturate});
myFunction2 (Ideal,Ideal):= opts -> (I,J)->(
    if opts.doSaturate then (saturate (I,J, opts.saturateOptions)) 
    else return I
    ); 