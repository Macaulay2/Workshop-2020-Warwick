restart
--installPackage "Graphs"
--installPackage "StatGraphs"
--installPackage "GraphicalModels"

loadPackage "Graphs"
loadPackage "StatGraphs"
loadPackage "GraphicalModels"
loadPackage "GraphicalModelsMLE"

debug needsPackage "GraphicalModelsMLE"

---------------
-- Testing noSaturate
---------------
restart
R = ZZ/32003[a..d];
I = ideal(a^3-b, a^4-c)
Ih = homogenize(I,d)

saturate(Ih,d)

J={1,3}
saturate (I,J)

J={a+1-a,b+3-b}

saturate (I, ideal J)

degree J_0


-----------------------
--Testing MLEmax function
-----------------------
restart
debug loadPackage "GraphicalModelsMLE"
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
V= sampleCovarianceMatrix U
--L={random(QQ^4,QQ^4), random(QQ^4,QQ^4), random(QQ^4,QQ^4), random(QQ^4,QQ^4)}
L={map(QQ^4,QQ^4,{{8/7, 5, 3, 7/5}, {9, 5/2, 1, 2/9}, {2, 2/5, 10/9,
      7/5}, {5/2, 7/3, 1/2, 5/9}}),map(QQ^4,QQ^4,{{1/5, 9/7, 3/2, 1/4},
      {5/4, 1/5, 1/2, 3/10}, {2, 1, 9/2, 4/3}, {5/2, 1, 3/4,
      3/4}}),map(QQ^4,QQ^4,{{1/7, 5, 9/7, 3/5}, {5/9, 1/3, 1/3, 6/7}, {9/10,
      4, 7/4, 8/9}, {3/5, 5, 1/3, 1}}),map(QQ^4,QQ^4,{{10/7, 1/3, 8/9, 7/2},
      {1/2, 6, 3/2, 1/3}, {1/5, 9/8, 10/7, 7/5}, {9/8, 6/7, 1/7,
      8/3}}),map(QQ^4,QQ^4,{{1/7, 5, 9/7, 3/5}, {5/9, 1/3, 1/3, 6/7}, {9/10,
      4, 7/4, 8/9}, {3/5, 5, 1/3, 1}})}
maxMLE(L,V)

    if #L1==0 then  error("No critical points to evaluate");
    if #L1==1 then  E:=inverse L_0;
    if #L1>=1 then 
    	eval:=for K in L1 list log det K- trace (S*K);
	evalReal:={}
	for pt in eval do (if isReal pt then evalReal=evalReal  | {pt})
	if #evalReal==0 then  error("No critical point evaluates to a real solution");
	indexOptimal:=position(eval, i ->i== max evalReal)
	E=inverse L1_indexOptimal
    return E;
    
    if #L==0 then  error("No critical points to evaluate")
    if #L==1 then  E:=inverse L_0
    if #L>=1 then 
    	eval:=for K in L list log det K- trace (S*K)
	evalReal:={}
	for pt in eval do (if isReal pt then evalReal=evalReal  | {pt})
	if #evalReal==0 then  error("No critical point evaluates to a real solution")
	indexOptimal=positions(eval, i ->i== max evalReal)
        E={};
	for i in indexOptimal do E=E | {inverse L_i};
    return E; 
--------------------------
-- Testing LMG function
--------------------------
-- Test 1: Input
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
U1=matrix{{1,2,1,-1},{2,1,3,0},{-1, 0, 1, 1},{-5, 3, 4, -6}}

-- Test 1: Confirm the function gives correct output
J=scoreEquations(R,U);
I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J===I)

V=sampleCovarianceMatrix U
Jcov=scoreEquationsFromCovarianceMatrix(R,V,sampleData=>false);
I===Jcov

-- Test1: Input covariance matrix directly
restart
debug loadPackage "GraphicalModelsMLE"
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
U=sampleCovarianceMatrix U

flag=false;

 L = directedEdgesMatrix R;
    -- K 
 K= undirectedEdgesMatrix R;
    -- Psi
 P = bidirectedEdgesMatrix R;
    
 d = numRows L;
 (F,lpR)=changeRing(d,R);
    -- create its fraction field
    FR = frac(lpR);
    
    -----------------------------------------------------
    -- Construct Omega
    -----------------------------------------------------
    -- Kinv
    K=substitute(K, FR);
    Kinv=inverse K;
    P=substitute(P,FR);
       
     --Omega
    if K==0 then W=P else (if P==0 then W=Kinv else W = directSum(Kinv,P));
    
    -- move to FR, the fraction field of lpR
    L= substitute (L,FR);
    
    -- Sigma
    if L==0 then S=W else (
	IdL = inverse (id_(FR^d)-L);
    	S = (transpose IdL) * W * IdL
	);
    if S == Kinv then Sinv= K else Sinv = inverse S; 
    
    -- Sample covariance matrix
    if flag then V = sampleCovarianceMatrix(U) else V=U
    
     V = sampleCovarianceMatrix(U) 
    -- Compute ideal J   
    C1 = trace(Sinv * V)
C1derivative := jacobianMatrixOfRationalFunction(C1);
    LL :=jacobianMatrixOfRationalFunction (det Sinv)*matrix{{1/det(Sinv)}} - C1derivative;
    LL=flatten entries(LL);
    denoms := apply(#LL, i -> lift(denominator(LL_i), lpR));

restart
debug needsPackage "GraphicalModelsMLE"
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U1=matrix{{1,2,1,-1},{2,1,3,0},{-1, 0, 1, 1},{-5, 3, 4, -6}}
J1=scoreEquations(R,U1);

J1=scoreEquations(U1,R);

U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J =scoreEquations(R,U);
J =scoreEquations(U,R);

I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J1===I)
assert(J===I)


J2=scoreEquations(R,U1,doSaturate=>false);
-- Test 1: Confirm no saturation works
JnoSat=scoreEquations(R,U1,doSaturate=>false);
dim JnoSat
degree JnoSat

--Test 1: Confirm saturation options work
JSatOpts=scoreEquations(R,U1,doSaturate=>true, saturateOptions => {DegreeLimit=>1, MinimalGenerators => false});
dim JSatOpts
degree JSatOpts

   
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
J=scoreEquationsFromCovarianceMatrix(R, U);
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
 debug loadPackage("GraphicalModelsMLE")
 t0= cpuTime();
 G=graph{{1,2},{2,5},{5,6},{2,4},{4,5},{3,4}}
 g=mixedGraph G
 R2=gaussianRing(g)
 U=matrix{{1, 2, 9, 6, 0, 0}, {2, 7, 7, 3, 2, 2}, {6, 3, 4, 1, 5, 5}, {5, 5, 8, 8, 7, 6}, {3, 2, 3, 8, 7, 5}, {8, 0, 5, 3, 8, 5}}
 --J2=time scoreEquationsFromCovarianceMatrix(R2,U)
 J3= scoreEquationsFromCovarianceMatrix(R,U,doSaturate =>false)
 
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
    C1derivative = jacobianMatrixOfRationalFunction(C1);
    LL =jacobianMatrixOfRationalFunction (det S)*matrix{{(-1/(2*det(S)))}} - (C1derivative);
    LL=flatten entries(LL);
    denoms = apply(#LL, i -> lift(denominator(LL_i), lpR))
    prod = product(denoms);
    J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR))
    t1=cpuTime();
    t1-t0
    -- running line- by-line 484.877033618
    
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
-------------------------------------------------
-- Tests for solverMLE
------------------------------------------------
restart
debug needsPackage "GraphicalModelsMLE"
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
solverMLE(G,U)

restart
debug needsPackage "GraphicalModelsMLE"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}


R.?graph

------------------------------------------------
-- Tests for gaussianRing, ...EdgesMatrix
------------------------------------------------
restart
debug loadPackage "GraphicalModelsMLE"

G=graph{{1,2},{1,3},{2,3}}
D=digraph{{1,6},{4,7}}
B=bigraph{{5,6},{6,7}}

R=gaussianRing G
R.?graph
R.?mixedGraph

R=gaussianRing D
R.?digraph
R.?mixedGraph
R.?gaussianRingData

restart
debug loadPackage "GraphicalModelsMLE"
B=bigraph{{5,6},{6,7}}
R=gaussianRing B

R.?bigraph
R.?mixedGraph

R.?gaussianRingData
restart
debug loadPackage "GraphicalModelsMLE"
G=graph{{1,2},{1,3},{2,3}}
D=digraph{{1,6},{4,7}}
B=bigraph{{5,6},{6,7}}
G1=mixedGraph(G,B)
R=gaussianRing G1

R.?graph
R.?mixedGraph
R.?gaussianRingData

G1=mixedGraph(B,G)
R=gaussianRing G1
R.?graph
R.?mixedGraph

G1=mixedGraph(G,D,B)
R=gaussianRing G1
R.?graph
R.?mixedGraph

G1=mixedGraph(D,G,B)
R=gaussianRing G1
R.?graph
R.?mixedGraph

restart
debug needsPackage "GraphicalModelsMLE"
G=graph{{1,2}}
D=digraph{{1,3},{3,2},{6,7},{7,8},{6,8}}
B=bigraph{{5,4}}


R.?graph
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

myFunction3 = method(Options =>  {test => true});
myFunction3 (ZZ):= opts -> (a)->(
     return opts
    ); 

myFunction3 (2)
myFunction3 (2, test=>true)
myFunction3 (2, test=>false)

myFunction4  = method( Options =>{doFactor => true});
myFunction4(ZZ) := opts ->(a) -> (
    if opts.doFactor then
    (    print "Enter if";
	a=factor a);
    return a;
    );

myFunction4 (2)
myFunction4 (2, doFactor=>true)
myFunction4 (78998989898798982, doFactor=>false)

restart
installPackage "GraphicalModelsMLE"
viewHelp checkPSD

restart
debug loadPackage "GraphicalModelsMLE"
myFunction=method(Options =>{doSaturate => true, saturateOptions => options saturate, sampleData=>true})
myFunction(Ring,List):=  opts ->(G,U) -> (
        if opts.doSaturate then (
	 argSaturate1:=opts.saturateOptions  >>newOpts-> args ->(args, newOpts);--,sampleData=>opts.sampleData
	 << argSaturate1(R,U);
         J:=scoreEquationsFromCovarianceMatrix(argSaturate1(R,U)););
    --else J= scoreEquationsFromCovarianceMatrix(R,U,doSaturate=>false, sampleData=>opts.sampleData);
     return J;
     );
 
 restart
 R = ZZ/32003[a..d];
 I = ideal(a^3-b, a^4-c)
 Ih = homogenize(I,d)
 doSaturate=true
 saturateOptions= options saturate
 argSaturate:=saturateOptions  >>newOpts-> args ->(args, newOpts);
 saturateOptions= options saturate
 argSaturate(Ih,d)
 saturate(argSaturate(Ih,d))


if doSaturate then 
    (   argSaturate:=opts.saturateOptions  >>newOpts-> args ->(args, newOpts);
    	 
	    J=saturate(argSaturate(J,denoms_i))); 
	);
    
-- Test 1: Input
G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
--U1=matrix{{1,2,1,-1},{2,1,3,0},{-1, 0, 1, 1},{-5, 3, 4, -6}}

-- Test 1: Confirm the function gives correct output
J=myFunction(R,U);
I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J===I)

-------------------------------
--scoreEquationsInternal
-------------------------------

restart
loadPackage "Graphs"
loadPackage "StatGraphs"
loadPackage "GraphicalModels"
loadPackage "GraphicalModelsMLE"

debug needsPackage "GraphicalModelsMLE"

scoreEquationsInternal={doSaturate => true, saturateOptions => options saturate, sampleData=>true}>>opts->(R,U)->(
    << opts;
    ----------------------------------------------------
    -- Previous checks 
    ----------------------------------------------------
    if R.?graph then return scoreEquationsUndir(R,U,opts);
    if not R.?mixedGraph then error "Expected a ring created with gaussianRing of a Graph or MixedGraph";
    if not #U==#vertices R.mixedGraph then error "Size of sample data does not match the MixedGraph.";
    ----------------------------------------------------
    -- Extract information about the graph
    ---------------------------------------------------- 
    -- Lambda
    L := directedEdgesMatrix R;
    -- Psi
    P := bidirectedEdgesMatrix R;
    
    -- If the mixedGraph only has undirected part, call specific function for undirected.
    if L==0 and P==0 then 
    return scoreEquationsUndir(R,U,opts);
       
    -- K 
    K := undirectedEdgesMatrix R;
    ----------------------------------------------------
    -- Create an auxiliary ring and its fraction field
    -- which do not have the s variables
    ----------------------------------------------------
    -- create a new ring, lpR, which does not have the s variables
    lpR:=coefficientRing(R)[gens R-set support covarianceMatrix R];
    -- create its fraction field
    FR := frac(lpR);
    
    -----------------------------------------------------
    -- Construct Omega
    -----------------------------------------------------
    -- Kinv
    K=sub(K, FR);
    Kinv:=inverse K;
    P=sub(P,FR);
       
     --Omega
    if K==0 then W:=P else (if P==0 then W=Kinv else W = directSum(Kinv,P));
    
    -- move to FR, the fraction field of lpR
    L= sub(L,FR);
    
    -- Sigma
    d:=numcols L;
    if L==0 then S:=W else (
	IdL := inverse (id_(FR^d)-L);
    	S = (transpose IdL) * W * IdL
	);
    Sinv := inverse S; 
      
    -- Sample covariance matrix
    if opts.sampleData then V := sampleCovarianceMatrix(U) else (
       --If the input is a list of lists we convert it into a list of matrices
       if class U_0===List then U=apply(#U, i -> matrix{U_i});
       -- convert list of matrices into a matrix
       V=U_0;
       for i from 1 to  #U-1 do V= V||U_i
       ); 
   
    -- Compute ideal J   
    C1 := trace(Sinv * V);
    C1derivative := jacobianMatrixOfRationalFunction(C1);
    LL :=jacobianMatrixOfRationalFunction (det Sinv)*matrix{{1/det(Sinv)}} - C1derivative;
    LL=flatten entries(LL);
    denoms := apply(#LL, i -> lift(denominator(LL_i), lpR));
    J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR));
    --Saturate
    if opts.doSaturate then (
        argSaturate:=opts.saturateOptions  >>newOpts-> args ->(args, newOpts);
    	for i from 0 to (#denoms-1) do (
	    if degree denoms_i =={0} then J=J else  
	    J=saturate(argSaturate(J,denoms_i))
	    ); 
	);
    return (J,Sinv);
);

G = mixedGraph(digraph {{1,2},{1,3},{2,3},{3,4}},bigraph {{3,4}})
R=gaussianRing(G)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
(J,SInv)=scoreEquationsInternal(R,U);
J
SInv
I=ideal(20*p_(3,4)+39,50*p_(4,4)-271,440104*p_(3,3)-742363,230*p_(2,2)-203,16*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J===I)

(J,SInv)=scoreEquationsInternal(R,U,doSaturate=>false)
