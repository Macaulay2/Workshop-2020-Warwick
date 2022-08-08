-- GRAPH 1
restart
load "functions.m2"
-- initialize all functions in the beginning of the file, then run this code
n=3
R=QQ[l_1..l_n]
K=matrix{{l_1,l_2,0,l_2},{l_2,l_1,l_3,0},{0,l_3,l_1,l_2},{l_2,0,l_2,l_1}}

(V,n,K2)=embeddedK(K);
--rk 3 component
BC4=boundaryComponents(K2,4,n) --2 degree 2 pols
BC3=boundaryComponents(K2,3,n) -- 2 linear, 1 degree 2
BC2=boundaryComponents(K2,2,n)
netList BC4
netList BC3

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
--p=1 for rk=1,2,3

-- GRAPH 7
restart
-- initialize all functions in the beginning of the file, then run this code
n=5
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_3},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_3,0,l_5,l_1}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --1 degree 4 pol
BC3=boundaryComponents(K2,3,n) --1 degree 4 pol
BC2=boundaryComponents(K2,2,n) --1 linear pol
BC1=boundaryComponents(K2,1,n) 
netList BC4
netList BC3
netList BC2

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList IG1_*
sub(IG1_0,Rtotal)==sub((BC4_0)_0,Rtotal)

--p=1 for rk=2,3, 0<p<1 for rk=1

-- GRAPH 9
restart
-- initialize all functions in the beginning of the file, then run this code
n=6
R=QQ[l_1..l_n]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) --2 degree 3 pol,1 quadratic
BC2=boundaryComponents(K2,2,n) --2 linear 
BC1=boundaryComponents(K2,1,n) 
netList BC3
netList BC2

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList IG1_* --2 min gens
sub(IG1_0,Rtotal)==sub((BC3_3)_0,Rtotal)
sub(IG1_1,Rtotal) % sub(BC3_0+BC3_1,Rtotal)

--p=1 for rk=2,3, we think NO for rk=1 based on different criteria
-- I) non-vanishing of one of the generators of IG1
-- II) no boundary components for rk 3 (complementary of rk 1 in the dual)

-- GRAPH 16
restart
-- initialize all functions in the beginning of the file, then run this code
n=7
R=QQ[l_1..l_n]
K=matrix{{l_1,l_4,0,l_7},{l_4,l_1,l_5,0},{0,l_5,l_2,l_6},{l_7,0,l_6,l_3}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) --didn't finish in a few minutes
BC2=boundaryComponents(K2,2,n)  
BC1=boundaryComponents(K2,1,n) 

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList (trim IG1)_* --7 minimal generators

--check changes of sign of the generators in IG1

use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--there are two generators that don't change sign

--p=1 for rk=2,3, NO for rk=1 according to Uhler's table
-- also consistent with our criteria
-- I) non-vanishing of at least one of the generators of IG1
-- II) no boundary components for rk 3 (complementary of rk 1 in the dual)


-- GRAPH 11
restart
n=5
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_4,0,l_5,l_2}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) -- 3 quadratic
BC2=boundaryComponents(K2,2,n) --empty
BC1=boundaryComponents(K2,1,n) 

netList BC3

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList (trim IG1)_* --principal

sub(IG1_0,Rtotal)==sub((BC3_0)_0,Rtotal)

--check changes of sign of the generators in IG1

use ring IG1
(L1,L2)=differentSign(IG1_0,3,10000,p,n,stats)
--there is a single generator that doesn't change sign

--p=1 for rk=2,3, we think NO for rk=1 based on different criteria:
-- I) non-vanishing of the single generators of IG1
-- II) no boundary components for rk 3 (complementary of rk 1 in the dual)



-- GRAPH 13
restart
n=6
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_6,0,l_5,l_2}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) --2 quartics, 1 quadratic
BC2=boundaryComponents(K2,2,n) --4 linear
BC1=boundaryComponents(K2,1,n)

netList BC3 
netList BC2

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList (trim IG1)_* --2 minimal generators

--check changes of sign of the generators in IG1

use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--both generators change signs

--p=1 for rk=2,3, 0<p<1
-- NOT CONSISTENT WITH OUR BOUNDARY COMPONENT CRITERIA BECAUSE BC4 IS EMPTY!!!!!
-- BUT STILL CONSISTENT WITH THE SIGN CHANGING APPROACH

-- GRAPH 14
restart
n=6
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_6,0,l_5,l_2}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) -- 3 quartics
BC2=boundaryComponents(K2,2,n) -- empty
BC1=boundaryComponents(K2,1,n)

netList BC3 

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList (trim IG1)_* -- 2 minimal generators

--check changes of sign of the generators in IG1

use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--one doesn't change signs

--p=1 for rk=2,3, we think NO for rk=1

-- GRAPH 15
restart
n=6
R=QQ[l_1..l_n]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_1,l_4,0},{0,l_4,l_1,l_5},{l_6,0,l_5,l_2}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) -- 1 quadratic, 1 octic
BC2=boundaryComponents(K2,2,n) -- 1 linear
BC1=boundaryComponents(K2,1,n)

netList BC3 
netList BC2

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList (trim IG1)_* -- 4 minimal generators

--check changes of sign of the generators in IG1

use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--all 4 minimal generators change signs

--p=1 for rk=2,3, 0<p<1
-- NOT CONSISTENT WITH OUR BOUNDARY COMPONENT CRITERIA BECAUSE BC4 IS EMPTY!!!!!
-- BUT STILL CONSISTENT WITH THE SIGN CHANGING APPROACH

-- GRAPH 17
restart
n=7
R=QQ[l_1..l_n]
K=matrix{{l_1,l_4,0,l_7},{l_4,l_2,l_5,0},{0,l_5,l_1,l_6},{l_7,0,l_6,l_3}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) --1 quartic
BC2=boundaryComponents(K2,2,n) --2 linear
BC1=boundaryComponents(K2,1,n) 

netList BC3
netList BC2

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList (trim IG1)_* -- 5 minimal generators
betti (trim IG1) 
betti (trim BC3_1)

--Check better here and in other examples:
--no gen of IG1 belongs to BC3_1

--check changes of sign of the generators in IG1

use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--there are 3 pols that don't change sign

--p=1 for rk=2,3, we think NO because:
-- I) non-vanishing of at least one of the generators of IG1
-- II) no boundary components for rk 3 (complementary of rk 1 in the dual)

-- GRAPH 18
restart
n=8
R=QQ[l_1..l_n]
K=matrix{{l_1,l_5,0,l_8},{l_5,l_2,l_6,0},{0,l_6,l_3,l_7},{l_8,0,l_7,l_4}}

load "functions.m2"
(V,n,K2)=embeddedK(K);
BC4=boundaryComponents(K2,4,n) --empty
BC3=boundaryComponents(K2,3,n) --didn't finish in a few minutes
BC2=boundaryComponents(K2,2,n)  
BC1=boundaryComponents(K2,1,n) 

(p,n,Rtotal,S)=coloredData(K)
stats=sub(suffStat(K),Rtotal)
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)))
netList (trim IG2)_* -- principal

IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)))
netList (trim IG1)_* -- 7 minimal generators

--check changes of sign of the generators in IG1 and IG2

use ring IG2
netList empiricalVanishingPolynomials(IG2,3,10000,p,n,stats)
--changes sign

use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--there are 4 gens that don't change sign

--p=1 for rk=3, NO for rk=1 according to Uhler's table
-- consistent with our vanishing criteria
-- I) non-vanishing of at least one of the generators of IG1
-- II) no boundary components for rk 3 (complementary of rk 1 in the dual)

-- 0<p<1 for rk=2 according to Uhler
--consistent with vanishing criteria, 
-- but don't know about empty boundary component for rk2 since it didn't finish
