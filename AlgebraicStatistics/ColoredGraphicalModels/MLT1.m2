-----------------------------------------------------------
--3-cycle with 3 vertices equal
-----------------------------------------------------------

restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_2,l_3},{l_2,l_1,l_4},{l_3,l_4,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Boundary components
(V,n,K2)=embeddedK(K);
BC2=boundaryComponents(K2,3,n)
netList BC2
BC1=boundaryComponents(K2,2,n)
netList BC1

--Elimination ideal
use Rtotal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
codim IG1,dim IG1,degree IG1
netList IG1_*
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

sub((BC2_0)_0,Rtotal)==sub(IG1_0,Rtotal)

-- Change of signs
use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--changes

empiricalMLEexistence(1,50,K)
empiricalMLEexistence(2,50,K)

--gradient ideal and Hessian
use R
f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI --3
resol=resolution(gradI)
rank(resol.dd_2) --3 (4-1)
betti(resol) --l=0
H=jacobian ideal flatten entries jacobian ideal f
rank H --4


-----------------------------------------------------------
--4-cycle with 4 vertices equal
-----------------------------------------------------------

restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_2,0,l_3},{l_2,l_1,l_4,0},{0,l_4,l_1,l_5},{l_3,0,l_5,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*


RS=QQ[support S]
gens RS
S=sub(S,RS)
J4=ideal 1_RS
dim J4, codim J4, degree J4 ---why?
J3=ideal det S;
dim J3, codim J3, degree J3
J2=minors(3,S);
dim J2, codim J2, degree J2
J1=minors(2,S);
dim J1, codim J1, degree J1


--Elimination ideal
use Rtotal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
codim IG1,dim IG1,degree IG1
netList IG1_*
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

-- Change of signs
use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--changes

time empiricalMLEexistence(1,1000,K)
-- used 2.54212 seconds
-- 100 out of 100
-- used 25.6145 seconds
-- 1000 out of 1000
empiricalMLEexistence(2,100,K)


-----------------------------------------------------------
--5-cycle with 5 vertices equal
-----------------------------------------------------------

restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_2,0,0,l_3},{l_2,l_1,l_4,0,0},{0,l_4,l_1,l_5,0},{0,0,l_5,l_1,l_6},{l_3,0,0,l_6,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
use Rtotal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
codim IG1,dim IG1,degree IG1
netList IG1_*
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

-- Change of signs
use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--changes

time empiricalMLEexistence(1,5,K)
-- used 2752.62 seconds
-- 1 out of 1
-- used 15378.3 seconds
-- 5 out of 5
empiricalMLEexistence(2,100,K)



-----------------------------------------------------------
--6-cycle with 6 vertices equal
-----------------------------------------------------------

restart
load "functions.m2"
R=QQ[l_1..l_7]
K=matrix{{l_1,l_2,0,0,0,l_3},{l_2,l_1,l_4,0,0,0},{0,l_4,l_1,l_5,0,0},{0,0,l_5,l_1,l_6,0},{0,0,0,l_6,l_1,l_7},{l_3,0,0,0,l_7,l_1}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

--Elimination ideal
use Rtotal
rk=1;
IG1=sub(rankProjection(stats,rk,S),ring(suffStat(K)));
codim IG1,dim IG1,degree IG1
netList IG1_*
IG2=sub(rankProjection(stats,2,S),ring(suffStat(K)));
IG2

-- Change of signs
use ring IG1
netList empiricalVanishingPolynomials(IG1,3,10000,p,n,stats)
--changes

time empiricalMLEexistence(1,5,K)
--Too many heap sections: Increase MAXHINCR or MAX_HEAP_SECTS

empiricalMLEexistence(2,100,K)
