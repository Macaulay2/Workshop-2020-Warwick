-- 3-cycle: 2 vertices and 2 edges including edge between the two vertices
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_3,l_3},{l_3,l_1,l_4},{l_3,l_4,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_22,d_23,d_33]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,3)
netList 

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,0)-Dadj_(1,1)
f2=Dadj_(0,1)-Dadj_(0,2)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33,t_3=>2*d_12+2*d_13,t_4=>2*d_23})
netList eq_*

MLEeq=ideal{f1}+ideal{f2}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
isPrime MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
netList minimalPrimes MLEeq1
isPrime MLEeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
netList (L_2)_*
codim L_2,degree L_2
det D % gb L_2


MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq2)
det D % gb MLEeq2d
--nice grobner basis

-- 3-cycle: 2 vertices and 2 edges not including edge between the two vertices
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_3,l_4},{l_3,l_1,l_4},{l_4,l_4,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_22,d_23,d_33]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,3)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,0)-Dadj_(1,1)
f2=Dadj_(1,2)-Dadj_(0,2)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33,t_3=>2*d_12,t_4=>2*d_13+2*d_23})
netList eq_*

MLEeq=ideal{f1}+ideal{f2}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
isPrime MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
L1=minimalPrimes MLEeq1
length L1
L1_0==MLEeq1
L1_0==radical(MLEeq1)
isPrime MLEeq1
isPrimary MLEeq1 --says false but doesn't make sense
det D % gb MLEeq1  ---0

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
isPrime MLEeq1d
betti (trim MLEeq1d)
netList (L_0)_*
det D % gb L_0 --0
isSubset(minors(2,D),MLEeq1d)

--infeasibility certificate
IG1=sub(rankProjection(stats,1,S),ring(suffStat(K)));
netList IG1_*
f=IG1_0
factor f

RL=QQ[t_1..t_4,u11,u12,u13,u22,u23,u33]
LT=matrix{{u11,0,0},{u12,u22,0},{u13,u23,u33}}
PSD=LT*transpose(LT)
f=sub(f,RL)
fu=sub(f,{t_1=>PSD_(0,0)+PSD_(1,1),t_2=>PSD_(2,2),t_3=>2*PSD_(0,1),t_4=>2*PSD_(0,2)+2*PSD_(1,2)})

loadPackage "SemidefiniteProgramming"
loadPackage "SumsOfSquares"
sol=solveSOS (-fu)
peek sol
s=sosPoly sol
peek s
coefs=s#coefficients
gene=s#generators
prod=apply(coefs,gene,(i,j)->i*j^2)
sos=sum(prod)
-fu==sos --true


-- 3-cycle with 3 edges equal
restart
load "functions.m2"
R=QQ[l_1..l_4]
K=matrix{{l_1,l_4,l_4},{l_4,l_2,l_4},{l_4,l_4,l_3}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_22,d_23,d_33]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,3)
netList 

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,1)-Dadj_(1,2)
f2=Dadj_(0,1)-Dadj_(0,2)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11,t_2=>d_22,t_3=>d_33,t_4=>2*d_12+2*d_13+2*d_23})
netList eq_*

MLEeq=ideal{f1}+ideal{f2}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
isPrime MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
netList minimalPrimes MLEeq1
isPrime MLEeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
length L
netList (L_1)_*
netList (L_0)_*
codim L_1,degree L_1
det D % gb L_1


MLEeq1s=trim eliminate({d_11, d_12, d_13, d_22, d_23, d_33},MLEeq1)
netList MLEeq1s_* 
MLEeq1s==minors(2,S)

MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq2)
det D % gb MLEeq2d

--nice grobner basis

-- 3-cycle with 2 vertices equal
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,l_4},{l_3,l_1,l_5},{l_4,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_22,d_23,d_33]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,3)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f=Dadj_(0,0)-Dadj_(1,1)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33,t_3=>2*d_12,t_4=>2*d_13,t_5=>2*d_23})
netList eq_*

MLEeq=ideal{f}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

mingens MLEeq
netList minimalPrimes MLEeq
isPrime MLEeq
codim MLEeq

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S)
netList (trim MLEeq1)_*
netList minimalPrimes MLEeq1

toString gens Raux
MLEeq1d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq1)
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
L=minimalPrimes MLEeq1d
length L
netList (L_1)_*
netList (L_0)_*
codim L_1,degree L_1
det D % gb L_1


MLEeq1s=trim eliminate({d_11, d_12, d_13, d_22, d_23, d_33},MLEeq1)
netList MLEeq1s_* 
MLEeq1s==minors(2,S)

MLEeq2=MLEeq+minors(3,S)
netList (trim MLEeq2)_*
netList minimalPrimes MLEeq2
isPrime MLEeq2
MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(2,2), s_(2,3), s_(3,3)},MLEeq2)
det D % gb MLEeq2d



--4-cycle: GRAPH 9
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,4)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(1,1)-Dadj_(3,3)
f2=Dadj_(0,1)-Dadj_(0,3)
f3=Dadj_(0,2)
f4=Dadj_(1,3)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11,t_2=>d_22+d_44,t_3=>d_33,t_4=>2*d_12+2*d_14,t_5=>2*d_23,t_6=>2*d_34})
netList eq_*

MLEeq=ideal{f1,f2,f3,f4}+eq
netList MLEeq_*    
netList (trim MLEeq)_*

--mingens MLEeq
--L=minimalPrimes MLEeq
--netList L
--codim MLEeq

S=sub(S,Raux)
MLEeq1=trim(MLEeq+minors(2,S))
netList MLEeq1_*
--L1=time minimalPrimes MLEeq1
--netList L1
--length L1
--netList (L1_0)_*
--netList (L1_1)_*
--netList (L1_2)_*
--netList (L1_3)_*
--netList (L1_4)_*
--netList (L1_5)_*
--det D % gb L1_5

---Alternative stuff----------------------------------------------------------------------
J=ideal{f1,f2}
reciprocal=time saturate(J,det D);
netList reciprocal_*
isPrime reciprocal
test3=trim(reciprocal+minors(4,D)) 
--are all rank 3 matrices are in the reciprocal variety?
-- true if and only if reciprocal is contained in det(D), i.e. satisfying the rk3 condition
--already implies satisfying the condition to be in the reciprocal
isSubset(reciprocal,minors(4,D))
Ltest3=time minimalPrimes test; 
netList Ltest3
netList (Ltest3_0)_*
netList (Ltest3_1)_*
test2=trim(reciprocal+minors(3,D)) 
isSubset(reciprocal,minors(3,D))
--not all rank 2 matrices are in the reciprocal variety
Ltest2=time minimalPrimes test2; 
netList Ltest2
netList (Ltest2_0)_*
netList (Ltest2_1)_*

test1=trim(reciprocal+minors(2,D))
isPrime test1
netList (test1)_*
isSubset(minors(2,D),reciprocal)
test1==minors(2,D) --all rank 1 matrices are in the reciprocal variety, 
--hence all rk 1 PSD matrices are in the closure of K_G^{-1}. Is that actually true?
--That would be saying that K_G=reciprocal variety intersected with the cone of PD matrices.
adjeq=ideal{f1,f2}
adjeq1=adjeq+minors(2,S)
adjeq1==MLEeq1
L1adj=time minimalPrimes adjeq1
length L1adj
netList (L1adj_0)_*
netList (L1adj_1)_*
det D % gb L1adj_0
det D % gb L1adj_1
det D % gb adjeq
det D % gb adjeq1
------------------------------------------------------------------------------------

MLEeq1d=time eliminate({s_(1,1),s_(1,2),s_(1,3),s_(1,4),s_(2,2),s_(2,3),
	s_(2,4),s_(3,3),s_(3,4),s_(4,4)},MLEeq1);
-- used 7.32812 seconds
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
--(4,3)
betti (trim MLEeq1d)
--1 quadric, 5 cubics, 1 quintic, 6 sextics, 8 septics, 3 octics
L1d=time minimalPrimes MLEeq1d;
-- used 192.828 seconds
length L1d
netList (L1d_0)_*
--d33=0
netList (L1d_1)_*
--d11=0
netList (L1d_2)_*
--d13^2-d_1d33
netList (L1d_3)_*
det D % gb L1d_3 --0
netList (L1d_4)_*
--d34^2-d33d44


--MLEeq1s=trim eliminate({d_11, d_12, d_13,d_14, d_22, d_23,d_24,d_33,d_34,d_44},MLEeq1)
--netList MLEeq1s_* 
--MLEeq1s==minors(2,S)

--MLEeq2=MLEeq+minors(3,S)
--netList (trim MLEeq2)_*
--netList minimalPrimes MLEeq2
--isPrime MLEeq2
--MLEeq2d=trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq2)
--det D % gb MLEeq2d


-- Frets' head
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_4,0,l_5,l_2}}
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,4)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,0)-Dadj_(1,1); --lambda1 from Uhler's paper
f2=Dadj_(2,2)-Dadj_(3,3); --lambda2
f3=Dadj_(0,3)-Dadj_(1,2); --lambda4
f4=Dadj_(0,2)
f5=Dadj_(1,3)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33+d_44,t_3=>2*d_12,t_4=>2*d_23+2*d_14,t_5=>2*d_34})
netList eq_*

MLEeq=ideal{f1,f2,f3,f4,f5}+eq;
netList (trim MLEeq)_*

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S);
MLEeq1d=time trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq1);
-- used 108.453 seconds
ideal(d_14^2*d_33-d_24^2*d_33-2*d_13*d_14*d_34+2*d_23*d_24*d_34+d_11*d_34^2-d_22*d_34^
      2+d_13^2*d_44-d_23^2*d_44-d_11*d_33*d_44+d_22*d_33*d_44,d_14*d_23*d_24-d_13*d_24^2-d_
      14*d_22*d_34+d_12*d_24*d_34+d_13*d_22*d_44-d_12*d_23*d_44,d_14^2*d_23-d_14*d_23^2-d_13
      *d_14*d_24+d_13*d_23*d_24+d_14*d_22*d_33-d_12*d_24*d_33-d_12*d_14*d_34-d_13*d_22*d_34+
      d_12*d_23*d_34+d_11*d_24*d_34+d_12*d_13*d_44-d_11*d_23*d_44,d_13*d_14*d_23-d_13^2*d_24
      -d_12*d_14*d_33+d_11*d_24*d_33+d_12*d_13*d_34-d_11*d_23*d_34,d_13^2*d_22-d_14^2*d_22-2
      *d_12*d_13*d_23+d_11*d_23^2+2*d_12*d_14*d_24-d_11*d_24^2+d_12^2*d_33-d_11*d_22*d_33-d_
      12^2*d_44+d_11*d_22*d_44,d_14^4+14*d_14*d_23^3+d_23^4+4*d_13*d_14^2*d_24-10*d_13*d_23^
      2*d_24+6*d_13^2*d_24^2-2*d_11*d_14*d_23*d_33-12*d_14*d_22*d_23*d_33-d_11*d_23^2*d_33-d
      _22*d_23^2*d_33+10*d_12*d_23*d_24*d_33+3*d_11*d_24^2*d_33-5*d_22*d_24^2*d_33+d_12^2*d_
      33^2-2*d_11*d_13*d_14*d_34+4*d_13*d_14*d_22*d_34-2*d_12*d_14*d_23*d_34+10*d_13*d_22*d_
      23*d_34-14*d_12*d_23^2*d_34-6*d_12*d_13*d_24*d_34-4*d_11*d_14*d_24*d_34-12*d_11*d_23*d
      _24*d_34+10*d_22*d_23*d_24*d_34+2*d_11*d_12*d_33*d_34+2*d_12*d_22*d_33*d_34+2*d_11^2*d
      _34^2+6*d_11*d_22*d_34^2-4*d_22^2*d_34^2+d_11*d_13^2*d_44-4*d_12*d_13*d_14*d_44-d_11*d
      _14^2*d_44-6*d_14^2*d_22*d_44-10*d_12*d_13*d_23*d_44+2*d_11*d_14*d_23*d_44-2*d_14*d_22
      *d_23*d_44+13*d_11*d_23^2*d_44-6*d_22*d_23^2*d_44+10*d_12*d_14*d_24*d_44-5*d_11*d_24^2
      *d_44-d_11^2*d_33*d_44+7*d_12^2*d_33*d_44-9*d_11*d_22*d_33*d_44+5*d_22^2*d_33*d_44+2*d
      _11*d_12*d_34*d_44+2*d_12*d_22*d_34*d_44-4*d_12^2*d_44^2+5*d_11*d_22*d_44^2)

netList MLEeq1d_*
--codim MLEeq1d, degree MLEeq1d
--betti MLEeq1d

L1d=time minimalPrimes MLEeq1d;
 -- used 1604.75 seconds
{ideal(d_24^2*d_33-2*d_23*d_24*d_34+d_22*d_34^2+d_23^2*d_44-d_22*d_33*d_44,d_14*d_24*d
      _33-d_14*d_23*d_34-d_13*d_24*d_34+d_12*d_34^2+d_13*d_23*d_44-d_12*d_33*d_44,d_14^2*d_
      33-2*d_13*d_14*d_34+d_11*d_34^2+d_13^2*d_44-d_11*d_33*d_44,d_14*d_23*d_24-d_13*d_24^2-
      d_14*d_22*d_34+d_12*d_24*d_34+d_13*d_22*d_44-d_12*d_23*d_44,d_14*d_23^2-d_13*d_23*d_24
      -d_14*d_22*d_33+d_12*d_24*d_33+d_13*d_22*d_34-d_12*d_23*d_34,d_14^2*d_23-d_13*d_14*d_
      24-d_12*d_14*d_34+d_11*d_24*d_34+d_12*d_13*d_44-d_11*d_23*d_44,d_13*d_14*d_23-d_13^2*d
      _24-d_12*d_14*d_33+d_11*d_24*d_33+d_12*d_13*d_34-d_11*d_23*d_34,d_14^2*d_22-2*d_12*d_
      14*d_24+d_11*d_24^2+d_12^2*d_44-d_11*d_22*d_44,d_13*d_14*d_22-d_12*d_14*d_23-d_12*d_13
      *d_24+d_11*d_23*d_24+d_12^2*d_34-d_11*d_22*d_34,d_13^2*d_22-2*d_12*d_13*d_23+d_11*d_23
      ^2+d_12^2*d_33-d_11*d_22*d_33,d_14^4+d_23^4+4*d_13*d_14^2*d_24+4*d_13*d_23^2*d_24+6*d_
      13^2*d_24^2-2*d_11*d_14*d_23*d_33+2*d_14*d_22*d_23*d_33-d_11*d_23^2*d_33-d_22*d_23^2*d
      _33-4*d_12*d_23*d_24*d_33+d_12^2*d_33^2-2*d_11*d_13*d_14*d_34+2*d_12*d_14*d_23*d_34-4*
      d_13*d_22*d_23*d_34-2*d_12*d_13*d_24*d_34-4*d_11*d_14*d_24*d_34-10*d_11*d_23*d_24*d_34
      +2*d_11*d_12*d_33*d_34+2*d_12*d_22*d_33*d_34+2*d_11^2*d_34^2-4*d_12^2*d_34^2+7*d_11*d_
      22*d_34^2+d_22^2*d_34^2+d_11*d_13^2*d_44-4*d_12*d_13*d_14*d_44-d_11*d_14^2*d_44-10*d_
      12*d_13*d_23*d_44+2*d_11*d_14*d_23*d_44-2*d_14*d_22*d_23*d_44+10*d_11*d_23^2*d_44-d_22
      *d_23^2*d_44-2*d_12*d_14*d_24*d_44+d_11*d_24^2*d_44-d_11^2*d_33*d_44+7*d_12^2*d_33*d_
      44-6*d_11*d_22*d_33*d_44+2*d_11*d_12*d_34*d_44+2*d_12*d_22*d_34*d_44+2*d_12^2*d_44^2-d
      _11*d_22*d_44^2), ideal(d_34-d_44,d_33-d_44,d_23-d_24,d_13-d_14,d_14^2+2*d_14*d_24+d_
      24^2-d_11*d_44-2*d_12*d_44-d_22*d_44),
      ideal(d_34+d_44,d_33-d_44,d_23+d_24,d_13+d_14,d_14^2-2*d_14*d_24+d_24^2-d_11*d_44+2*d_
      12*d_44-d_22*d_44), ideal(d_14-d_24,d_13-d_23,d_12-d_22,d_11-d_22,d_23^2+2*d_23*d_24+d
      _24^2-d_22*d_33-2*d_22*d_34-d_22*d_44),
      ideal(d_14+d_24,d_13+d_23,d_12+d_22,d_11-d_22,d_23^2-2*d_23*d_24+d_24^2-d_22*d_33+2*d_
      22*d_34-d_22*d_44), ideal(d_33-d_44,d_14-d_23,d_13-d_24,d_11-d_22,d_23^2*d_24-d_24^3-d
      _22*d_23*d_34+d_12*d_24*d_34-d_12*d_23*d_44+d_22*d_24*d_44,4*d_23^3-2*d_23*d_24^2-d_22
      *d_24*d_34-d_12*d_24*d_44,2*d_24^4-4*d_12*d_23^2*d_34+3*d_22*d_23*d_24*d_34-2*d_12*d_
      24^2*d_34+d_22^2*d_34^2-4*d_22*d_23^2*d_44+3*d_12*d_23*d_24*d_44-2*d_22*d_24^2*d_44+2*
      d_12*d_22*d_34*d_44+d_12^2*d_44^2)}

length L1d
netList (L1d_0)_*
netList (L1d_1)_*
netList (L1d_2)_*
netList (L1d_3)_*
netList (L1d_4)_*
netList (L1d_5)_*

det D % gb L1d_0 --0
det D % gb L1d_1 --0
det D % gb L1d_2 --0
det D % gb L1d_3 --0
det D % gb L1d_4 --0
det D % gb L1d_5 --not 0

-- compute list of ppal minors 
principalMinors=M->(
--M=sub(M,QQ[support M]);
mm:=numcols M;
ss:=drop(subsets mm,-1);
for s in ss list submatrix'(M,s,s)
);
-- compute list of determinants of ppal minors 
detPrincipalMinors=M->(
L:=principalMinors M;
for l in L list det l
);
netList unique detPrincipalMinors(D)
netList apply(unique detPrincipalMinors(D),i->i % gb L1d_5)
Rd=QQ[support L1d_5]
comp5=sub(L1d_5,Rd);
dim comp5, degree comp5
--(4,8)
netList comp5_*
D=sub(D,Rd)
sub(comp5,{d_11=>1,d_22=>1,d_33=>1,d_44=>1,d_12=>0,d_13=>0,d_14=>0,d_23=>0,d_24=>0,d_34=>0})
--Identity matrix belongs to this component

--GRAPH 11

restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_4,0,l_5,l_2}}
--2+2 opposite vertices, 2 opposite edges

(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,4)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f0=Dadj_(0,0)-Dadj_(2,2) --l1
f1=Dadj_(1,1)-Dadj_(3,3) --l2
f2=Dadj_(1,2)-Dadj_(0,3) --l4
f3=Dadj_(0,2)
f4=Dadj_(1,3)

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_33,t_2=>d_22+d_44,t_3=>2*d_12,t_4=>2*d_14+2*d_23,t_5=>2*d_34})
netList eq_*

MLEeq=ideal{f0,f1,f2,f3,f4}+eq;

S=sub(S,Raux)
MLEeq1=trim(MLEeq+minors(2,S))
netList MLEeq1_*

MLEeq1d=time trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq1);
-- used 301.25 seconds
ideal(d_12^2+d_14^2-d_11*d_22+2*d_14*d_23+d_23^2-d_22*d_33-2*d_12*d_34+d_34^2-d_11*d_
      44-d_33*d_44,d_14*d_23*d_24-d_13*d_24^2-d_14*d_22*d_34+d_12*d_24*d_34+d_13*d_22*d_44-d
      _12*d_23*d_44,d_14^2*d_23-d_14*d_23^2-d_13*d_14*d_24+d_13*d_23*d_24+d_14*d_22*d_33-d_
      12*d_24*d_33-d_12*d_14*d_34-d_13*d_22*d_34+d_12*d_23*d_34+d_11*d_24*d_34+d_12*d_13*d_
      44-d_11*d_23*d_44,d_13*d_14*d_23-d_13^2*d_24-d_12*d_14*d_33+d_11*d_24*d_33+d_12*d_13*d
      _34-d_11*d_23*d_34,d_14^2*d_22-2*d_12*d_14*d_24+d_11*d_24^2-d_24^2*d_33+2*d_23*d_24*d_
      34-d_22*d_34^2-d_14^2*d_44-2*d_14*d_23*d_44-2*d_23^2*d_44+2*d_22*d_33*d_44+2*d_12*d_34
      *d_44-d_34^2*d_44+d_11*d_44^2+d_33*d_44^2,d_13^2*d_22-2*d_12*d_13*d_23+d_11*d_23^2-2*d
      _14^2*d_33-2*d_14*d_23*d_33-d_23^2*d_33+d_22*d_33^2+2*d_13*d_14*d_34+2*d_12*d_33*d_34-
      d_11*d_34^2-d_33*d_34^2-d_13^2*d_44+2*d_11*d_33*d_44+d_33^2*d_44)

--netList MLEeq1d_*
--codim MLEeq1d, degree MLEeq1d
--betti MLEeq1d

L1=time minimalPrimes MLEeq1d;
-- used 17.5 seconds
toString L1
L1={ideal(d_12^2+d_14^2-d_11*d_22+2*d_14*d_23+d_23^2-d_22*d_33-2*d_12*d_34+d_34^2-d_11*d_
      44-d_33*d_44,d_24^2*d_33-2*d_23*d_24*d_34+d_22*d_34^2+d_23^2*d_44-d_22*d_33*d_44,d_14*
      d_24*d_33-d_14*d_23*d_34-d_13*d_24*d_34+d_12*d_34^2+d_13*d_23*d_44-d_12*d_33*d_44,d_14
      ^2*d_33-2*d_13*d_14*d_34+d_11*d_34^2+d_13^2*d_44-d_11*d_33*d_44,d_14*d_23*d_24-d_13*d_
      24^2-d_14*d_22*d_34+d_12*d_24*d_34+d_13*d_22*d_44-d_12*d_23*d_44,d_14*d_23^2-d_13*d_23
      *d_24-d_14*d_22*d_33+d_12*d_24*d_33+d_13*d_22*d_34-d_12*d_23*d_34,d_14^2*d_23-d_13*d_
      14*d_24-d_12*d_14*d_34+d_11*d_24*d_34+d_12*d_13*d_44-d_11*d_23*d_44,d_13*d_14*d_23-d_
      13^2*d_24-d_12*d_14*d_33+d_11*d_24*d_33+d_12*d_13*d_34-d_11*d_23*d_34,d_14^2*d_22-2*d_
      12*d_14*d_24+d_11*d_24^2-d_14^2*d_44-2*d_14*d_23*d_44-d_23^2*d_44+d_22*d_33*d_44+2*d_
      12*d_34*d_44-d_34^2*d_44+d_11*d_44^2+d_33*d_44^2,d_13*d_14*d_22-d_12*d_14*d_23-d_12*d_
      13*d_24+d_11*d_23*d_24-d_14^2*d_34-2*d_14*d_23*d_34-d_23^2*d_34+d_22*d_33*d_34+2*d_12*
      d_34^2-d_34^3+d_11*d_34*d_44+d_33*d_34*d_44,d_13^2*d_22-2*d_12*d_13*d_23+d_11*d_23^2-2
      *d_14*d_23*d_33-d_23^2*d_33+d_22*d_33^2-2*d_13*d_14*d_34+2*d_12*d_33*d_34+d_11*d_34^2-
      d_33*d_34^2+d_13^2*d_44+d_33^2*d_44),
      ideal(d_14-d_23,d_24*d_33-d_13*d_44,d_22*d_33-d_11*d_44,d_22*d_23-d_12*d_24+d_24*d_34-
      d_23*d_44,d_13*d_22-d_11*d_24,d_12*d_13-d_11*d_23+d_23*d_33-d_13*d_34,d_12^2-d_11*d_22
      +4*d_23^2-2*d_12*d_34+d_34^2-2*d_11*d_44-d_33*d_44,2*d_23*d_24*d_34+d_22*d_34^2+2*d_13
      *d_24*d_44-2*d_12*d_34*d_44+d_34^2*d_44-3*d_11*d_44^2-d_33*d_44^2,2*d_13*d_23*d_34-2*d
      _12*d_33*d_34+d_11*d_34^2+d_33*d_34^2+2*d_13^2*d_44-3*d_11*d_33*d_44-d_33^2*d_44,2*d_
      23*d_24^2+d_22*d_24*d_34-3*d_12*d_24*d_44+4*d_24*d_34*d_44-4*d_23*d_44^2,2*d_13*d_24^2
      +d_12*d_24*d_34-3*d_24*d_34^2-3*d_11*d_24*d_44+4*d_23*d_34*d_44-d_13*d_44^2,2*d_12*d_
      24^2+d_22^2*d_34-2*d_24^2*d_34-3*d_12*d_22*d_44+3*d_22*d_34*d_44-d_12*d_44^2,2*d_11*d_
      24^2+d_12*d_22*d_34-3*d_22*d_34^2-3*d_11*d_22*d_44-2*d_13*d_24*d_44+3*d_12*d_34*d_44-d
      _34^2*d_44+2*d_11*d_44^2+d_33*d_44^2,2*d_23^2*d_24+d_12*d_24*d_34-d_24*d_34^2-2*d_12*d
      _23*d_44-d_11*d_24*d_44+2*d_23*d_34*d_44-d_13*d_44^2,2*d_13*d_23*d_24+d_11*d_24*d_34-3
      *d_11*d_23*d_44-d_23*d_33*d_44+d_13*d_34*d_44,2*d_12*d_23*d_24+d_12*d_22*d_34-3*d_11*d
      _22*d_44+8*d_23^2*d_44+2*d_13*d_24*d_44-3*d_12*d_34*d_44+2*d_34^2*d_44-7*d_11*d_44^2-2
      *d_33*d_44^2,2*d_11*d_23*d_24+d_11*d_22*d_34-3*d_11*d_12*d_44-2*d_13*d_23*d_44-d_12*d_
      33*d_44+3*d_11*d_34*d_44,2*d_13^2*d_24+d_11*d_23*d_34+3*d_23*d_33*d_34-2*d_13*d_34^2-3
      *d_11*d_13*d_44-d_13*d_33*d_44,2*d_11*d_13*d_24+d_11*d_12*d_34+3*d_12*d_33*d_34-3*d_11
      *d_34^2-d_33*d_34^2-3*d_11^2*d_44-2*d_13^2*d_44+2*d_11*d_33*d_44+d_33^2*d_44,2*d_13*d_
      23^2-2*d_12*d_23*d_33+d_11*d_23*d_34+d_23*d_33*d_34-d_11*d_13*d_44-d_13*d_33*d_44,2*d_
      11*d_23^2+6*d_23^2*d_33+d_11*d_12*d_34-d_12*d_33*d_34-d_11*d_34^2+d_33*d_34^2-3*d_11^2
      *d_44-4*d_11*d_33*d_44-d_33^2*d_44,2*d_13^2*d_23-3*d_11*d_23*d_33-d_23*d_33^2+d_11*d_
      13*d_34+d_13*d_33*d_34,2*d_11*d_13*d_23-3*d_11*d_12*d_33-2*d_13*d_23*d_33-d_12*d_33^2+
      d_11^2*d_34+3*d_11*d_33*d_34)};

length L1
netList (L1_0)_*
det D % gb L1_0 --0
netList (L1_1)_*
det D % gb L1_1 --not 0

-- compute list of ppal minors 
principalMinors=M->(
--M=sub(M,QQ[support M]);
mm:=numcols M;
ss:=drop(subsets mm,-1);
for s in ss list submatrix'(M,s,s)
);
-- compute list of determinants of ppal minors 
detPrincipalMinors=M->(
L:=principalMinors M;
for l in L list det l
);
netList apply(unique detPrincipalMinors(D),i->i % gb L1_1)
sub(L1_1,{d_11=>1,d_22=>1,d_33=>1,d_44=>1,d_12=>0,d_13=>0,d_14=>0,d_23=>0,d_24=>0,d_34=>0})
--Does the identity matrix belong to this component? No!


J=time saturate(ideal{f0,f1,f2,f3,f4},ideal{det D});
J==ideal{f0,f1,f2,f3,f4}
isSubset(ideal{f0,f1,f2,f3,f4},J)
I=J+eq+minors(2,sub(S,Raux));
time I==MLEeq1 --false
-- used 222.641 seconds
time isSubset(MLEeq1,I) --true
Id=time trim eliminate({s_(1,1),s_(1,2),s_(1,3),s_(1,4),s_(2,2),s_(2,3),s_(2,4),
	s_(3,3),s_(3,4),s_(4,4)},I);
netList Id_*
LId=time minimalPrimes Id;
length LId
netList (trim LId_0)_*
det D % gb LId_0

netList (trim LId_1)_*
det D % gb LId_1
netList apply(unique detPrincipalMinors(D),i->i % gb LId_1)

time LId_0==L1_0 --false
isSubset(L1_0,LId_0) --true The nice component improves!
det D % gb LId_0
det D % gb L1_0

time LId_1==L1_1 --true The ugly component does not look nicer!


LId={ideal(d_14-d_23,d_24*d_33-d_13*d_44,d_22*d_33-d_11*d_44,2*d_23*d_24-d_22*d_34-d_12*d_
      44,2*d_13*d_24-d_12*d_34-d_34^2-d_11*d_44+d_33*d_44,2*d_23^2-d_12*d_34+d_34^2-d_11*d_
      44-d_33*d_44,d_22*d_23-d_12*d_24+d_24*d_34-d_23*d_44,2*d_13*d_23-d_12*d_33-d_11*d_34,d
      _12*d_23-d_11*d_24+d_23*d_34-d_13*d_44,d_13*d_22-d_11*d_24,d_12*d_13-d_11*d_23+d_23*d_
      33-d_13*d_34,d_12^2-d_11*d_22-d_34^2+d_33*d_44),
      ideal(d_14-d_23,d_24*d_33-d_13*d_44,d_22*d_33-d_11*d_44,d_22*d_23-d_12*d_24+d_24*d_34-
      d_23*d_44,d_13*d_22-d_11*d_24,d_12*d_13-d_11*d_23+d_23*d_33-d_13*d_34,d_12^2-d_11*d_22
      +4*d_23^2-2*d_12*d_34+d_34^2-2*d_11*d_44-d_33*d_44,2*d_23*d_24*d_34+d_22*d_34^2+2*d_13
      *d_24*d_44-2*d_12*d_34*d_44+d_34^2*d_44-3*d_11*d_44^2-d_33*d_44^2,2*d_13*d_23*d_34-2*d
      _12*d_33*d_34+d_11*d_34^2+d_33*d_34^2+2*d_13^2*d_44-3*d_11*d_33*d_44-d_33^2*d_44,2*d_
      23*d_24^2+d_22*d_24*d_34-3*d_12*d_24*d_44+4*d_24*d_34*d_44-4*d_23*d_44^2,2*d_13*d_24^2
      +d_12*d_24*d_34-3*d_24*d_34^2-3*d_11*d_24*d_44+4*d_23*d_34*d_44-d_13*d_44^2,2*d_12*d_
      24^2+d_22^2*d_34-2*d_24^2*d_34-3*d_12*d_22*d_44+3*d_22*d_34*d_44-d_12*d_44^2,2*d_11*d_
      24^2+d_12*d_22*d_34-3*d_22*d_34^2-3*d_11*d_22*d_44-2*d_13*d_24*d_44+3*d_12*d_34*d_44-d
      _34^2*d_44+2*d_11*d_44^2+d_33*d_44^2,2*d_23^2*d_24+d_12*d_24*d_34-d_24*d_34^2-2*d_12*d
      _23*d_44-d_11*d_24*d_44+2*d_23*d_34*d_44-d_13*d_44^2,2*d_13*d_23*d_24+d_11*d_24*d_34-3
      *d_11*d_23*d_44-d_23*d_33*d_44+d_13*d_34*d_44,2*d_12*d_23*d_24+d_12*d_22*d_34-3*d_11*d
      _22*d_44+8*d_23^2*d_44+2*d_13*d_24*d_44-3*d_12*d_34*d_44+2*d_34^2*d_44-7*d_11*d_44^2-2
      *d_33*d_44^2,2*d_11*d_23*d_24+d_11*d_22*d_34-3*d_11*d_12*d_44-2*d_13*d_23*d_44-d_12*d_
      33*d_44+3*d_11*d_34*d_44,2*d_13^2*d_24+d_11*d_23*d_34+3*d_23*d_33*d_34-2*d_13*d_34^2-3
      *d_11*d_13*d_44-d_13*d_33*d_44,2*d_11*d_13*d_24+d_11*d_12*d_34+3*d_12*d_33*d_34-3*d_11
      *d_34^2-d_33*d_34^2-3*d_11^2*d_44-2*d_13^2*d_44+2*d_11*d_33*d_44+d_33^2*d_44,2*d_13*d_
      23^2-2*d_12*d_23*d_33+d_11*d_23*d_34+d_23*d_33*d_34-d_11*d_13*d_44-d_13*d_33*d_44,2*d_
      11*d_23^2+6*d_23^2*d_33+d_11*d_12*d_34-d_12*d_33*d_34-d_11*d_34^2+d_33*d_34^2-3*d_11^2
      *d_44-4*d_11*d_33*d_44-d_33^2*d_44,2*d_13^2*d_23-3*d_11*d_23*d_33-d_23*d_33^2+d_11*d_
      13*d_34+d_13*d_33*d_34,2*d_11*d_13*d_23-3*d_11*d_12*d_33-2*d_13*d_23*d_33-d_12*d_33^2+
      d_11^2*d_34+3*d_11*d_33*d_34)};

sub(LId_1,{d_11=>1,d_22=>1,d_33=>1,d_44=>1,d_12=>0,d_13=>0,d_14=>0,d_23=>0,d_24=>0,d_34=>0})
--Does the identity matrix belong to this component? NO!

--With Cholesky decomposition and reciprocal variety
restart
load "functions.m2"
R=QQ[l_1..l_5]
K=matrix{{l_1,l_3,0,l_4},{l_3,l_2,l_4,0},{0,l_4,l_1,l_5},{l_4,0,l_5,l_2}}
--2+2 opposite vertices, 2 opposite edges

(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_(1,1),d_(1,2),d_(1,3),d_(1,4),d_(2,2),d_(2,3),d_(2,4),d_(3,3),d_(3,4),d_(4,4)];
gens Raux
D=genericSymmetricMatrix(Raux,d_(1,1),4)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f0=Dadj_(0,0)-Dadj_(2,2) --l1
f1=Dadj_(1,1)-Dadj_(3,3) --l2
f2=Dadj_(1,2)-Dadj_(0,3) --l4
f3=Dadj_(0,2)
f4=Dadj_(1,3)
J=time saturate(ideal{f0,f1,f2,f3,f4},ideal{det D}); --ideal of the reciprocal variety

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_(1,1)+d_(3,3),t_2=>d_(2,2)+d_(4,4),t_3=>2*d_(1,2),t_4=>2*d_(1,4)+2*d_(2,3),t_5=>2*d_(3,4)})
netList eq_*

R2=QQ[d_(1,1),d_(1,2),d_(1,3),d_(1,4),d_(2,2),d_(2,3),d_(2,4),d_(3,3),d_(3,4),d_(4,4),
    u_(1,1),u_(1,2),u_(1,3),u_(1,4),u_(2,2),u_(2,3),u_(2,4),u_(3,3),u_(3,4),u_(4,4),
    s_(1,1),s_(1,2),s_(1,3),s_(1,4),s_(2,2),s_(2,3),s_(2,4),s_(3,3),s_(3,4),s_(4,4),
    x_1,x_2,x_3,x_4]
LT=matrix{{u_(1,1),0,0,0},{u_(1,2),u_(2,2),0,0},{u_(1,3),u_(2,3),u_(3,3),0},{u_(1,4),u_(2,4),u_(3,4),u_(4,4)}}
PSD=LT*transpose(LT)
subd=flatten toList apply(1..4,i->toList apply(i..4,j-> d_(i,j)=>PSD_(i-1,j-1)))
Ju=trim sub(sub(J,R2),subd);
netList Ju_*
equ=trim sub(sub(eq,R2),subd);
netList equ_*
PSD1=transpose(matrix{{x_1,x_2,x_3,x_4}})*matrix{{x_1,x_2,x_3,x_4}}
rank PSD1
subs=flatten toList apply(1..4,i->toList apply(i..4,j-> s_(i,j)=>PSD1_(i-1,j-1)))
equs=sub(equ,subs)
I=Ju+equs;
netList I_*
toString support I
Iu=time eliminate({x_1, x_2, x_3, x_4},I);
-- used  seconds
netList Iu_*
det sub(sub(D,R2),subs)  % gb Iu
Ru=QQ[d_(1,1),d_(1,2),d_(1,3),d_(1,4),d_(2,2),d_(2,3),d_(2,4),d_(3,3),d_(3,4),d_(4,4)]
Lu=time minimalPrimes sub(Iu,Ru);


-- GRAPH 13
restart
load "functions.m2"
R=QQ[l_1..l_6]
K=matrix{{l_1,l_3,0,l_6},{l_3,l_1,l_4,0},{0,l_4,l_2,l_5},{l_6,0,l_5,l_2}}
-- 2+2 consecutive vertices, no edges equal
(p,n,Rtotal,S)=coloredData(K)
suffStat(K)
stats=sub(suffStat(K),Rtotal);
netList stats_*

Raux=QQ[gens Rtotal,d_11,d_12,d_13,d_14,d_22,d_23,d_24,d_33,d_34,d_44]
gens Raux
D=genericSymmetricMatrix(Raux,d_11,4)

load "setupReciprocalVarieties.m2"
Dadj=adjugate D
f1=Dadj_(0,0)-Dadj_(1,1); --lambda1 
f2=Dadj_(2,2)-Dadj_(3,3); --lambda2

stats=sub(stats,Raux)
eq=trim sub(stats,{t_1=>d_11+d_22,t_2=>d_33+d_44,t_3=>2*d_12,t_4=>2*d_23,t_5=>2*d_34,t_6=>2*d_14})
netList eq_*

MLEeq=ideal{f1,f2}+eq
netList (trim MLEeq)_*

S=sub(S,Raux)
MLEeq1=MLEeq+minors(2,S);
MLEeq1d=time trim eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4),s_(2,2), s_(2,3),s_(2,4), s_(3,3),s_(3,4),s_(4,4)},MLEeq1);
-- used 0.125 seconds
netList MLEeq1d_*
codim MLEeq1d, degree MLEeq1d
--4,74
betti MLEeq1d
--2 cubics, 1 quartic, 1 quadric

L1d=time minimalPrimes MLEeq1d;
length L1d

