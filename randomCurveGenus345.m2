

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 3
restart;
--kk = QQ;
kk = ZZ;
R = kk[x_0..x_2];
randomQuartic = random(4,R);
randomQuartic

-- random number with height
tally apply(100,i->random(ZZ,Height=>100))
tally apply(100,i->random(ZZ)
h = 100

-- random plane quartic with coefficients of height h
randomQuarticHeight = sum apply(flatten entries matrix basis(4,R), b -> ((random(ZZ,Height=>2*h)-h)*b))
    
randomQuarticHeight = (sub (random(kk^1,kk^15,Height=>h),R)*transpose matrix basis(4,R))_(0,0)



-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 4
restart;
kk = ZZ;
R = kk[x_0..x_3]
-- up to PGL(3), we choose the quadric
quadric = ideal( x_0*x_3-x_1*x_2 )
singQuadric = ideal (x_0^2-x_1*x_2)

-- Te be continued

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 5


-- We construct canonically embedded curve 
-- using Mukai's description

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 6
restart;
kk = ZZ;
R = ZZ[x_0..x_5];
RQ = QQ[x_0..x_5]

-- random choice 
--M' = random(R^5, R^{5:-1}, Height=>100)
--M = M' - transpose M'
-- 
M' = matrix{{0,x_0,x_1,x_2,x_3},
    	    {0,0,x_4,x_5,x_0+x_1},
	    {0,0,0,x_1+x_2,x_2+x_4},
	    {0,0,0,0,x_3+x_5},
	    {0,0,0,0,0}}
M = M' - transpose M'
-- smooth del Pezzo of degree 5
delPezzo = pfaffians(4,M) 
--dim sub(saturate ideal singularLocus delPezzo,RQ)

-- random quadric hypersurface
quadricHyp = sub(ideal random(2, R, Height => 100),RQ)

canCurve = quadricHyp + delPezzo;

dim canCurve, degree canCurve, genus canCurve, dim singularLocus canCurve

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 7
restart;
kk = QQ;
-- construction of the orthogonal Grassmannian in IP^15
R = kk[x_0,x_12..x_15,x_23..x_25,x_34,x_35,x_45,x_2345,x_1345,x_1245,x_1235,x_1234];
M = genericSkewMatrix(R,x_12,5)
v = transpose matrix{{x_2345,x_1345,x_1245,x_1235,x_1234}};
M'=ideal(M*v);
-- betti(resM'=res M')
-- betti(resM'4=resM'.dd_4)
OG=ideal(resM'4*transpose matrix{{x_0,1}});
betti res OG, dim OG, degree OG


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 8
restart;
kk = QQ
R = kk[x_0..x_7]

G26 = Grassmannian(1,5,CoefficientRing=>QQ);
RG26 = ring G26;
projToP7 = map(RG26, R, random(RG26^1,RG26^{8:-1}, Height => 100)); 
linSectionG26 = G26 + ideal random(RG26^1,RG26^{7:-1}, Height => 100);
dim linSectionG26, degree linSectionG26, genus linSectionG26
time canCurve = preimage_projToP7(linSectionG26);
dim canCurve, degree canCurve, genus canCurve

-- check smoothness over a finite field
S = ZZ/12343[flatten entries vars R];
canCurveS = sub(canCurve,S);
dim canCurveS, degree canCurveS, genus canCurveS
projToP3 = map(S, (coefficientRing S)[z_0..z_3], random(S^1,S^{4:-1}));
time spaceCurve = preimage_projToP3(canCurveS);
dim singularLocus spaceCurve

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 9
restart;
kk = QQ;
S = kk[x_0..x_5,y_0..y_5,u,z];
Mx = genericSymmetricMatrix(S,x_0,3);
My = genericSymmetricMatrix(S,y_0,3);

-- defining equations of LG(3,6)
LG = ideal(u*matrix{toList(y_0,-y_1,y_2,y_3,-y_4,y_5)}+mingens minors(2,Mx)) +
     ideal(u*z*matrix id_(S^3)-Mx*My) +
     ideal(z*matrix{toList(x_0,-x_1,x_2,x_3,-x_4,x_5)}+mingens minors(2,My));
betti res LG, dim LG, degree LG

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- random canonical curve of genus 10


restart;
--loadPackage "RandomSpaceCurves"
uninstallPackage "RandomSpaceCurves"
installPackage "RandomSpaceCurves"

loadPackage "RandomSpaceCurves"
viewHelp
-- space curve of degree 6 and genus 0 
syzLLL = method()
syzLLL (Matrix) := M ->(
    syzM := syz M;
    -- works only in a single degree
    -- one has to work degree by degree 
    (monM, coeffM) := coefficients syzM;
    coeffLLL := LLL lift(coeffM,ZZ);
    syzMLLL = map(target syzM, source syzM, monM * coeffLLL);
    -- TEST: Are these still syzygies?
    assert(M * syzMLLL == 0); 
    -- TEST: Is homogeneous?
    assert(isHomogeneous syzMLLL);
    syzMLLL)
betti coeffM

(g,d) = (0,6);
expectedBetti(0,3,6)
R = ZZ[x_0..x_3]
M1 = random(R^3,R^{9:-1},Height=>3)
betti(syzM1 = syz M1)
betti(syzM1 = syzLLL(M1))
degrees syzM1
M2 = syzM1 * random(source syzM1, R^{6:-2,1:-3},Height=>3);
betti M2
tally apply(flatten entries M2, isHomogeneous)
isHomogeneous matrix entries M2 
isHomogeneous M2
betti(syzM2 = syz M2)
I = trim ideal syzM2;
coefficient I_0

RQ = QQ[gens R]
IQ = trim sub(I,RQ);
IQ_0

I = (random spaceCurve)(d,g,R);
numgens R
betti res I
(gens I)_{0}
tally apply(100, i-> random(QQ))
random(19)-10
methods random
random(R^2,R^{2:-1},Height=>100)

options random
about random

M = sub(random(R^2,R^{4:-1}, Height=>10),R)
minors(2,M)
syzM = syz M
ideal mingens ideal syzM == minors(2,M)
(mon, coeff) = coefficients syzM;
coeffLLL = LLL lift(coeff,ZZ)
M * mon * coeffLLL
max flatten entries coeffLLL
max flatten entries coeff
max flatten entries (coeffLLL_{0})
max flatten entries (coeff_{0})

resM = res coker M
resM.dd_2
