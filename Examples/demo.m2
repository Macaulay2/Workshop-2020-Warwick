-- Proj-Uniform Parametrization
-- 06/05/2020
--
restart;
load "syzLLLtests.m2"
loadPackage("RandomSpaceCurves", Reload => true) 
printWidth=50

-- %%%%%%%%%%
-- Example 1:
-- %%%%%%%%%%
--
-- Construction of a curve of degree 7 and genus 2 in IP^3 
--
-- using RandomSpaceCurves package over QQ
time I = (random spaceCurve)(7,2,RQQ);
betti res I
height gens I
-- computation over ZZ is not possible
time IZZ = (random spaceCurve)(5,1,RZZ);
betti res IZZ
--
-- computation by hand 
h = 10;
-- random 3x10 matrix with linear entries 
M = random(RZZ^{3:0},RZZ^{10:-1},Height => h);
height M
-- we compare syzygies over QQ and ZZ
syzMQQ = syz sub(M,RQQ);
betti (syzM = syz M)
height syzMQQ, height syzM
-- 
-- construction of the curve 
-- we pick 8 linear syzygies of M
betti (N = syzM * random( source syzM, RZZ^{8:-2}, Height => h))
-- the entries of the syzygy matrix of N defines the curve
betti (syzN = syz N) 
I1 = sub(ideal syzN, RQQ);
-- a curve of degree 7 and genus 2
dim I1, degree I1, genus I1
log(height gens I1)/log(10)
-- 
-- we are doing the same computation using LLL
betti (syzMLLL = syzLLL1 M) 
height syzMLLL
betti (NLLL = syzMLLL * random( source syzMLLL, RZZ^{8:-2}, Height => h))
betti (syzNLLL = syzLLL1 NLLL) 
ILLL = sub(ideal matrixLLL syzNLLL,RQQ);
dim ILLL, degree ILLL, genus ILLL
-- comparing logarithmic height
(log(height gens I)/log(10), log(height gens I1)/log(h), log(height gens ILLL)/log(h))

-- %%%%%%%%%%
-- Example 2:
-- %%%%%%%%%%
--
-- Construction of a curve of degree 11 and genus 11 in IP^13 
--
-- using RandomSpaceCurves package over QQ
time J = (random spaceCurve)(11,11,RQQ);
betti res J
log(height gens J)/log(10)
height gens J
--
-- computation by hand 
h = 10;
M = random(RZZ^{3:0},RZZ^{10:-1},Height => h);
height M
betti (syzM = syz M)
betti (N = syzM * random( source syzM, RZZ^{7:-2,1:-3}, Height => h))
betti (syzN = syz N); 
J1 = sub(ideal syzN, RQQ);
dim J1, degree J1, genus J1
-- computation using LLL
betti (syzMLLL = syzLLL1 M) 
height syzMLLL
betti (NLLL = syzMLLL * random( source syzMLLL, RZZ^{8:-2}, Height => h))
betti (syzNLLL = syzLLL1 NLLL) 
JLLL = sub(ideal matrixLLL syzNLLL,RQQ);
dim JLLL, degree JLLL, genus JLLL
(log(height gens J)/log(10), log(height gens J1)/log(10), log(height gens JLLL)/log(10))

-- routine for computing several examples and comparing logarithmic height
-*
apply(10, i -> (
     J := (random spaceCurve)(11,11,RQQ);
     betti (N = syzM * random( source syzM, RZZ^{7:-2,1:-3}, Height => h));
     betti (syzN = syz N); 
     J1 = sub(ideal syzN, RQQ);
     betti (syzMLLL = syzLLL1 M); 
     betti (NLLL = syzMLLL * random( source syzMLLL, RZZ^{8:-2}, Height => h));
     betti (syzNLLL = syzLLL1 NLLL); 
     JLLL = sub(ideal matrixLLL syzNLLL,RQQ);
     print(log(height gens J)/log(10), log(height gens J1)/log(10), log(height gens JLLL)/log(10));
     ));
*-

-- %%%%%%%%%%
-- Example 3:
-- %%%%%%%%%%

-- Construction of a curve of degree 12 and genus 12 in IP^3 
--
-- Using the RandomSpaceCurves package does not finish after 45 min!!!
-- time C = (random spaceCurve)(12,12,RQQ);
--
-- Betti table of a curve of degree 12 and genus 12 
-- |total: 1 7 8 2|
-- |    0: 1 . . .|
-- |    1: . . . .|
-- |    2: . . . .|
-- |    3: . . . .|
-- |    4: . 7 5 .|
-- \    5: . . 3 2|

-- Computation by hand using LLL
betti (M12 = random(RZZ^{2:0},RZZ^{3:-1,5:-2}))
time betti (syzM12 = syzLLL1 M12)
betti (N12 = syzM12 * random( source syzM12, RZZ^{7:-3}))
time betti (syzN12 = syz(sub(N12,RQQ), DegreeLimit=>8)) -- used 135.296 seconds
C = ideal matrixLLL syzN12;
dim C, degree C, genus C
betti res C
log(height gens C)/log(10)
