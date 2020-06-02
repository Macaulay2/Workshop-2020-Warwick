-- random curve of genus 3
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
-- random canonical curve of genus 8
restart;
kk = ZZ
R = ZZ[x_0..x_7]

G26 = Grassmannian(1,5);
RG26 = ring G26;
mapp = map(R, RG26, random(R^1,R^{15:-1}, Height => 100)); 



