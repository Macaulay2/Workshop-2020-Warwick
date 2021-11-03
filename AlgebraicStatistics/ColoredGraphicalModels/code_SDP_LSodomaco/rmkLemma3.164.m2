-- The Horn matrix is a 5x5 copositive matrix which does not satisfy the P + N condition of Lemma 3.164.
-- We check it with an SDP program.
restart
R = QQ[z_1..z_15,n_1..n_10]

H = matrix{{1,-1,1,1,-1},
           {-1,1,-1,1,1},
	   {1,-1,1,-1,1},
	   {1,1,-1,1,-1},
	   {-1,1,1,-1,1}}

P = genericSymmetricMatrix(R,z_1,5)
N = matrix{{0,n_1,n_2,n_3,n_4},
    	   {n_1,0,n_5,n_6,n_7},
	   {n_2,n_5,0,n_8,n_9},
	   {n_3,n_6,n_8,0,n_10},
	   {n_4,n_7,n_9,n_10,0}}

I = ideal flatten entries (H-P-N)

Q = P++diagonalMatrix(toList(n_1..n_10))

L = {};
var = {};
for i in 1..15 do if (z_i%I)!=z_i then L = append(L, z_i => z_i%I) else var = append(var, z_i);
for i in 1..10 do if (n_i%I)!=n_i then L = append(L, n_i => n_i%I) else var = append(var, n_i);
L
var

Q' = sub(Q,L)

objFun = 1_R

needsPackage "SemidefiniteProgramming"
P = sdp(var, Q', objFun)
(XX,yy,zz,stat) = optimize P
-- dual infeasible, so no solution
