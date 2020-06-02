-*
Examples from the book Groebner Bases and Convex Polytopes by Bernd Sturmfels.
Chapter 11 "Canonical Subalgebra Bases"
*-

-- example 11.2
--complete?

-- example 11.7
restart
needsPackage "SubalgebraBases"
R= QQ[t_(1,1)..t_(3,3),MonomialOrder=>Lex]
T=transpose genericMatrix(R,3,3)
S32=sort subsets(3,2)
kf = matrix{flatten apply(S32, S1 -> (
	    apply(S32, S2 -> (
		    (i1, j1) := (first S1, last S1);
		    (i2, j2) := (first S2, last S2);
		    det T_{i1,j1}^{i2,j2}
		    )
    		)
	    )
	)}
A = subring kf
sbA = subalgebraBasis A
peek A.cache

-- no finite subalg basis
restart
R=QQ[x,y,MonomialOrder=>Lex]
(options R).MonomialOrder
needsPackage "SubalgebraBases"
methods subring
A = subring {x+y,x*y,x*y^2}
elapsedTime sB = subalgebraBasis(A, Limit=>30);
elapsedTime sB = subalgebraBasis(A, Limit=>30); -- this should take less time

-- example 11.9
restart
needsPackage "SubalgebraBases"
w={1,2,3,4,5,1,4,9,16,25}
R= QQ[t_(1,1)..t_(2,5),
    MonomialOrder=>{Weights=>w}
    ]
T=transpose genericMatrix(R,5,2)

kA25 = matrix{apply(sort subsets(5,2), S -> (
	(i, j) := (first S, last S);
	det T_{i,j}
	)
    )}
leadTerm kA25
exponentMatrix = M -> matrix(first \ exponents \ (flatten entries leadTerm M))
A = transpose exponentMatrix kA25
liftedWeightVector = flatten entries(matrix{w}*A)
S25={(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)}
apply(S25,(i,j)->i+j^2)-liftedWeightVector
S=QQ[apply(S25,s->X_(first s,last s)),MonomialOrder=>{Weights=> liftedWeightVector}]
phi = map(R,S,kA25)
netList flatten entries leadTerm(1,ker phi) -- lead terms of Plucker relations (equation 11.9)

-- example 11.19
R=QQ[t_1,t_2]
kA=matrix{{t_1^2,t_1*t_2,t_2^2}}
M = matrix{{t_1^2,t_1*t_2}}
S=QQ[gens R|{x_0,x_1,x_2},MonomialOrder => (getMO R) | {Eliminate 2}]
phi = map(R,S,vars R | kA)
kAS = S / ker phi
xus = sub(M, kAS)
l=selectInSubring(1,syz xus)
phi lift(l,S)

-- example 11.22
R=QQ[t_1,t_2,t_3]
kA=matrix{{t_1*t_2*t_3,t_1^2*t_2,t_1*t_2^2,t_1^2*t_3,t_2^2*t_3,t_2*t_3^2}}
S=QQ[gens R|toList(x_0..x_(numcols kA-1)),MonomialOrder => (getMO R) | {Eliminate numgens R}]
phi = map(R,S,vars R | kA)
IA=ideal selectInSubring(1,gens ker phi)
kAS = S / ker phi -- note: it is cached
i=8
M = matrix{{
	t_1^(2*i)*t_2^(2*i)*t_3^(2*i),
	t_1^(3*i+2)*t_2*t_3^(3*i)}}
xus = sub(M, kAS)
l=selectInSubring(1,syz xus)
syzM = phi lift(l,S) -- enough?
