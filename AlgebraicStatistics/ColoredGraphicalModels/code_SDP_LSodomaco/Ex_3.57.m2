-- exercise 3.57
restart
R = QQ[q_1..q_55,gamma,x,y,z]

f = x^4+y^4+z^4-4*x*y*z+2*x+3*y+4*z -- want to find the global minimum of this polynomial function

-- PART ONE: SOS relaxation
X = symmetricPower(2,matrix{{1,x,y,z}}) -- f has degree 4=2*2, so we need to consider all monomials of degree <=2 in x,y,z
Q = genericSymmetricMatrix(R,q_1,10) 

F = f-gamma-(X*Q*transpose(X))_(0,0) -- we impose that f-gamma is SOS

-- now we impose that all coefficients of F are zero
X' = symmetricPower(4,matrix{{1,x,y,z}})
I = ideal first entries sub(contract(X', F), {x=>0,y=>0,z=>0})

L = apply(55, i-> q_(i+1) => q_(i+1)%I)

Q' = sub(Q, L) -- in our SDP program we want Q'>= 0

objFun = gamma
var = {gamma,q_11,q_12,q_13,q_20,q_21,q_22,q_23,q_24,q_28,q_29,q_30,q_31,q_32,q_33,q_41,q_42,q_46,q_47,q_48,q_53}
-- I have manually created the vector var with all variables appearing in Q'
M = Q'

needsPackage "SemidefiniteProgramming"
P = sdp(var, M, -objFun) -- this is for maximizing lambda
(X1,y1,Z1,stat) = optimize P
y1_(0,0) -- this is f_sos


-- PART TWO: check that f_sos is the global minimum
S = QQ[x,y,z]
f = x^4+y^4+z^4-4*x*y*z+2*x+3*y+4*z

-- we find all critical points of f
jacf = diff(matrix{{x,y,z}},f)
Igrad = ideal first entries jacf
needsPackage "NumericalAlgebraicGeometry"
sol = solveSystem first entries jacf
coordsol = apply(sol, i-> coordinates i)
#coordsol
for i in 0..#coordsol-1 do print (i,sub(f,{x=>(coordsol#i)#0,y=>(coordsol#i)#1,z=>(coordsol#i)#2}))

globalminimizer = coordsol#3 -- check the index, it might not be the same always
sub(f,{x=>globalminimizer#0,y=>globalminimizer#1,z=>globalminimizer#2}) -- this value coincides with f_sos

-- double check that the Hessian is positive definite, not necessary
Hessf = diff(transpose matrix{{x,y,z}},jacf)
eigenvalues sub(Hessf,{x=>globalminimizer#0,y=>globalminimizer#1,z=>globalminimizer#2}) 
-- hessian positive definite, so f is convex at the point, that is a local minimum.
