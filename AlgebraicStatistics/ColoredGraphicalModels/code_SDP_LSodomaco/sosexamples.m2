restart
R = QQ[z_1..z_21,x_1,x_2]

f = (x_1^2-x_1*x_2+3)^2+(3-2*x_1-5*x_2^2)^2+25*(1-x_1^2+7*x_2)^2
f = x_1^4+x_2^4

Q = genericSymmetricMatrix(R,z_1,6)

X = symmetricPower(2,matrix{{1,x_1,x_2}})

F = f-(X*Q*transpose(X))_(0,0)

X' = symmetricPower(4,matrix{{1,x_1,x_2}})
I = ideal first entries sub(contract(X', F), {x_1=>0,x_2=>0})

L = apply(21, i-> z_(i+1) => z_(i+1)%I)

Q' = sub(Q, L)

objFun = 1_R
--var = toList(z_1..z_21)
var = {z_7,z_8,z_12,z_13,z_14,z_19}
M = Q'

needsPackage "SemidefiniteProgramming"
P = sdp(var, M, objFun)
(X1,y1,Z1,stat) = optimize P
L = apply(first entries transpose y1, i-> promote(i,QQ))

Q'' = sub(Q', apply(6, i-> (var#i)=>(L#i)))

f-(X*Q''*transpose(X))_(0,0) -- = 0












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
eq = first entries jacf

needsPackage "Bertini"
makeB'InputFile(
    storeBM2Files,
    AffVariableGroup=>{x,y,z},
    B'Polynomials=>eq
    )

runBertini(storeBM2Files)
readFile(storeBM2Files)

realSolutions = importSolutionsFile(storeBM2Files,NameSolutionsFile=>"real_finite_solutions")

valuesf = apply(realSolutions, s -> sub(f,{x=>s#0,y=>s#1,z=>s#2}))

min valuesf
minPos = minPosition(valuesf)

Hessf = diff(transpose matrix{{x,y,z}},jacf)
eigenvalues sub(Hessf,{x=>(realSolutions#minPos)#0,y=>(realSolutions#minPos)#1,z=>(realSolutions#minPos)#2}) 
-- hessian positive definite, so f is convex at the point, that is a local minimum.





-- exercise 3.132
restart
R = QQ[z_1..z_21,s1,s2,c_1..c_6,x,y]

f = x+y^3-2
g = 1-x^2-y^2

X = symmetricPower(2,matrix{{1,x,y}})
t = (X*transpose(matrix{{c_1..c_6}}))_(0,0)

Q = genericSymmetricMatrix(R,z_1,6)
s0 = (X*Q*transpose(X))_(0,0)

F = t*f+s0+s1*g1+s2*g2
degree(x_1,F),degree(x_2,F)

X' = symmetricPower(4,matrix{{1,x_1,x_2}})
I = ideal first entries sub(contract(X', F), {x_1=>0,x_2=>0})

G = first entries gens gb I

L = apply(21, i-> z_(i+1) => z_(i+1)%I)

Q' = sub(Q, L)

objFun = gamma
var = {c_1,c_2,c_3,c_4,c_5,c_6,s_1,s_2,z_1,z_7,z_8,z_12,z_13,z_14,z_19}
M = Q'

needsPackage "SemidefiniteProgramming"
P = sdp(var, M, -objFun) -- this is for maximizing lambda
(X1,y1,Z1,stat) = optimize P
y1_(0,0) -- this is f_sos
