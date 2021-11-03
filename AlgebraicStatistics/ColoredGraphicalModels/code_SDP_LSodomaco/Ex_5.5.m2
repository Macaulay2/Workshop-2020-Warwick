-- EXERCISE 4 --
-- Maximize and minimize the linear function 13x + 17y + 23z over the spectrahedron S in Example 12.4.
-- n=3
restart
needsPackage "SemidefiniteProgramming"
-- This constructs an SDP in dual form

R = QQ[x,y,z]
b = sub(matrix{{2,2,2}},R)

-- INPUTS:
-- var, a list, of variables
var = gens R
-- M, a matrix, constraint matrix with affine-linear entries
M = matrix{{1,x,y},{x,1,z},{y,z,1}}
-- objFun, a ring element, linear function to be minimized
objFun = (b*transpose(basis(1,R)))_(0,0)

P1 = sdp(var, M, objFun) -- this for minimizing the linear functional
P2 = sdp(var, M, -objFun) -- this for maximizing

--X, an n×n matrix, primal variable (not available if Solver=>"M2")
--y, an m×1 matrix, dual variable
--Z, an n×n matrix, dual variable
(X1,y1,Z1,stat1) = optimize P1 -- (-.5,-.5,-.5)
(sub(b,RR)*y1)_(0,0) -- this is the minimum
(X2,y2,Z2,stat2) = optimize P2 -- (1,1,1)
(sub(b,RR)*y2)_(0,0) -- this is the maximum


-- n=4
restart
needsPackage "SemidefiniteProgramming"
-- This constructs an SDP in dual form

R = QQ[x_1..x_6]
b = sub(matrix{{2,2,2,2,2,2}},R)

-- INPUTS:
-- var, a list, of variables
var = gens R
-- M, a matrix, constraint matrix with affine-linear entries
M = matrix{{1,x_1,x_2,x_3},{x_1,1,x_4,x_5},{x_2,x_4,1,x_6},{x_3,x_5,x_6,1}}
-- objFun, a ring element, linear function to be minimized
objFun = (b*transpose(basis(1,R)))_(0,0)

P1 = sdp(var, M, objFun) -- this for minimizing the linear functional
P2 = sdp(var, M, -objFun) -- this for maximizing

--X, an n×n matrix, primal variable (not available if Solver=>"M2")
--y, an m×1 matrix, dual variable
--Z, an n×n matrix, dual variable
(X1,y1,Z1,stat1) = optimize P1 -- (-.33,...,-.33)
(sub(b,RR)*y1)_(0,0) -- this is the minimum
(X2,y2,Z2,stat2) = optimize P2 -- (1,...,1)
(sub(b,RR)*y2)_(0,0) -- this is the maximum
