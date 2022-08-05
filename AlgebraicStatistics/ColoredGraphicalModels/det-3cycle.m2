restart
R=QQ[l1,l2,l3,l4,l5,l6]
K=matrix{{l1,l2,l3},{l2,l1,l4},{l3,l4,l5}} -- 2 vertices
K=matrix{{l1,l2,l3},{l2,l1,l4},{l3,l4,l1}} -- 3 vertices
K=matrix{{l1,l2,l2},{l2,l3,l4},{l2,l4,l5}} -- 2 edges
K=matrix{{l1,l2,l2},{l2,l3,l2},{l2,l2,l4}} -- 3 edges
K=matrix{{l1,0,0},{0,l1,0},{0,0,l1}} -- 3 vertices, no edges
K=matrix{{l1,0,0},{0,l2,0},{0,0,l1}} -- 2 vertices, no edges
K=matrix{{l1,0,0},{0,l2,0},{0,0,l3}} -- uncolored vertices, no edges
f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI

hesf= rank jacobian gradI
-- 2 vertices - codim 3, rank Hessian=5
-- 3 vertices - codim 3, rank Hessian=4
-- 2 edges - codim 3, rank Hessian=5
-- 3 edges - codim 3, rank Hessian=4
-- 3 vertices, no edges - codim 1
-- 2 vertices, no edges - codim 1
-- uncolored vertices, no edges - codim 2
--linear rank??


--K4
restart
R=QQ[l_1..l_10]
K=matrix{{l_1,l_2,l_3,l_4},{l_2,l_1,l_5,l_6},{l_3,l_5,l_7,l_8},{l_4,l_6,l_8,l_9}} -- 2 vertices
K=matrix{{l_1,l_2,l_3,l_4},{l_2,l_1,l_5,l_6},{l_3,l_5,l_1,l_8},{l_4,l_6,l_8,l_9}} -- 3 vertices
K=matrix{{l_1,l_2,l_2,l_4},{l_2,l_3,l_5,l_6},{l_2,l_5,l_7,l_8},{l_4,l_6,l_8,l_9}} -- 2 edges

f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI

hesf= rank jacobian gradI


--4-cycle
K=matrix{{l_1,l_2,0,l_4},{l_2,l_10,l_5,0},{0,l_5,l_7,l_8},{l_4,0,l_8,l_9}} -- uncolored
K=matrix{{l_1,l_2,0,l_4},{l_2,l_1,l_5,0},{0,l_5,l_7,l_8},{l_4,0,l_8,l_9}} -- 2 vertices
K=matrix{{l_1,l_2,0,l_2},{l_2,l_10,l_5,0},{0,l_5,l_7,l_8},{l_2,0,l_8,l_9}} -- 2 edges
f=det K
gradI= ideal flatten entries jacobian ideal f
codim gradI

hesf= rank jacobian gradI
