restart
-- 2 vertices in 3-cycle
R=QQ[sigma_11..sigma_33,s_11..s_33]
I=ideal(sigma_11+sigma_22-s_11-s_22,sigma_12-s_12,sigma_13-s_13,sigma_23-s_23,sigma_33-s_33,
    sigma_22*sigma_33-sigma_23^2-sigma_11*sigma_33+sigma_13^2)
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}
J=ideal(s_11*s_22-s_12^2,s_11*s_33-s_13^2,s_22*s_33-s_23^2, det S)
K=trim eliminate(toList(s_11,s_12,s_13,s_22,s_23,s_33),I+J)
Sigma=genericSymmetricMatrix(R,sigma_11,3)
--ideal(det(Sigma))+K

J=minors(2,S)
J=ideal(det S)

-- ideal(d_23^2-d_22*d_33,d_13*d_23-d_12*d_33,d_13^2-d_11*d_33,d_12^2*d_
--      33-d_11*d_22*d_33)

restart
-- 2 edges  in 3-cycle
R=QQ[sigma_11..sigma_33,s_11..s_33]
I=ideal(sigma_12+sigma_13-s_12-s_13,sigma_11-s_11,sigma_22-s_22,sigma_23-s_23,sigma_33-s_33,
    sigma_13*sigma_23-sigma_12*sigma_33-sigma_12*sigma_23+sigma_13*sigma_22)
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}
J=ideal(s_11*s_22-s_12^2,s_11*s_33-s_13^2,s_22*s_33-s_23^2, det S)
K=trim eliminate(toList(s_11,s_12,s_13,s_22,s_23,s_33),I+J)
minPrimesK=minimalPrimes K
K=trim eliminate(toList(s_11,s_12,s_13,s_22,s_23,s_33),I+ideal(det S))

restart
-- 3 edges  in 3-cycle
R=QQ[sigma_11..sigma_33,s_11..s_33]
I=ideal(sigma_12+sigma_13+sigma_23-s_12-s_13-s_23,sigma_11-s_11,sigma_22-s_22,sigma_33-s_33,
    sigma_13*sigma_23-sigma_12*sigma_33-sigma_12*sigma_23+sigma_13*sigma_22,sigma_13*sigma_23-sigma_12*sigma_33-sigma_12*sigma_13+sigma_11*sigma_23)
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}
J=ideal(s_11*s_22-s_12^2,s_11*s_33-s_13^2,s_22*s_33-s_23^2, det S)
K=trim eliminate(toList(s_11,s_12,s_13,s_22,s_23,s_33),I+J)
minPrimesK=minimalPrimes K
K=trim eliminate(toList(s_11,s_12,s_13,s_22,s_23,s_33),I+ideal(det S))