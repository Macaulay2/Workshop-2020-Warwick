restart
-- 2 vertices in 3-cycle
R=QQ[sigma_11..sigma_33,s_11..s_33]
I=ideal(sigma_11+sigma_22-s_11-s_22,sigma_12-s_12,sigma_13-s_13,sigma_23-s_23,sigma_33-s_33,
    sigma_22*sigma_33-sigma_23^2-sigma_11*sigma_33+sigma_13^2)
S=matrix{{s_11,s_12,s_13},{s_12,s_22,s_23},{s_13,s_23,s_33}}
J=ideal(s_11*s_22-s_12^2,s_11*s_33-s_13^2,s_22*s_33-s_23^2, det S)
K=trim eliminate(toList(s_11,s_12,s_13,s_22,s_23,s_33),I+J)
minPrimesK=minimalPrimes K
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

restart
R1=QQ[sigma_11..sigma_33,x_1..x_9]
use R1
XGeneral=matrix{{x_1..x_3},{x_4..x_6},{x_7..x_9}}
S=transpose(XGeneral)*XGeneral
I=ideal(sigma_12+sigma_13+sigma_23-S_(0,1)-S_(0,2)-S_(1,2),
    sigma_11-S_(0,0),sigma_22-S_(1,1),sigma_33-S_(2,2),
    sigma_13*sigma_23-sigma_12*sigma_33-sigma_12*sigma_23+sigma_13*sigma_22,
    sigma_13*sigma_23-sigma_12*sigma_33-sigma_12*sigma_13+sigma_11*sigma_23)

J=minors(2,S)
K=trim eliminate (toList(x_1..x_9),I+J)

X=matrix{{x_1..x_3},{x_4..x_6}}
SHat=transpose(X)*X
Adj=exteriorPower (2,SHat)