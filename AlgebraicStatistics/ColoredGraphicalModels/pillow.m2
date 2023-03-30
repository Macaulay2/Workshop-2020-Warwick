restart

dualVariety=(IX,n,x,u)->(
    c:=codim IX;
    JacX:=diff(matrix{toList(x_1..x_n)}, transpose gens IX);
    AugJacX:=matrix{toList(u_1..u_n)}||JacX;
    SingX:=IX+minors(c,jacobian IX);
    conormalX:=saturate(IX+minors(c+1,AugJacX),SingX);
    dualX:=eliminate(toList(x_1..x_n),conormalX);
    dualX
    )

boundaryComponents=(K,p)->(
    I:=minors(p,K);
    minPrimes:=minimalPrimes I;
    m:= length minPrimes;
    allComponents:=for i to  m-1 list (dualComponent:=dualVariety(minPrimes_i,n,l,t),numgens trim dualComponent);
    boundaryComponents:= for i in allComponents list (if i_1==1 then i_0)
    )
algBoundary=(K)->(
    s=numgens target K;       
    delete (null, flatten (for p from 1 to s list boundaryComponents(K,p))) 
    )

n=3
R=QQ[l_1..l_n,t_1..t_n]
K=matrix{{1,l_1,0,l_1},{l_1,1,l_2,0},{0,l_2,1,l_3},{l_1,0,l_3,1}}
I3=minors(3,K)
gens gb I3
I4=minors(4,K)
gens gb I4

minPrimes3=minimalPrimes I3 
dualVariety(minPrimes3_0,n,l,t)
dualVariety(minPrimes3_1,n,l,t)
minPrimes4=minimalPrimes I4 
dualVariety(minPrimes4_0,n,l,t)
boundaryComponents(K,4)
boundaryComponents(K,3)
