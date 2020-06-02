restart
loadPackage "TateOnProducts"

(S,E)=productOfProjectiveSpaces {2,1}
psi=random(E^{{-1,0}},E^{{-2,-1}})
phi=beilinson psi

beilinson(E^{{-1,0}})
beilinson(E^{{0,-1}})

T=chainComplex psi
C=beilinson T
prune HH C

restart
loadPackage "TateOnProducts"
(S,E)=productOfProjectiveSpaces{1,2}
M=S^1
low={-5,-5},high={5,5}

cohomologyMatrix(M,low,high)
T=tateResolution(M,low,high)
betti T
cohomologyMatrix(T,2*low,2*high)
betti T

