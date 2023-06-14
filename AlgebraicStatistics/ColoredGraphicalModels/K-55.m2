restart
n=3
m=ceiling (n/2)
k=floor(n/2)
L=flatten for i from 1 to n list
 for j from i to n list (i,j)
ElimLL=flatten  for i from 1 to m-1 list
 for j from i+1 to  m list (i,j)
ElimLR=flatten  for i from m+1 to n-1 list
 for j from i+1 to  n list (i,j)
ElimL= ElimLL|ElimLR
AllVar=for i in L list s_i

R=QQ[AllVar] 
ElimVar= for i in ElimL list s_i
S = genericSymmetricMatrix(R,n)
Igcr=eliminate(ElimVar,minors(k+1,S)) --expected non-zero because two minors survive
Imlt=eliminate(ElimVar,minors(k,S))

restart
n=4
m=ceiling (n/2)
k=floor(n/2)
L=flatten for i from 1 to n list
 for j from i to n list (i,j)
ElimLL=flatten  for i from 1 to m-1 list
 for j from i+1 to  m list (i,j)
ElimLR=flatten  for i from m+1 to n-1 list
 for j from i+1 to  n list (i,j)
ElimL= ElimLL|ElimLR
AllVar=for i in L list s_i

R=QQ[AllVar] 
ElimVar= for i in ElimL list s_i
S = genericSymmetricMatrix(R,n)
Igcr=eliminate(ElimVar,minors(k+1,S))
Imlt=eliminate(ElimVar,minors(k,S))

restart
n=5
m=ceiling (n/2)
k=floor(n/2)
L=flatten for i from 1 to n list
 for j from i to n list (i,j)
ElimLL=flatten  for i from 1 to m-1 list
 for j from i+1 to  m list (i,j)
ElimLR=flatten  for i from m+1 to n-1 list
 for j from i+1 to  n list (i,j)
ElimL= ElimLL|ElimLR
AllVar=for i in L list s_i

R=QQ[AllVar] 
ElimVar= for i in ElimL list s_i
S = genericSymmetricMatrix(R,n)
Igcr=eliminate(ElimVar,minors(k+1,S)) 
Imlt=eliminate(ElimVar,minors(k,S))

restart
n=6
m=ceiling (n/2)
k=floor(n/2)
L=flatten for i from 1 to n list
 for j from i to n list (i,j)
ElimLL=flatten  for i from 1 to m-1 list
 for j from i+1 to  m list (i,j)
ElimLR=flatten  for i from m+1 to n-1 list
 for j from i+1 to  n list (i,j)
ElimL= ElimLL|ElimLR
AllVar=for i in L list s_i

R=QQ[AllVar] 
ElimVar= for i in ElimL list s_i
S = genericSymmetricMatrix(R,n)
Igcr=eliminate(ElimVar,minors(k+1,S))
Imlt=eliminate(ElimVar,minors(k,S))


restart
KK=ZZ/113;
R=KK[s_(0,0)..s_(0,9),s_(1,1)..s_(1,9),s_(2,2)..s_(2,9),s_(3,3)..s_(3,9),s_(4,4)..s_(4,9),s_(5,5)..s_(5,9),s_(6,6)..s_(6,9),s_(7,7)..s_(7,9),s_(8,8),s_(8,9),s_(9,9)]
gens R
X=matrix{{3,5,7,11,13,17,19,23,29,31},{37,41,43,47,53,59,61,67,71,73},{79,83,89,97,101,103,107,109,113,127},{131,137,139,149,151,157,163,167,173,179}}
Snum=transpose(X)*X
S=matrix{{s_(0,0),Snum_(0,1),Snum_(0,2),Snum_(0,3),Snum_(0,4),}}

