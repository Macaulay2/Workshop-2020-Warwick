--This file contains some of the commands used during the Macaulay2
--introduction

1+1
3!
binomial(4,2)
factor(60)

QQ
ZZ/3
ZZ/32749
GF(4)

RR -- you probably don't want this one!
CC -- or this one!



--Matrices

A = matrix {{1,1,1,1},{0,1,2,3}}

rank A

kernel A

gens kernel A

entries transpose gens kernel A

B = transpose(A)*A;

det(B)

rank(B)

R=QQ[x,y,z]

A=matrix {{1,x,x^2},{1,y,y^2},{1,z,z^2}}

entries A
flatten entries A

det(A)

factor(det(A))

help det

viewHelp det


--Ideals
R=QQ[x_0,x_1,x_2,x_3]

I=ideal(x_0*x_2-x_1^2, x_0*x_3-x_1*x_2, x_1*x_3-x_2^2)

gb I
gens gb I
leadTerm I

isPrime I

dim(I)

degree(I)

hilbertPolynomial(I,Projective=>false)


--Maps between rings

R=QQ[x,y]
S=QQ[t]
f = map(S,R,{t,t^2})

f(x^2+y^2)

kernel f

--First break 

--if/then
a = 7;
if a<8 then <<"a is less than 8"<<endl;   --  << is the print command

--loops (three versions)

scan(10,i->(<<i^2<<endl;))   -- loops start at 0

apply(10,i->i^2)

for i from 0 to 9 do
    <<i^2<<endl;

--Functions

--evaluate a polynomial

p = n ->(n^2+2*n+3);

p(3)
p(x)

--Test using "assert"

assert(p(1)==6)

assert(p(1)==5)
    
--Functions can take multiple inputs, and take up multiple lines

--Compute the minor of a matrix A indexed by the lists of indices in rowindex and columnindex

minorOfMatrix = (A,rowindex,columnindex)->(
    if #rowindex != #columnindex then error("The minor needs to be square");
    if min(rowindex)<0 or min(columnindex)<0 then error("The input range is wrong");
    if max(rowindex)> rank target A -1 or max(columnindex)> rank source A -1 then error("The input range is wrong");
    B:=A^rowindex_columnindex;
    return(det(B));   
);

minorOfMatrix = method();

minorOfMatrix(Matrix,List,List):= (A,rowindex,columnindex)->(
    if #rowindex != #columnindex then error("The minor needs to be square");
    if min(rowindex)<0 or min(columnindex)<0 then error("The input range is wrong");
    if max(rowindex)> rank target A -1 or max(columnindex)> rank source A -1 then error("The input range is wrong");
    B:=A^rowindex_columnindex;
    return(det(B));   
);


--Return all minors with a given row subset
minorOfMatrix(Matrix,List):=(A,rowindex)->(
    d:=#rowindex;
    m:=rank source A;
    apply(subsets(m,d),s->(minorOfMatrix(A,rowindex,s)))
);    
   
    
A=matrix{{1,1,1,1},{0,2,5,11}}    

minorOfMatrix(A,{0,1},{2,3})

minorOfMatrix(A,{0,1})   

viewHelp method



end

restart

load "M2examples.m2"