-- example 3.129
restart
R = QQ[z_1..z_21,s1,s2,c_1..c_6,x_1,x_2]

f = x_1^2+x_2^2-1
g1 = 3*x_2-x_1^3-2
g2 = x_1-8*x_2^3
-- want to certify that the system {f=0,g1>=0,g2>=0} has no solutions over RR^2,
-- using a Positivstellensatz certificate
-- need to find a polynomial t and sos polynomials s0,s1,s2 such that t*f+s0+s1*g1+s2*g2 = -1
-- we fix degree bound 4, in particular deg(t)=2, deg(s0)=4, deg(s1)=deg(s2)=0

X = symmetricPower(2,matrix{{1,x_1,x_2}})
t = (X*transpose(matrix{{c_1..c_6}}))_(0,0)

Q0 = genericSymmetricMatrix(R,z_1,6)
s0 = (X*Q0*transpose(X))_(0,0)

F = t*f+s0+s1*g1+s2*g2+1
degree(x_1,F),degree(x_2,F)

X' = symmetricPower(4,matrix{{1,x_1,x_2}})
I = ideal first entries sub(contract(X', F), {x_1=>0,x_2=>0})

L = {};
var = {};
for i in 1..21 do if (z_i%I)!=z_i then L = append(L, z_i => z_i%I) else var = append(var, z_i);
for i in 1..6 do if (c_i%I)!=c_i then L = append(L, c_i => c_i%I) else var = append(var, c_i);
if (s1%I)!=s1 then L = append(L, s1 => s1%I) else var = append(var, s1);
if (s2%I)!=s2 then L = append(L, s2 => s2%I) else var = append(var, s2);

Q0' = sub(Q0,L)
Q = Q0'++matrix{{s1%I}}++matrix{{s2%I}}

objFun = 1_R

needsPackage "SemidefiniteProgramming"
P = sdp(var, Q, objFun)
(XX,yy,zz,stat) = optimize P

prom = apply(first entries transpose yy, i-> promote(i,QQ))
subs = apply(#var, i-> (var#i)=>(prom#i))

Q0'' = sub(Q0', subs)
s0' = (X*Q0''*transpose(X))_(0,0)
s1' = sub(s1, subs)
s2' = sub(s2, subs)
t' = sub(t, subs)

t'*f+s0'+s1'*g1+s2'*g2 -- if this is equal to -1, then you're done!

needsPackage "SumsOfSquares"

sol1 = solveSOS s0'
sosPoly sol1






-- exercise 3.132
restart
R = QQ[z_1..z_21,w_1..w_6,c_0..c_2,x,y]

f = x+y^3-2
g = 1-x^2-y^2
-- want to certify that the system {f=0,g=0} has no solutions over RR^2.
-- Instead, we show that the system {f=0,g>=0} has no solutions over RR^2,
-- using a Positivstellensatz certificate
-- need to find a polynomial t and sos polynomials s0,s1,s2 such that t*f+s0+s1*g = -1
-- we fix degree bound 4, in particular deg(t)=1, deg(s0)=4, deg(s1)=2

t = c_0+c_1*x+c_2*y

X1 = symmetricPower(2,matrix{{1,x,y}})
X2 = matrix{{1,x,y}};

Q0 = genericSymmetricMatrix(R,z_1,6)
Q1 = genericSymmetricMatrix(R,w_1,3)
s0 = (X1*Q0*transpose(X1))_(0,0)
s1 = (X2*Q1*transpose(X2))_(0,0)

F = t*f+s0+s1*g+1
degree(x,F),degree(y,F)

X' = symmetricPower(4,matrix{{1,x,y}})
I = ideal first entries sub(contract(X', F), {x=>0,y=>0})

L = {};
var = {};
for i in 1..21 do if (z_i%I)!=z_i then L = append(L, z_i => z_i%I) else var = append(var, z_i);
for i in 1..6 do if (w_i%I)!=w_i then L = append(L, w_i => w_i%I) else var = append(var, w_i);
for i in 0..2 do if (c_i%I)!=c_i then L = append(L, c_i => c_i%I) else var = append(var, c_i);

Q0' = sub(Q0, L)
Q1' = sub(Q1, L)
Q = Q0'++Q1'

objFun = 1_R

needsPackage "SemidefiniteProgramming"
P = sdp(var, Q, objFun)
(xx,yy,zz,stat) = optimize P

prom = apply(first entries transpose yy, i-> promote(i,QQ))
subs = apply(#var, i-> (var#i)=>(prom#i))

Q0'' = sub(Q0', subs)
Q1'' = sub(Q1', subs)
s0' = (X1*Q0''*transpose(X1))_(0,0)
s1' = (X2*Q1''*transpose(X2))_(0,0)
t' = sub(t, subs)

t'*f+s0'+s1'*g -- if this is equal to -1, then you're done!

needsPackage "SumsOfSquares"

sol1 = solveSOS s0'
sosPoly sol1

sol2 = solveSOS s1'
sosPoly sol2
