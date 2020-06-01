TOP_PRECRR = 128; -- implicit precision toRRp precision
-- Numbers
toQQ = (x,n) -> ((round(x*10^n))_QQ)/(10^n) -- convert RR number to QQ number with precision n digits after the decimal point
toRRp = (x) -> toRR(TOP_PRECRR,x) -- preset precision toRR
-- Linear algebra
mat = matrix -- matrix shortcut
inv = inverse -- matrix shortcut
rdim = A -> numRows(A) -- number of rows
cdim = A -> numColumns(A) -- number of columns
dims = A -> (rdim(A),cdim(A)) -- matrix size
rnum = A -> rdim(A)
cnum = A -> cdim(A)
ones = (m,n) -> trn(mat({toList(m:1)}))*mat({toList(n:1)}) -- matrix of ones
zeros = (m,n) ->0*ones(m,n) -- matrix of zeros
vo = A -> trn(matrix({flatten(entries(trn(A)))})) -- matrix vectorization
tra = A -> trace(A) -- trace
trn = A -> transpose(A) -- transposition shortcut
dia = x -> diagonalMatrix(x) -- diagonal matrix shortcut
eye = n -> dia(ones(n,1))
diag = A -> (n:=min(dims(A)); apply(toList(0..n-1),i->A_(i,i))) -- extract the diagonal from amatrix
wide = A -> (if rdim(A)>cdim(A) then trn(A) else A) -- make matrices wide
tall = A -> (if cdim(A)>rdim(A) then trn(A) else A) -- make matrices tall
fliplr = A -> A_(toList(reverse(0..cdim(A)-1))) -- flip matrix A horizontally
flipud = A -> A^(toList(reverse(0..rdim(A)-1))) -- flip matrix A verically
N2 = x -> first(ent(trn(x)*x)); -- squared euclidean norm
Real = A -> matf(A,a->realPart(a))
Imag = A -> matf(A,a->imaginaryPart(a))
CSum = A -> ones(1,rnum(A))*A
CAll = A -> matf(ones(1,rnum(A))*spy(A),a->if a==rnum(A) then 1 else 0)
Not  = A -> ones(rnum(A),cnum(A))-spy(A)
Find = A -> positions(ent(spy(A)),a->a!=0)
ReColsIx = (A,e) -> Find(matf(matcf(matf(Imag(A),a->abs(a)),c->max(c)),x->if x<e then 1 else 0))
ReCols = (A,e) -> Real(A_(ReColsIx(A,e))) 
PosColsIx = A -> Find(CAll(matf(A,a->if a>0 then 1 else 0)))
PosCols = A -> A_(PosColsIx(A))
xx = x -> (if (class(x)===Matrix) then mat{{0,-x_(2,0),x_(1,0)},{x_(2,0),0,-x_(0,0)},{-x_(1,0),x_(0,0),0}}
                                  else mat{{0,-x_(2),x_(1)},{x_(2),0,-x_(0)},{-x_(1),x_(0),0}}) -- 3x3 antisymmetric marix
ll = x -> mat{{-x_(1,2)},{x_(0,2)},{-x_(0,1)}} -- extract line coordinates from the skew symmetric matrix
elm  = (A,B) -> (local a; local b; m:=matrix({apply(flatten(entries(A)),flatten(entries(B)),(a,b)->a*b)}); trn(reshape((ring(m))^(cdim(A)),(ring(m))^(rdim(A)),m)) ) -- elemenetwise matrix division
matf = (A,f) -> ( m:=matrix({apply(flatten(entries(A)),f)}); trn(reshape((ring(m))^(cdim(A)),(ring(m))^(rdim(A)),m)) ) -- apply function to a matrix
matcf = (A,f) -> mat{apply(entries(trn(A)),f)} -- apply function f to folumns of f
mmf = (A,B,f) -> (local a; local b; m:=matrix({apply(flatten(entries(A)),flatten(entries(B)),(a,b)->f(a,b))}); trn(reshape((ring(m))^(cdim(A)),(ring(m))^(rdim(A)),m)) ) -- apply elementwise function to a pair of matrices
coln = A -> ones(1,rdim(A))*elm(A,A) -- matrix column norms
nul1 = (A) -> (local i; a:=wide(A); m:=trn(mat({apply(toList(0..cnum(a)-1),i->(-1)^i*det(submatrix'(a,,{i})))})); if rdim(A)>cdim(A) then trn(m) else m )
adj = A -> (local i; local j; matrix table(cnum(A),cnum(A),(i,j)->((-1)^(i+j))*det submatrix'(A,{j},{i}))) -- adjugate matrix
M2L = A -> (local i; apply(cnum(A),i->A_{i})) -- matrix to list of column vectors
L2M = L -> matrix({L}) -- list of column vectors to matrix
ent = A -> flatten(entries(trn(A))) -- list of elements of a matrix
Minors = (n,A) -> ent(gens(minors(n,A))) -- List of n x n minors of matrix A 
maxMinors = A -> ent(gens(minors(min(dims(A)),A))) -- list of maxial minors of matrix A 
maxMinorVec = A -> (X:=tall(A); trn(mat({apply(subsets(rdim(X),cdim(X)),i->det(X^i))}))) -- vectors of maximal minors
spy = A -> matf(A,x->(if x==0 then 0 else 1))
--- List manipulation
LIntersect = (A,B) -> (local x; select(A,x->member(x,B)))
LMinus = (A,B) -> (local x; select(A,x->not(member(x,B))))
LZip = (A,B) -> apply(A,B,(a,b)->{a,b})
mean = method()
mean(List) := A -> sum(A)/(#A)
mean(Matrix) := A -> mean(ent A)
var = method()
var(List) := A -> if #A==1 then 0 else sum(A-toList(#A:mean(A)),x->x^2)/(#A-1)
var(Matrix) := A -> var(ent A) 
std = A -> sqrt(var(A))
median = method()
median(List) := A -> (A = sort(A); if odd(#A) then i:={lift((#A-1)/2,ZZ)} else i={lift(#A/2-1,ZZ),lift(#A/2,ZZ)}; mean A_i)
median(Matrix) := A -> median(ent A)
--- Polyomial manipulation
eqnmaxabs = e -> apply(e,e->1/max(ent(gens(content(e)))/abs)*e) -- normalize coeffients of a list of polynomials to be in [-1,1]
primpart = f -> lcm((ent(gens(content(f))))/denominator)*f -- primitive part of a polynomial in ZZ

-- Solve 0-dim system of algebraic equations with QQ coefficients over CC
-- Solving polynomial equations via multiplicative matrix of a random pomynomial in R/I of a radical ideal I for a monomial basis b
-- returns {sols, vars, equation residuals, evaluated monomials, basis b}
SolvePolyEqByEigV = {DebugLevel => 0} >> opt -> X ->
(
    X = sequence(X); -- variable input must work even for single input
    I:={}; F:={};
    if isIdeal(X#0) then (I=X#0; F=ent gens I) else (F=X#0; I=ideal F); -- Ideal I and the list F of its generators
    R := ring(I); -- ring
    -- handle variable number of input arguments
    r := {}; if #X>3 then r = X#3;
    b := {}; if #X>2 then b = X#2;
    f := 0_R; if #X>1 then f = X#1;
    -- 
    if opt.DebugLevel>0 then print("Equations F = " | toString(F));
    if f==0 then (f = ((mat(randomMutableMatrix(1,cdim(vars(R)),0.0,100))*trn(vars(R)))_(0,0))_R + 1_R;) else (f = promote(f,R););-- random linear multiplication polynomial;
    if opt.DebugLevel>0 then print("Multiplier f = " | toString(f));
    -- f must not be constant
    if (support(f)=={}) then error("f must not be constant") else if opt.DebugLevel>0 then print("OK - f not constant");
    -- Compute the Groebner basis and install it to the list of known GBs
    -* try (J = IGeneratedByGBofI I;) 
    else (J = gens gb I;
          forceGB J; 
          J = ideal J);
    *-
    J := ideal gens gb I;
    -- Find the dimension and the degree of V(I) and check that the dimension=0
    dm := dim(J);
    dg := degree(J);
    if (dm=!=0) then error("dim I(F) must be 0") else if opt.DebugLevel>0 then print("OK - dim(I) = 0, deg(I) = " | toString(dg));
    mJ := rsort(unique(flatten(apply(ent(gens(J)),f->ent(monomials(f)))))); -- monmials of J
    s := {{},{},{}}; vr := {}; vs := {}; re := {}; 
    if all(flatten(mJ/degree),d->d<2) then ( -- linear system
	(B,A) := coefficients(gens(J),Monomials => mat({mJ})); A = trn(A); -- A x = 0 system
	a := -A_{cdim(A)-1}; -- right hand side
	A = A_{0..cdim(A)-2}; -- matrix of the system
	vs = solve(A,a); -- solve th elinear syste A x = a
	vr = mJ_{0..#mJ-2}; -- unknowns
	re = apply({ent(vs)},v->apply(F,f->sub(f,apply(vr,v,(b,v)->b=>v)))); -- equation residuals
	s = {ent(vs),vr,re,a2h(vs),trn(B)} -- costruct solutions	
    ) else (
        -- Factor ring
	A = R/J;
	use(R);
	-- The standard monomial basis of A
	B = sub(sort(basis(A),MonomialOrder=>Descending),R);
	if b=={} then (b = B;) else (b = mat({b}););
	if opt.DebugLevel>0 then print("Basis b = " | toString(b));
	-- Find thransformation from B to b (b may have more elements than B)
	C := (coefficients(b%J,Monomials => B))_1;
	if numcols B <= numcols b then (iC := trn(C)*inv(C*trn(C));) else (iC = inv(trn(C)*C)*trn(C););
	if eye(cdim(C))=!=sub(iC*C,ZZ) then (error("Bad change matrix from b to B");) else if opt.DebugLevel>0 then (print("OK - the change matrix from b to B is good"););
	-- Reduce f multiple of B by J
	if r=={} then (r = f*B;) else (r = mat({r}););
	if opt.DebugLevel>0 then print("Reduced pols r = " | toString(r));
	-- Multiplicaton matrix of y in R/I wrt B
	MfB := ((coefficients(r%I,Monomials => B))_1); -- wrt the standard monomial basis B
	Mf := iC*MfB*C; -- wrt the the given basis b
	if opt.DebugLevel>1 then print("Mf = " | toString(Mf));
	-- Check numerical solutions if computing in QQ
	if (coefficientRing(R)===QQ) then (
	    MR := matf(sub(Mf,QQ),x->toRRp(x)); -- to inexact real number
	    if opt.DebugLevel>1 then print("Mf_RR = " | toString(MR));
	    (e,v) := eigenvectors(trn(MR)); -- eigenvalues can be computer only in inexact complex numbers
	    if cdim(b)>1 and last(ent(b))==1 then ( -- normalize monomial vectors
		v = mmf(v,ones(cdim(v),1)*v^{rdim(v)-1},(x,y)->x/y); -- normalize to get evaluation of B in the eigenvectors 
		ix := positions(ent(b),m->(degree(m))#0==1); -- linear monomials
		if set((ent(b))_ix) === set(ent(vars(R))) then ( -- all variables present, then recover solutions
		    vix := v^ix; -- values of variables
		    vr = ent(trn(b_ix)); -- variables 
		    vs = trn(entries(vix)); -- solution list 
		    re = apply(vs,v->apply(F,f->sub(f,apply(vr,v,(b,v)->b=>v)))); -- equation residuals
		    s = {vs,vr,re}; -- solution structure
		);
	    );
	    s = join(s,{v,trn(b)}) -- return solutions and evauations on the basis
	) 
        else (error("Ring =!= QQ");)
    )
)
-- Gauss-Jordan elimination
gjelim = A -> (
    x := local x;
    T := ring(A)[x_1 .. x_(cdim(A))]; 
    X := trn(vars(T));
    G := ideal(sub(A,T)*X); 
    G  = mingens(gb(G)); 
    G  = sort(G,MonomialOrder=>Descending);
    (y,B) := coefficients(G,Monomials=>trn(X));
    sub(trn(B),ring(A)))
-- Remove fractions by row mutiplications
RowDeFract = A -> (
    r := M2L(trn(A)); -- list of row vectprs
    r  = apply(r,r->lcm((ent(r))/denominator)*r);
    trn(L2M(r))
    )
end ------------------

restart
-- loadPackage "NumericalAlgebraicGeometry"
needs "tomas-eigensolving.m2"
QQ[x,y]
I = ideal (x^2-1,y^3-1)
SolvePolyEqByEigV I
