-- A file intended for development only

pathToPackage = "./SubalgebraBases.m2"
installPackage(
    "SubalgebraBases",
    FileName=>pathToPackage,
    RunExamples => false,
    RerunExamples => false
    )
debug Core;


i = 2;
gndR = QQ[symbol t_1, symbol t_2, symbol t_3];
A := {t_1*t_2*t_3,
      t_1^2*t_2,
      t_1*t_2^2,
      t_1^2*t_3,
      t_1*t_3^2,
      t_2^2*t_3,
      t_2*t_3^2};
G := matrix {{t_1^(2*i)*t_2^(2*i)*t_3^(2*i), t_1^((3*i)+2)*t_2*t_3^(3*i)}}
subR = sagbi subring A;
assert((set first entries gens subR) === (set A)); 
tsyz := toricSyz(subR, G);
assert(tsyz * (transpose G) == 0);

ans1 = mingensSubring(subR, tsyz);
ans2 = extrinsicBuchberger(subR, tsyz);

assert(ans1 == ans2);

degs = matrix apply((entries (ans1//subR)), x -> apply(x, y -> first degree y));


print("---------------------");
print("result:");
print(ans1);
print("---------------------");
print("expected degrees:"|toString(i+1));
print("actual degrees:");
print(degs)
print("---------------------");
print("num. syzygies:")
print("Expected:"|toString( (2*i)+1 ));
print("  Actual:"|toString(numrows ans1));
print("---------------------");


error "stop";


-*
gndR = QQ[x_1..x_25];
M = genericMatrix(gndR, x_1, 5,5)
subR = subring gens minors(3,M);
ans = (det M )//subR
x = subR#"PresRing"#"FullSub"(ans)
error "stop"
*-
--setRandomSeed("randseed1");
moTest := (subR, n, maxDeg) -> (
    subMap := subR#"PresRing"#"Substitution";
    fullSub := subR#"PresRing"#"FullSub";
    incBase := subR#"PresRing"#"InclusionBase";

    print("-------------------------------------------"); 
    print("-- moTest:");

    for i from 1 to n do(
	
    	testElt :=  sum for deg from 0 to maxDeg list(
	    random(deg, ambient subR)
	    )//subR;
    	print("-------------------------------------------");
    	-- Possible ways to choose a lead term.
    	A := leadTerm testElt;	 
    	B := leadTerm(subR, testElt);
    	C := leadTerm subMap testElt;
    	D := leadTerm fullSub testElt;
	
    	print("LHS:");
    	print(leadTerm fullSub A);
    	print("RHS:");
    	print(leadTerm fullSub B);
	
	-- MO of the upper variables == MO induced by substitution map. 
    	--assert(leadTerm fullSub A == leadTerm fullSub B);
	
	-- If two terms of testElt have the same lead term under the substitution map
	-- and that term happens to be D, this will fail if they happen to cancel out.
	-- The chances of this happening seem to be very low, but we could cherry pick 
	-- examples using the function toricSyz.
	assert(leadTerm fullSub B == D);
    	-- MO of the lower variables == MO of the ambient ring.
    	assert(fullSub C == D);
    	);
    print("-------------------------------------------"); 
    print("-- End moTest"); 
    print("-------------------------------------------"); 
    );


-- The reason I use a finite field is because that's what the original programmers used
-- when they wrote the tests 20 years ago. I assume the reason they did this was so that
-- they wouldn't accidentally write code that relies on an infinite field.  
kk = ZZ/101

-- This function is from the tests, but here it is slightly modified.
genericminors = (minorsize,rowsize,colsize) -> (
    -- (minorsize) by (minorsize) minors of a generic (rowsize) by (colsize) matrix
    matdim := rowsize * colsize;
    x = symbol x;
    R := kk[x_1 .. x_matdim];
    Temp := genericMatrix(R,x_1,rowsize,colsize);
    print(Temp);
    gens minors(minorsize,Temp)
    );

print("toricSyz test (Sturmfels example 11.19)");
t = symbol t;
R = kk[symbol t_1, symbol t_2]
subR = sagbi subring matrix(R, {{t_1^2, t_1*t_2, t_2^2}});
M = matrix(R, {{t_1^2, t_1*t_2}});
M = subR#"PresRing"#"InclusionBase"(M);
ans = matrix(R,{{-t_2^2, t_1*t_2}, {-t_1*t_2, t_1^2}});
--ans = subR#"PresRing"#"InclusionBase"(ans);
time assert (toricSyz(subR, M) == ans);

--subR2 := value get "./test1";
--"./test1" << toString(subR) << close;


------------------------------------------
------------------------------------------
-- Sturmfels example 11.22
------------------------------------------
------------------------------------------
i = 2;
R = kk[symbol t_1, symbol t_2, symbol t_3];
A := {t_1*t_2*t_3,
      t_1^2*t_2,
      t_1*t_2^2,
      t_1^2*t_3,
      t_1*t_3^2,
      t_2^2*t_3,
      t_2*t_3^2};
G := matrix {{t_1^(2*i)*t_2^(2*i)*t_3^(2*i), t_1^((3*i)+2)*t_2*t_3^(3*i)}}
subR = sagbi subring A;
assert((set first entries gens subR) === (set A)); 
tsyz := toricSyz(subR, G);
assert(tsyz * (transpose G) == 0);

ans1 = mingensSubring(subR, tsyz);
ans2 = extrinsicBuchberger(subR, tsyz);

-- Subroutine 11.24 is guarenteed to return a (non-reduced) GB but subroutine 11.18 is only 
-- guarenteed to return some set of generators. So, it may be a coincidence that ans1 and ans2
-- coincide for all values of i that I have tried.
-- If tsyz generates syz(G), then it generates in(syz G) because in(syz G) = syz(G)?
assert(ans1 == ans2);

-*
print("i="|toString(i));
print("num syzygies:");
print(numrows ans1);
print("expected of degree "|toString(i+1)|":");
print((2*i)+2);
*-
-- the degree of each row is 3*(i+1) because all the generators have degree 3
degs := apply(entries (ans1//subR), row -> degree sum row)
assert(length (positions(degs, d -> first d == 3*(i+1))) == (2*i)+2);

------------------------------------------
------------------------------------------
-- Sturmfels example 11.25 
------------------------------------------
------------------------------------------

print("Sturmfels chapter 11 example 11.25.");
M = genericminors(2, 2, 5)

-- This monomial order is such that all the 2x2 determinants ab-cd satisfy ab > cd (I.e., the
-- positive term comes first.)
BaseRing := kk[x_1..x_10, MonomialOrder => {Weights=> {1,1,2,4,3,9,4,16,5,25}}]
N = sub(M, BaseRing);
subR := sagbi subring(N);
pres := subR#"PresRing";
tense = pres#"TensorRing";

-- The variables of the ambient ring are the variables of a generic matrix.
-- They are structured as follows, where the entries are subscripts of x,
-- the variable of the ambient ring:
--
-- | 1 3 5 7 9  |
-- | 2 4 6 8 10 |
--
-- The generators are the binomial(5,2)=10 2x2 determinants:
--
-- p_10 = [12]
-- p_11 = [13]
-- p_12 = [23]
-- p_13 = [14]
-- p_14 = [24]
-- p_15 = [34]
-- p_16 = [15]
-- p_17 = [25]
-- p_18 = [35]
-- p_19 = [45]

g1 = ((p_10-p_12)*(p_12-p_15))_tense;
g2 = ((p_11 + p_13 + p_16)*(p_16 + p_17 + p_18))_tense
G = matrix({{g1, g2}})
ltG = leadTerm pres#"Substitution" G;
-- g_3 in Sturmfels.
f = ((p_16*p_18*g1) - (p_12*p_15*g2))_tense

gensSyz := toricSyz(subR, ltG);

-- The result of toricSyz is really a KA module. 
KA = sagbi subring(leadTerm gens subR);
tenseKA = KA#"PresRing"#"TensorRing";

-- ans1 contains ans2. 
ans1 = extrinsicBuchberger(KA, gensSyz);
--ans2 = mingensSubring(KA, gensSyz);


error "stop";




dot = (v1, v2) -> (
    trace ((diagonalMatrix v1)*(diagonalMatrix v2))
    );
-- Returns a set of generators, not a Sagbi basis.
invariantsSO3 = vecList -> (
    subs1 := subsets(vecList, 2);
    dots := apply(subs1, S -> (
	    dot(first S, last S)
	    ));
    dotsDiag := apply(vecList, vec -> dot(vec, vec)); 
    dots = dots | dotsDiag;
    
    subs2 := subsets(vecList, 3);
    dets := apply(subs2, S -> (
	    det ((transpose S#0)|(transpose S#1)|(transpose S#2))
	    ));
    dots | dets 
    );

-*
R = QQ[(t_1..t_3)|(w_1..w_3)|(v_1..v_3), 
    MonomialOrder => {Eliminate 3, Lex}];
G = vars R;
eqns = invariantsSO3({G_{3..5}, G_{6..8}})
sag1 = sagbi eqns
t = transpose G_{0..2}
w = transpose G_{3..5}
v = transpose G_{6..8}
plucker := w||v
T = genericSkewMatrix(R, G_(0,0), 3);
zed = T-T;  
I = id_(source zed)
translation := matrix({{I, zed},{T, I}})*plucker
sag2 = sagbi transpose translation;
ans = sag2#"PresRing"#"FullSub"((gens sag1)//sag2)
ans = selectInSubring(1, gens sagbi ans)

print("Translation subaction:");
debugPrintGens(sag2);
print("Rotation subaction:");
debugPrintGens(sag1);
print("Intersection:");
error "stop";
*-


------------------------------------------------------------------------------------
-- Screw invariants
------------------------------------------------------------------------------------
-- Notice how we need an elimination order on the AMBIENT ring. 
gndR = QQ[(t_1..t_3)
         |(w_11,w_12,w_13)|(w_21,w_22,w_23)|(w_31,w_32,w_33)
         |(v_11,v_12,v_13)|(v_21,v_22,v_23)|(v_31,v_32,v_33),
     MonomialOrder => {Eliminate 3, Lex}];

W = matrix({{w_11,w_12,w_13},{w_21,w_22,w_23},{w_31,w_32,w_33}})
V = matrix({{v_11,v_12,v_13},{v_21,v_22,v_23},{v_31,v_32,v_33}})
t1 = matrix({{t_1, t_2, t_3}})
w1 = W^{0} 
w2 = W^{1}
w3 = W^{2}
v1 = V^{0}
v2 = V^{1}
v3 = V^{2}

-- the plücker coordinates of the screws
p1 = transpose (w1|v1)
p2 = transpose (w2|v2)
p3 = transpose (w3|v3)

-- Translation subaction
T = genericSkewMatrix(gndR, t1_(0,0), 3);
I = id_(source T)
translation := matrix({{I, T-T},{T, I}})
G := (flatten entries (translation*(p1|p2|p3)));
--sag1 = sagbi(G, PrintLevel=>3, Limit=>6)

sag1 = sagbi(G, PrintLevel=>2, Limit=>9)
sag1 = sagbi selectInSubring(1, gens sag1)
error "stop";



-- Rotation subaction
rot = sagbi(invariantsSO3({w1,v1,w2,v2}), PrintLevel=>2);

-- The verbose mode reveals some information here that agrees with Theorem 5.1:
-- 1. The 4 S-polynomials of degree 3 correspond to the binomial(4,3) 3x3 determinants.
-- 2. The following code verifies that the matrix of scalar products has 20 unique 2x2 
-- minors, which correspond to the 20 degree 4 S-polynomials:
H = (w1||v1||w2||v2)
H = exteriorPower_2 (H*(transpose H))
print("-- number of 2x2 minors of the matrix of scalar products:");
print(numgens autosubduce(subring flatten entries H));


print("-- Translation invariants:");
print(transpose gens sag1)
print("-- Intersection with rotation invariants:");
-- This is exactly what it's supposed to be according to the paper.
-- (These are f_1, f_2 and f_3, which they refer to as "the middle three invariants")
print(transpose ((gens sag1)//rot))

mid := compress rot#"PresRing"#"FullSub"((gens sag1)//rot)
other := matrix {invariantsSO3({w1,w2})}
ans = mid|other

print("-- Invariants of the adjoint action of SO(3) acting on two screws:");
print(transpose ans);

error "stop";
--=====================================================--
























R =  transpose genericMatrix(gndR, first gens gndR, 3, 3)
A = (R*(transpose R))-I;
B := det(R) - 1;
IG := gb ideal flatten((entries A)|{B})
test = matrix({{R, R-R},{R-R, R}});

error "stop";
O1 = genericSkewMatrix(gndR, (W^{0})_(0,0), 3);
O2 = genericSkewMatrix(gndR, (V^{0})_(0,0), 3);

--M := (R*O*(transpose R))%IG
A1 = (R*O1*(transpose R))%IG
A2 = (R*O2*(transpose R))%IG

test = matrix({{R, R-R},{R-R, R}});

error "stop";


testW = transpose W^{0}

test1 =  transpose genericSkewMatrix(gndR, first gens gndR, 3)
test2 =  transpose genericSkewMatrix(gndR, (W_(0,0)), 3)

error "stop"

M = transpose genericMatrix(gndR, (vars gndR)_(0,0),3,3)
A = (M*(transpose M))-(id_(source M))
A = (flatten entries A);
B = (det M) - 1 
eqns = A|{B}
IG = ideal eqns
sag2 = sagbi eqns
error "stop";

M = t||(W^{1})||(V^{1})
A = (M*(transpose M))-(id_(source M))
B = (det M) - 1 
eqns = eqns | (flatten entries A)|{B}
sag2 = sagbi eqns

ans0 = (sag1#"PresRing"#"FullSub")((gens sag2)//sag1)
ans1 = selectInSubring(1, gens sagbi ans0)





error "stop"












-- Invariants of the rotation subgroup of SE(3):
M =  transpose genericMatrix(R, first gens R, 3, 3)
A = (M*(transpose M))-(id_(source M))
B = (det M) - 1 
eqns := (flatten entries A)|{B}
sag1 = sagbi eqns

-- Invariants of the translation subgroup of SE(3):
G = vars R;
t = transpose G_{0..2}
w = transpose G_{3..5}
v = transpose G_{6..8}
plucker := w||v
T = genericSkewMatrix(R, G_(0,0), 3);
zed = T-T;  
I = id_(source zed)
translation := matrix({{I, zed},{T, I}})*plucker
sag2 = sagbi transpose translation;

-- This normal form calculation is computing an intersection of subrings.
ans0 = sag2#"PresRing"#"FullSub"((gens sag1)//sag2)
-- The sagbi call removes the constants/zeros and verifies that it's a sagbi basis. 
-- selectInSubring computes an intersection with a polynomial ring. 
ans1 = selectInSubring(1, gens sagbi ans0)


error "stop";


gndR = QQ[(r_1..r_9)|(t_1..t_3)|(w_1..w_3)|(v_1..v_3), 
    MonomialOrder => {Eliminate 12, Lex}];
assert(t_3 > w_1);
r1 := transpose (vars gndR)_{0..2}
r2 := transpose (vars gndR)_{3..5}
r3 := transpose (vars gndR)_{6..8}
t0 := transpose (vars gndR)_{9..11}
w0 := transpose (vars gndR)_{12..14}
v0 := transpose (vars gndR)_{15..17}

R =  transpose genericMatrix(gndR, first gens gndR, 3, 3)
O = genericSkewMatrix(gndR, w0_(0,0), 3);

I := id_(source R);
Z := I - I;

-- Construct the ideal that defines the Lie group SO(3).
A = R*(transpose R)-I;
B := det(R) - 1;
IG := ideal flatten((entries A)|{B})

-- The adjoint representation of SO(3) on its Lie algebra.
-- The nicest form of this matrix (i.e., the one of 3.22) is obtained by using
-- GRevLex order.
M := (R*O*(transpose R))%IG
----------------------------------------------------------------

-- The adjoint representation of SE(3) on its Lie aglebra in terms of Plucker
-- coordinates. (3.25)
-- a2( (R,t), (w,v) ) = A2*plucker (see 3.4.2.)
-- where: (R,t) \in SE(3,R)    	   (we use the skew matrix T to stand in for t)
--     	  (w,v) \in se(3,R)

T = genericSkewMatrix(gndR, t0_(0,0), 3);
A2 = (matrix({{R, Z},{T*R, R}}))
plucker = w0||v0;

test = transpose ((A2*plucker)||(transpose gens IG));
test =(transpose (A2*plucker));


sag = sagbi(autosubduce subring test, PrintLevel=>1);

ans = transpose ((gens sag)%IG)
extractEntries(ans, transpose plucker)
error "stop";

-- Elements of se(3, R) are of the form (B, b) where B is 3x3 skew symmetric
-- and b is a 1x3 vector. Hence, it only requires 3 coordinates.


----------------------------------------------------------------





error "stop"


print("---------------------------------------------------");
print("---------------------------------------------------");
print("---------------------------------------------------");
print("---------------------------------------------------");
print("---------------------------------------------------");

N = 6
M = 12
--gndG = QQ[(x_1..x_N)|(a_1..a_M), MonomialOrder => {Eliminate N}];
--gndG = QQ[(x_1..x_N)|(a_1..a_M), MonomialOrder => {Lex}];
gndG = QQ[(w_1,w_2,w_3,v_1,v_2,v_3,r_1,r_2,r_3,r_4,r_5,r_6,r_7,r_8,r_9,t_1,t_2,t_3), 
    MonomialOrder => {Eliminate N, Lex}];

rotG = matrix entries transpose (vars gndG)_{N..(8+N)}
transG = matrix entries transpose (vars gndG)_{(9+N)..(11+N)}
R = transpose genericMatrix(gndG, rotG_(0,0), 3, 3)
T = genericSkewMatrix(gndG, transG_(0,0), 3);
zed = T - T;  
I = id_(source zed)
adjAction = matrix({{R, zed},{T*R, R}});

-- gensRot generates the ideal of the rotation subgroup of G.
ROrtho := (R*(transpose R)) - I;
RUnit := (det R) - 1;
gensRot := flatten ((entries ROrtho)|{RUnit})
-- The ideal that defines G.
IG = ideal gensRot

w0 := transpose (vars gndG)_{0..2}
v0 := transpose (vars gndG)_{3..5}
plucker := w0||v0

subR = sagbi transpose (adjAction*plucker)

M = (gens subR)%IG

print("screw variables:");
print(transpose  plucker);

print("translation:");
print(transpose  transG);
print("rotation:");
print(transpose rotG);

print("sagbi:")
print(transpose M)

-- Not sure about the significance of these.
confused0 = extractEntries(transpose gens subR, matrix({{w_1,w_2,w_3}}))
confused1 = extractEntries(transpose M, matrix({{w_1,w_2,w_3}}))
confused2 = extractEntries(transpose M, matrix({{v_1,v_2,v_3}}))




error "stop";

------------------------------------------------------------------------------------
------------------------------------------------------------------------------------

------------------------------------------
------------------------------------------
-- Example 1, Stillman and Tsai 
-- Falsely succeeds
------------------------------------------
------------------------------------------
gndR = kk[x,y, MonomialOrder => Lex]
I = ideal(x^2 - x*y)
Q = gndR/I
subR = sagbi(subring {x}, PrintLevel => 2)
elt = ((gens subR)_(0,0))^1

J = gens subR#"PresRing"#"SyzygyIdeal"

for j from 0 to (numcols J)-1 do(
    pj := J_(0,j);
    
    --error "stop";
    );

G = lift(gens subR, ambient ambient subR)
subR2 = subring G;
extra = subR2#"PresRing"#"FullSub"((presentation Q)//subR2)
newElts = sub(extra, Q)

subR = sagbi(subring ((gens subR)|newElts), PrintLevel => 0)

error "stop";


------------------------------------------
------------------------------------------
-- Normal form demo
------------------------------------------
------------------------------------------
 
f = f_tense;
print("-- f:");
print(f);
print("-- f % subR:");
output1 := f % subR;
print(output1);
print("-- f // subR:");
output2 := f // subR;
print(output2);

-- Notice how the flags are set when sagbi fails.
R = kk[x, y]
F = matrix{{x, x*y-y^2, x*y^2}}
subR = sagbi(F,Limit=>30) 
peek subR

-- Notice how the flags are set when sagbi succeeds.
R = kk[y, x]
F = matrix{{x, x*y-y^2, x*y^2}}
subR = sagbi(F) 
peek subR

error "stop"
