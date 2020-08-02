-- This file contains:
-- - useful debugging functions (debugPrintMap and debugPrintAllMaps) 
-- - A new example from Sturmfels chapter 11
-- - Usage of new features
compactMatrixForm = true;

pathToPackage = "./SubalgebraBases.m2"
installPackage(
    "SubalgebraBases",
    FileName=>pathToPackage,
    RunExamples => false,
    RerunExamples => false
    )
export {
    "debugPrintMap",
    "debugPrintAllMaps"
    }
debugPrintMap := f -> (
    a := gens source f;
    for i from 0 to (length a)-1 do(
	elt := f(a_i);
	print("maps "|toString(a_i)|" to "|toString(elt));
	);    
    );
-- This is by far the easiest way to find the map that you need.
debugPrintAllMaps := subR -> (
    pres := makePresRing(subR);
    print("--------------------------------");
    print("-- PresRing map info dump:");
    print("--------------------------------");
    print("-- ProjectionInclusion:");
    debugPrintMap (pres#"ProjectionInclusion");
    print("-- ProjectionBase:");
    debugPrintMap (pres#"ProjectionBase");
    print("-- InclusionBase:");
    debugPrintMap (pres#"InclusionBase");
    print("-- Substitution:");
    debugPrintMap (pres#"Substitution");
    print("-- FullSub:");
    debugPrintMap (pres#"FullSub");
    print("-- ProjectionPres:");
    debugPrintMap (pres#"ProjectionPres");
    print("-- InclusionPres:");
    debugPrintMap (pres#"InclusionPres");
    print("-- SubstitutionPres:");
    debugPrintMap (pres#"SubstitutionPres");
    print("--------------------------------");
    print("-- End PresRing map info dump.");
    print("--------------------------------");
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
    --print(Temp);
    gens minors(minorsize,Temp)
    );


print("toricSyz test (Sturmfels example 11.19)");
t = symbol t;
R = kk[symbol t_1, symbol t_2]
subR = subalgebraBasis subring matrix(R, {{t_1^2, t_1*t_2, t_2^2}});
M = matrix(R, {{t_1^2, t_1*t_2}});
M = subR#"PresRing"#"InclusionBase"(M);
ans = matrix(R,{{-t_2^2, t_1*t_2}, {-t_1*t_2, t_1^2}});
--ans = subR#"PresRing"#"InclusionBase"(ans);
time assert (toricSyz(subR, M) == ans);


-- Sturmfels example 11.22
i = 2;
R = kk[symbol t_1, symbol t_2, symbol t_3];
A := {t_1*t_2*t_3, t_1^2*t_2, t_1*t_2^2, t_1^2*t_3, t_1*t_3^2, t_2^2*t_3, t_2*t_3^2};
B := matrix {{t_1^(2*i)*t_2^(2*i)*t_3^(2*i), t_1^((3*i)+2)*t_2*t_3^(3*i)}}
subR = subalgebraBasis subring A;
assert((set first entries gens subR) === (set A)); 
-- The algorithm was never guarenteed to generate a minimal set of generators.
-- In this case, the generators are redundant.
result := transpose toricSyz(subR, B)

tense := subR#"PresRing"#"TensorRing";
subst := subR#"PresRing"#"Substitution";
proj := subR#"PresRing"#"ProjectionBase";
pres := subR#"PresRing"#"PresentRing";

result2 := gens gb result;

debugPrintAllMaps subR;
error "stop";






print("Sturmfels chapter 11 example 11.25.");
M = genericminors(2, 2, 5)
-- (I don't understand where that weight vector comes from.)
BaseRing := kk[x_1..x_10, MonomialOrder => {Weights=> {1,1,2,4,3,9,4,16,5,25}}]
N = sub(M, BaseRing);
subR := subring(N);
pres := subR#"PresRing";
tense = pres#"TensorRing";
g1 = ((p_10-p_12)*(p_12-p_15))_tense;
g2 = ((p_15 + p_13 + p_16)*(p_16 + p_17 + p_18))_tense
G = matrix({{g1, g2}})
-- These agree with what Sturmfels says they should be
print("-- lead terms of G:");
print(leadTerm G)
-- g_3 in Sturmfels.
f = ((p_16*p_18*g1) - (p_12*p_15*g2))_tense

subRSagbi = subalgebraBasis subR;
presSagbi = subRSagbi#"PresRing";
tenseSagbi = presSagbi#"TensorRing";

assert(ambient subRSagbi === ambient subR);
assert(tenseSagbi =!= tense);
assert(subR == subRSagbi);
assert(G%subRSagbi == 0);



------------------------------------------
------------------------------------------
-- intrinsicBuchberger test
------------------------------------------
------------------------------------------
m1 := map(tenseSagbi, tense, gens tenseSagbi);
G = m1 G;

result = intrinsicBuchberger(subRSagbi, G)

result = result // subRSagbi;
h = intrinsicReduce(subRSagbi, result, f);

ltH = leadTerm(subRSagbi,h);
ltResult = leadTerm(subRSagbi, result);
print("is h zero? - "|toString(h == 0));

------------------------------------------
------------------------------------------
-- end intrinsicBuchberger test
------------------------------------------
------------------------------------------


-- Normal form demo
 
f = f_tense;
print("-- f:");
print(f);
print("-- f % subR:");
output1 := f % subR;
print(output1);
print("-- f // subR:");
output2 := f // subR;
print(output2);

-- Notice how the flags are set when subalgebraBasis fails.
R = kk[x, y]
F = matrix{{x, x*y-y^2, x*y^2}}
subR = subalgebraBasis(F,Limit=>30) 
peek subR

-- Notice how the flags are set when subalgebraBasis succeeds.
R = kk[y, x]
F = matrix{{x, x*y-y^2, x*y^2}}
subR = subalgebraBasis(F) 
peek subR

error "stop"
