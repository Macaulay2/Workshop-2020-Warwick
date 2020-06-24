-- Example 7.7 from Kaveh-Manon 

needsPackage "Tropical"

R = QQ[x1, x2, x3];
S = QQ[e1, e2, e3, y];

f = map(R, S, {x1 + x2 + x3, 
x1*x2 + x1*x3 + x2*x3,
x1*x2*x3,
(x1 - x2)*(x1 - x3)*(x2 - x3)});

I = kernel f;
toString mingens I

tropI = tropicalVariety I;

rays tropI
maxCones tropI

-- Write down a vector in the interior of each maximal cone
maxConeWeights = for cone in maxCones tropI list (
	for index from 0 to ambDim tropI - 1 list (
		val = 0;
		for rayIndex in cone do val = val + (rays tropI)_(index, rayIndex);
		val
	)
);

-- Note that leadTerm uses max convention so negate all weight vectors



maxConeIdeals = for w in maxConeWeights list (
	SOrdered = QQ[e1, e2, e3, y, Weights => -w,Global=>false];
	J = sub(I, SOrdered);
	inJ = leadTerm(1, J);
	print concatenate("Weight: ", toString w, 
		", Generators: ", toString inJ, 
		", isPrime: ", toString isPrime ideal inJ);
	J
);

needsPackage "SubalgebraBases"
A = subring matrix f
-- write fd as a polynomial in generators of A
T = QQ[gens R|gens S, MonomialOrder=> Eliminate 3]
G = sub(gens A,T)
I= ideal apply(drop(gens T, numgens R), toList(0..numgens A-1), (g, i) -> g - G_(0,i))
selectInSubring(1,gens gb I)
d=1
f = x1^(d+1)*x2^d+x2^(d+1)*x3^d+x3^(d+1)*x1^d
presentation A

gbIGens = gens gb I
fd % gbIGens -- fd = (fd%gbIGens) + (fd // gbIGens)*gbIGens

