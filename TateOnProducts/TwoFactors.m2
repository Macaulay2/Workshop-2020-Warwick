restart
loadPackage "TateOnProducts"
loadPackage "Complexes"
kk = ZZ/101

relevantindex = {4,4}
topindex = {0,0}

totaltopindex = 0
totalrelevantindex = sum relevantindex

--Use 6 if using the old ChainComplex extend function
low = topindex
high = relevantindex

S = kk [x_0..x_1, y_0..y_1, Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
P1234 = ideal (x_0*x_1 - x_0 ^ 2 - x_1 ^ 2)
P5 = ideal random (S^1, S^{1:{-1,-1},1:{-1,-2}})

I = saturate intersect (P1234, P5)

--Set up all the relevant packages and rings. The constant relevantindex is the choice of index to use when calculating the action on global sections. Low and high are the bounds to use when calculating the cohomology table.

(R,E) = productOfProjectiveSpaces ({1,1}, CoefficientField=>kk)

I = sub (I, vars R)

ExactSequence = R^{1:{0,1}} ** res I

T1 = naiveTruncation (complex tateResolution (C = R^{1:{0,1}} ** module I,low,high), -(totalrelevantindex + 1), totaltopindex);
T2 = naiveTruncation (complex tateResolution (B=ExactSequence_1,low,high), -(totalrelevantindex + 1), totaltopindex);
T3 = naiveTruncation (complex tateResolution (A=ExactSequence_2,low,high), -(totalrelevantindex + 1), totaltopindex);

-- Find Tate resolutions of the modules involved.

N = netList {betti T3, betti T2, betti T1}

phi = ExactSequence.dd_2;
psi = ExactSequence.dd_1;
truncatedmap = inducedMap (C, truncate(relevantindex+{1,1},C));
inclusion = inducedMap (target psi, C);
psi = psi // inclusion;

-- This creates a collection of maps in the exact sequence, though the only relevant ones to use are phi and psi. The  issue is that the free resolution is not quite the same as the exact sequence we want.

d1 = psi * inducedMap (B, truncate (relevantindex+{1,1},B)) // truncatedmap;
--basis (relevantindex+1, psi)
--
--inducedMap (C, truncate(relevantindex+1,C))
d2 = phi * inducedMap (A, truncate (relevantindex+{1,1},A)) // inducedMap (B, truncate(relevantindex+{1,1},B));
--basis (relevantindex+1, phi)
-- We find the induced map on truncations, which in large degree is the induced map on global sections.

str = map (kk, R)

transfer1 = (str ** d1) ** E^{relevantindex+{1,1}};
transfer2 = (str ** d2) ** E^{relevantindex+{1,1}};

-- We need to coerce an integer matrix over R to one over E. This is a hugely silly way of doing that

e1 = extend (T1, T2, transfer1);
e2 = extend (T2, T3, transfer2);

--Now we extend the morphisms on global sections to one on the whole chain complex.

T1cover = cone e2;

coneQuasi = extend (T1, naiveTruncation (T1cover, -(totalrelevantindex + 1), -low_0), e1_-(totalrelevantindex + 1));

--The quasi-isomorphism between T1 and the cone over e2

coneInverseProjection = extend (T2, T1cover, id_(T2_-(totalrelevantindex + 1)));
coneInverseInclusion = extend (naiveTruncation (T1cover[1], -(totalrelevantindex + 1), totaltopindex), T3, gens (ker coneInverseProjection[1])_-(totalrelevantindex + 1))[-1];

--These two maps are the two maps which are not chain complex morphisms

prT1 = inducedMap (HH^(totalrelevantindex + 1) T1, T1_-(totalrelevantindex + 1));
prT1cover = inducedMap (HH^(totalrelevantindex + 1) T1cover, T1cover_-(totalrelevantindex + 1));

coneQuasiInverse = map (naiveTruncation (T1cover, -(totalrelevantindex + 1), -low_0), naiveTruncation (T1, -(totalrelevantindex + 1), -low_0) ,i -> (extend (naiveTruncation (T1cover, -(totalrelevantindex + 1), -low_0), T1, ((inverse (HH^(totalrelevantindex + 1) coneQuasi)) * prT1) // prT1cover))_i);

-- This is a quasi-isomorphism in middle degree, we lose some information when we use truncation. In fact even prune HH (id_T1) is not the action of the identity on cohomology because of truncation issues!

coneInclusion = map (naiveTruncation(T1cover, -(totalrelevantindex + 1), totaltopindex), naiveTruncation (T2,-(totalrelevantindex + 1),totaltopindex), i-> T1cover_i_(positions (toList(0..rank source coneInverseProjection_(i) -1), j -> not((coneInverseProjection_(i))_{j} == 0))));
coneProjection = map (T3[-1], naiveTruncation (T1cover, -(totalrelevantindex + 1), totaltopindex), i-> T1cover_i^(positions (toList(0..rank target coneInverseInclusion_(i) -1), j -> not((coneInverseInclusion_(i))^{j} == 0))));
boundaryMap = coneProjection * coneQuasiInverse;

--One should expect that e1 is the same as coneInclusion composed with coneQuasi, but in general it seems not to be. They agree on cohomology, so it's a problem in the extend method used to produce them. The most canonical objects to work with are e1, e2 and boundaryMap

(U,F) = productOfProjectiveSpaces ({1}, CoefficientField=>kk)

str = map (F, E,DegreeMap=> i -> {i_0})

isHomogeneous str

He1 = HH str(e1);
He2 = HH str(e2);
HboundaryMap = HH str(boundaryMap);

bettilist = {}

for i in (concentration source He1)_0..(concentration source He1)_1 do (
    bettilist = append (bettilist, HboundaryMap_i);
    bettilist = append (bettilist, He1_i);
    bettilist = append (bettilist, He2_i))

LES = (chainComplex bettilist);

betti LES

betti T3
betti T2
betti T1

subcomplexByDegrees = method ()

subcomplexByDegrees (ChainComplex, ZZ) := (C, d) -> (
    minC = min C;
    shiftC = C [-minC];
    return chainComplex apply (1 .. max C - minC, i -> submatrixByDegrees (shiftC.dd_i, d,d))[minC])

subcomplexByDegrees (ChainComplex, List) := (C, degs) -> (
    minC = min C;
    shiftC = C [-minC];
    return complex chainComplex apply (1 .. max C - minC, i -> submatrixByDegrees (shiftC.dd_i, degs,degs))[minC])

singlestrand = subcomplexByDegrees (LES, {0})[3*totalrelevantindex+6]

apply ({-10,-10}.. {10,10}, i-> if length subcomplexByDegrees (LES, i) == 0 then "" else (length subcomplexByDegrees (LES, i)+1, i))

-- The composite e1*e2 will not be zero, since it should know the boundary map, but its action on cohomology will vanish in middle degrees where we have truncated appropriately. The other composites however are zero. 
--coneQuasi and coneQuasiInverse are actually inverses to one another on the level of cohomology, it's not obvious since in a subquotient the identity map need not be the identity matrix.

LES

test = LES.dd_5

sourcepres = presentation source test
targetpres = presentation target test

symExt (sourcepres, U);
symExt (targetpres, U);


LES

HH^0 sheaf prune coker symExt (presentation HH^-6 LES, U)
HH^0 sheaf prune coker symExt (presentation HH^-1 LES, U)
HH^0 sheaf prune coker symExt (presentation HH^-2 LES, U)


betti LES

symExt (presentation source test, U)
symExt (presentation target test, U)

matrix test

test = prune symExt ( symExt (presentation l, E), R)

test

presentation R
