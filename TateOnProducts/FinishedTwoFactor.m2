loadPackage "TateOnProducts"
loadPackage "Complexes"
kk = ZZ/101

S = kk [x_0..x_2, y_0..y_1, Degrees=>{{1,0},{1,0},{1,0},{0,1},{0,1}}]
P1234 = ideal (x_0*x_1 - x_0*x_2, x_0*x_1 - x_1 *x_2)
P5 = ideal random (S^1, S^{1:{-2,-1}, 1:{-1,0}})

I = saturate intersect (P1234, P5)

--Set up all the relevant packages and rings. The constant relevantindex is the choice of index to use when calculating the action on global sections. Low and high are the bounds to use when calculating the cohomology table.

(R,E) = productOfProjectiveSpaces ({2,1}, CoefficientField=>kk)

J = module sub (I, vars R) ** R^{1:{0,0}}

ExactSequence = res J

A=ker gens J
B=source gens J
C=J

relevantindex = toList apply (0..1, i -> max ({(coarseMultigradedRegularity A)_i, (coarseMultigradedRegularity B)_i, (coarseMultigradedRegularity C)_i}))
topindex = {-1,-1}

totaltopindex = 1
totalrelevantindex = sum relevantindex

--Use 6 if using the old ChainComplex extend function
low = topindex
high = relevantindex

T1 = naiveTruncation (complex tateResolution (C,low,high), -(totalrelevantindex + 1), totaltopindex);
T2 = naiveTruncation (complex tateResolution (B,low,high), -(totalrelevantindex + 1), totaltopindex);
T3 = naiveTruncation (complex tateResolution (A,low,high), -(totalrelevantindex + 1), totaltopindex);

-- Find Tate resolutions of the modules involved.

N = netList {betti T3, betti T2, betti T1}

phi = inducedMap (B,A);
psi = gens C;

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

diagonalOnesMatrix=method()

diagonalOnesMatrix(Ring,ZZ,ZZ) := (E,a,b)->(diagonalMatrix(E,a, b, toList apply (1.. min({a,b}), i-> 1)))

coneInverseProjection = extend (T2, T1cover, id_(T2_-(totalrelevantindex + 1)));
--coneInverseInclusion = extend (naiveTruncation (T1cover[1], -(totalrelevantindex + 1), totaltopindex), T3, gens (ker coneInverseProjection[1])_-(totalrelevantindex + 1))[-1];

--These two maps are the two maps which are not chain complex morphisms

prT1 = inducedMap (HH^(totalrelevantindex + 1) T1, T1_-(totalrelevantindex + 1));
prT1cover = inducedMap (HH^(totalrelevantindex + 1) T1cover, T1cover_-(totalrelevantindex + 1));

quasiInverseMap = extend (naiveTruncation (T1cover, -(totalrelevantindex + 1), -low_0), T1, ((inverse (HH^(totalrelevantindex + 1) coneQuasi)) * prT1) // prT1cover);

coneQuasiInverse = map (naiveTruncation (T1cover, -(totalrelevantindex + 1), -low_0), naiveTruncation (T1, -(totalrelevantindex + 1), -low_0) ,i -> quasiInverseMap_i);

-- This is a quasi-isomorphism in middle degree, we lose some information when we use truncation. In fact even prune HH (id_T1) is not the action of the identity on cohomology because of truncation issues!

truncT1cover = naiveTruncation(T1cover, -(totalrelevantindex + 1), totaltopindex);
truncT2 = naiveTruncation (T2,-(totalrelevantindex + 1),totaltopindex);

--oldconeInclusion = map (truncT1cover, truncT2, i-> T1cover_i_(positions (toList(0..rank source coneInverseProjection_(i) -1), j -> not((coneInverseProjection_(i))_{j} == 0))));
--oldconeProjection = map (T3[-1], truncT1cover, i-> T1cover_i^(positions (toList(0..rank target coneInverseInclusion_(i) -1), j -> not((coneInverseInclusion_(i))^{j} == 0))));

coneInclusion = map (truncT1cover, truncT2, i-> map (truncT1cover_i, truncT2_i, transpose coneInverseProjection_i));
---newconeProjection = map (T3[-1], truncT1cover, i-> map ((T3[-1])_i, truncT1cover_i, transpose coneInverseInclusion_i));

coneProjection = complex map (chainComplex T3[-1], chainComplex truncT1cover, i-> diagonalOnesMatrix (E, rank (T3[-1])_i, rank (truncT1cover_i)));

--coneInclusion = map (truncT1cover, truncT2, i-> map (truncT1cover_i, truncT2_i, transpose coneInverseProjection_i));
--newconeProjection = map (T3[-1], truncT1cover, i-> map ((T3[-1])_i, truncT1cover_i, transpose coneInverseInclusion_i));

boundaryMap = coneProjection * coneQuasiInverse;

--One should expect that e1 is the same as coneInclusion composed with coneQuasi, but in general it seems not to be. They agree on cohomology, so it's a problem in the extend method used to produce them. The most canonical objects to work with are e1, e2 and boundaryMap

goodColumns=method()
goodColumns(Module,List,Sequence) := (F,c,IJK) -> (
    degF:=degrees F;
    --sensitive to signs
    --select(#degF,g-> goodDegree(degF_g,c,IJK))
    select(#degF,g-> goodDegree(-degF_g,c,IJK))
    )

goodDegree=method()
goodDegree(List,List,Sequence) := (d,c,IJK) -> (
    (I,J,K) := IJK;
    #select(I,i-> d_i<c_i)==#I and #select(J,j->d_j==c_j)==#J and #select(K,k->d_k>=c_k)==#K
    )

inBeilinsonWindow = method(Options => options beilinson)
-- inBeilinsonWindow(deg, E)
--   returns a boolean value.
--   returns (beilinson(E^{-deg}) != 0),
--    i.e. whether the 'deg' values are 'in range'.
inBeilinsonWindow(List, Ring) := opts -> (deg, E) -> (
    (t,v,n,varsList,irrList) := ringData E;
    if opts.BundleType === PrunedQuotient or opts.BundleType === QuotientBundle then (
    for i from 0 to #deg-1 do (
        if deg#i < 0 or deg#i > n#i
        then return false;
        );
    true
    )
    else if opts.BundleType === FreeBundle then (
	for i from 0 to #deg-1 do (
            if deg#i < 0 or deg#i > (n#i)+1
            then return false;
        );
    	true
    )
    else if opts.BundleType === SubBundle  then (
	for i from 0 to #deg-1 do (
            if deg#i < 0 or deg#i > (n#i)
            then return false;
        );
    	true
    )    
    else if opts.BundleType === MapsBetweenFreeBundles then (
	for i from 0 to #deg-1 do (
	    if deg#i < (-1) or deg#i > (n#i)+1
	    then return false;
	    );
	true
	)
    else if opts.BundleType === DummyQuotientBundle then(
	for i from 0 to #deg-1 do (
        if deg#i < (-1) or deg#i > n#i
        then return false;
        );
    true
	)
)

ringData = method()
ringData Ring := E -> if not E.?TateRingData then E.TateRingData = (
  differentDegrees := unique last degrees vars E;
  varsLists := apply(differentDegrees, deg -> select (gens E, x-> degree x == deg));
  t := #varsLists;
  irrList := apply(varsLists, L -> ideal(L));
  v := varsLists/(L->#L);
  n := apply(v, i-> i-1);
  (t,v,n,varsLists,irrList)
  ) else E.TateRingData

ringData Module := M -> ringData ring M


removeZeroTrailingTerms = method()
removeZeroTrailingTerms(ChainComplex) := W -> (
    E := ring W;
    mi := nonzeroMin W;
    ma := nonzeroMax W;
    W' := W[mi];
    if mi==ma then (return (chainComplex({map(E^0,W'_0,0),map(W'_0,E^0,0)}))[-mi+1]) else
    (chainComplex apply(toList(1..ma-mi),i->W'.dd_i))[-mi]
    )

projectionMapOnEs=method()
projectionMapOnEs(Module,List) := (M,I)->(
    S := ring M;
    t := #degree S_0;
    v := apply(unique degrees source vars S, d->
	#select(degrees source vars S,e->e==d));    
    if not all(I,i->0 <= i and i <= t-1) then error "expected a sublist of {0,..,t-1}";
--    J := select(toList(0..t-1),j-> not member(j,I));
    nI := apply(I,i-> v_i-1);
    kk := coefficientRing S;
    (SI,EI) := productOfProjectiveSpaces (nI,CoefficientField=>kk);
    a:= null;
    phi1:= matrix {flatten apply(t,i->
	if member(i,I) then ( 
	    l:=position(I,j->j==i);
	    a=sum(l,j->nI_j+1);
	    apply(v_i,k->EI_(a+k))
	    )
	else apply(v_i,k->0))};
    E := (tateData S)#Rings_1;
    map(EI,E,phi1)
    )


nonzeroMin = method()
nonzeroMin(ChainComplex) := C -> (
    --assert( not C==0);
    if C==0 then return min C;
    m:= min C;
    while C_m==0 do (m=m+1);
    m)


nonzeroMax = method()
nonzeroMax(ChainComplex) := C -> (
    --assert( not C==0);
    if C==0 then return max C;
    m:= max C;
    while C_m==0 do (m=m-1);
    m)

projdegrees = {1};

t = #degree E_0;
J = select(toList(0..t-1),j-> not member(j,projdegrees));
str = projectionMapOnEs(module E,J);
F = target str;

sT1 = strand (chainComplex T1, toList (t:0), projdegrees);
sT2 = strand (chainComplex T2, toList (t:0), projdegrees);
sT3 = strand (chainComplex T3, toList (t:0), projdegrees);

sT1W = removeZeroTrailingTerms beilinsonWindow sT1;
sT2W = removeZeroTrailingTerms beilinsonWindow sT2;
sT3W = removeZeroTrailingTerms beilinsonWindow sT3;

--- complex seems to mess up everything, using chainComplex instead

sT1toT1 = map (chainComplex T1, sT1, i -> T1_i_(goodColumns (T1_i, toList (t:0), ({},projdegrees,{}))));
T1tosT1 = map (sT1, chainComplex T1, i -> T1_i^(goodColumns (T1_i, toList (t:0), ({},projdegrees,{}))));
sT2toT2 = map (chainComplex T2, sT2, i -> T2_i_(goodColumns (T2_i, toList (t:0), ({},projdegrees,{}))));
T2tosT2 = map (sT2, chainComplex T2, i -> T2_i^(goodColumns (T2_i, toList (t:0), ({},projdegrees,{}))));
sT3toT3 = map (chainComplex T3, sT3, i -> T3_i_(goodColumns (T3_i, toList (t:0), ({},projdegrees,{}))));
T3tosT3 = map (sT3, chainComplex T3, i -> T3_i^(goodColumns (T3_i, toList (t:0), ({},projdegrees,{}))));

sT1tosT1W = map (sT1W, sT1, i -> sT1_i^(positions (degrees sT1_i, a -> inBeilinsonWindow (a, E))));
sT1WtosT1 = map (sT1, sT1W, i -> sT1_i_(positions (degrees sT1_i, a -> inBeilinsonWindow (a, E))));
sT2tosT2W = map (sT2W, sT2, i -> sT2_i^(positions (degrees sT2_i, a -> inBeilinsonWindow (a, E))));
sT2WtosT2 = map (sT2, sT2W, i -> sT2_i_(positions (degrees sT2_i, a -> inBeilinsonWindow (a, E))));
sT3tosT3W = map (sT3W, sT3, i -> sT3_i^(positions (degrees sT3_i, a -> inBeilinsonWindow (a, E))));
sT3WtosT3 = map (sT3, sT3W, i -> sT3_i_(positions (degrees sT3_i, a -> inBeilinsonWindow (a, E))));

se1 = sT1tosT1W * T1tosT1 * chainComplex e1 * sT2toT2 * sT2WtosT2;
se2 = sT2tosT2W * T2tosT2 * chainComplex e2 * sT3toT3 * sT3WtosT3;
sboundaryMap = sT3tosT3W * T3tosT3[-1] * chainComplex boundaryMap * sT1toT1 * sT1WtosT1;

mi = min sT1W; ma=max sT1W;
W1 = new ChainComplex;
W1.ring = F;
apply(toList(mi..ma),i-> W1_i = F^(-apply(degrees sT1W_i,d->d_J))); 
apply(toList(mi+1..ma),i->W1.dd_i = map(W1_(i-1),W1_i,str(sT1W.dd_i)));

mi = min sT2W; ma=max sT2W;
W2 = new ChainComplex;
W2.ring = F;
apply(toList(mi..ma),i-> W2_i = F^(-apply(degrees sT2W_i,d->d_J))); 
apply(toList(mi+1..ma),i->W2.dd_i = map(W2_(i-1),W2_i,str(sT2W.dd_i)));

mi = min sT3W; ma=max sT3W;
W3 = new ChainComplex;
W3.ring = F;
apply(toList(mi..ma),i-> W3_i = F^(-apply(degrees sT3W_i,d->d_J))); 
apply(toList(mi+1..ma),i->W3.dd_i = map(W3_(i-1),W3_i,str(sT3W.dd_i)));

We1 = map (W1, W2, i -> matrix str(se1)_i);
We2 = map (W2, W3, i -> matrix str(se2)_i);
WboundaryMap = map (W3[-1], W1, i -> matrix str(sboundaryMap)_i);

FR1 = beilinson W1;
FR2 = beilinson W2;
FR3 = beilinson W3;

FRe1 = map (FR1, FR2, i -> beilinson We1_i);
FRe2 = map (FR2, FR3, i -> beilinson We2_i);
FRboundaryMap = map (FR3[-1], FR1, i -> beilinson WboundaryMap_i);

bettilist = {};

for i in (min source FRe1)..(max source FRboundaryMap) do (
    bettilist = append (bettilist, (HH FRboundaryMap)_i);
    bettilist = append (bettilist, (HH FRe1)_i);
    bettilist = append (bettilist, (HH FRe2)_i));

LES = complex bettilist;

isExact LES

prune LES
