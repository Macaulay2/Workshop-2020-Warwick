loadPackage "TateOnProducts"
kk = ZZ/101

S = kk [x_0..x_2]
P1234 = ideal (x_0*x_1 - x_0 * x_2, x_0 * x_2 -x_1 * x_2)
P5 = ideal random (S^1, S^{2:-1})

I = saturate intersect (P1234, P5)

--set up five points in P2

(R,E) = productOfProjectiveSpaces ({2}, CoefficientField=>kk)

I = sub (I, vars R)

--substitute them into the productofProjectiveSpaces format

ExactSequence = res I

--calculate a free resolution, from where we will get the exact sequence we will study

T1 = tateResolution (C = module I,{-10},{10})
T2 = tateResolution (B=ExactSequence_1,{-10},{10})
T3 = tateResolution (A=ExactSequence_2,{-10},{10})

--Find our three Tate resolutions

N = netList {betti T3, betti T2, betti T1}

--In case we want to inspect the cohomology table later

phi = ExactSequence.dd_2
psi = ExactSequence.dd_1

phi = ExactSequence.dd_2
psi = ExactSequence.dd_1
inclusion = inducedMap (target psi, C)
psi = psi // inclusion

--Obtain the morphisms in our exact sequence. Note that the differential ExactSequence.dd_1 is not quite right, 
-- it ends up inside the free hull of I, not inside I itself. So you need to lift along the inclusion.

d1 = psi * inducedMap (B, truncate (5,B)) // inducedMap (C, truncate(5,C))
d2 = phi * inducedMap (A, truncate (5,A)) // inducedMap (B, truncate(5,B))

-- These give two maps from free modules of the correct dimensions, with constant coefficients

source d2 == truncate(5,A)
target d2 == truncate(5,B)

source d1 == truncate(5,B)
target d1 == truncate(5,C)

d1*d2

T3_(-5)
matrix d1;
matrix d2;

--------------------------

phi = ExactSequence.dd_2
psi = ExactSequence.dd_1
inclusion = inducedMap (target psi, C)
psi = psi // inclusion

-- This creates a collection of maps in the exact sequence, though the only relevant ones to use are phi and psi. The  issue is that the free resolution is not quite the same as the exact sequence we want.

relevantindex=10
d1 = psi * inducedMap (B, truncate (relevantindex+1,B)) // inducedMap (C, truncate(relevantindex+1,C))
d2 = phi * inducedMap (A, truncate (relevantindex+1,A)) // inducedMap (B, truncate(relevantindex+1,B));

matrix d1;
matrix d2;
N
-- We find the induced map on truncations, which in large degree is the induced map on global sections.

str = map (kk, R)

transfer1 = (str ** d1) ** E^{relevantindex+1};
transfer2 = (str ** d2) ** E^{relevantindex+1};

-- We need to coerce an integer matrix over R to one over E. This is a hugely silly way of doing that
betti T1
betti T1[-relevantindex-1]
e1 = extend (T1[-relevantindex-1], T2[-relevantindex-1], transfer1);
e2 = extend (T2[-relevantindex-1], T3[-relevantindex-1], transfer2);

--Now we extend the morphisms on global sections to one on the whole chain complex.
T1shifted=T1[-relevantindex-1];
T2shifted=T2[-relevantindex-1];
T3shifted=T3[-relevantindex-1];

betti T1shifted
betti T2shifted
betti T3shifted

Nshifted=netList({betti T3shifted, betti T2shifted, betti T1shifted})

tempIndex=4;
assert(T1shifted.dd_(tempIndex+1) * e1_(tempIndex+1)==e1_tempIndex * T2shifted.dd_(tempIndex+1))
-- Caveat: not very sure it works properly in general; due to a basis choice issue

-- test function we may use:
needsPackage "ChainComplexExtras"
isChainComplexMap (e1)
isChainComplexMap (e2)


tempIndex2=9;
betti e1_tempIndex2
betti e2_tempIndex2

matrix{{0,0,0},{0,0,2},{2,0,0}}

betti e2_9
betti e1_9

betti e1_10
betti T1shifted.dd_10
betti e1_10

viewHelp Tor
viewHelp "ChainComplexExtras"
viewHelp HH

HH_8 e1
HH_9 e1
HH_10 e1

k=E^1/(ideal vars E)
e1k=e1**k;
HH_8 e1k
HH_9 e1k
HH_10 e1k

e2k=e2**k;
HH_8 e2k
HH_9 e2k
HH_10 e2k

source (HH_9 e2k)
target (HH_9 e2k)

residueMap = map(kk[t], E, apply(numgens E, i->E_i=>0))
betti residueMap (T1shifted)
betti residueMap T2shifted
betti residueMap T3shifted
reducede1=residueMap e1;
reducede2=residueMap e2;

HH_8 reducede1
HH_9 reducede1
HH_10 reducede1

HH_8 reducede2
HH_9 reducede2
HH_10 reducede2

kere1=(ker e1);
betti kere1
pruneKere1=prune kere1
pruneKere1_0
pruneKere1_1
pruneKere1_2

kere2=(ker e2);
pruneKere2=prune kere2

cokere1=(coker e1);
pruneCokere1=prune cokere1

prune ((ker e1)_10/(image e2)_10)
inducedMap((ker e1)_5,(image e2)_5)

betti cone e1
bConee1=prune HH beilinson(cone e1)

betti cone e2
bConee2=(prune HH beilinson(cone e2))

apply(toList(8..12), i->bConee2_i)
apply(toList(8..12), i->betti res bConee2_i)
betti I
saturate ann bConee2_10
ideal syz transpose presentation bConee2_11 == I

-- T3->T2->Cone(e2)->T3[1]

-- WANT: When we have a short exact seq. of chain complexes,
-- we want to compute explicit long exact sequences in homology.
-- In particular, we need to compute 'explicit' connecting homomorphisms.

-- Question : Two methods 'basis' and 'truncate' picks up bases, but they do not coincide in general.
-- Question : basis issue from Koszul complexes, basis of the universal bundle, ... sometimes such methods choose non-standard bases.
