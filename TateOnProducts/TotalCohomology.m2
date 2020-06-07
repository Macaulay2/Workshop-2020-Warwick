restart
loadPackage "TateOnProducts"
loadPackage "Complexes"
kk = ZZ/101

relevantindex = 5
--Use 6 if using the old ChainComplex extend function
low = {-6}
high = {relevantindex}

S = kk [x_0..x_2]
P1234 = ideal (x_0*x_1 - x_0 * x_2, x_0 * x_2 -x_1 * x_2)
P5 = ideal random (S^1, S^{2:-1})

I = saturate intersect (P1234, P5)

--Set up all the relevant packages and rings. The constant relevantindex is the choice of index to use when calculating the action on global sections. Low and high are the bounds to use when calculating the cohomology table.

(R,E) = productOfProjectiveSpaces ({2}, CoefficientField=>kk)

I = sub (I, vars R)

ExactSequence = res I

T1 = naiveTruncation (complex tateResolution (C = module I,low,high), -high_0, 10)
T2 = naiveTruncation (complex tateResolution (B=ExactSequence_1,low,high), -high_0, 10)
T3 = naiveTruncation (complex tateResolution (A=ExactSequence_2,low,high), -high_0, 10)

-- Find Tate resolutions of the modules involved.

N = netList {betti T3, betti T2, betti T1}

phi = ExactSequence.dd_2
psi = ExactSequence.dd_1
--inclusion = inducedMap (target psi, image psi)
truncatedmap = inducedMap (C, truncate(relevantindex+1,C))
inclusion = inducedMap (target psi, C)
psi = psi // inclusion

-- This creates a collection of maps in the exact sequence, though the only relevant ones to use are phi and psi. The  issue is that the free resolution is not quite the same as the exact sequence we want.

d1 = psi * inducedMap (B, truncate (relevantindex+1,B)) // truncatedmap
--basis (relevantindex+1, psi)
--
--inducedMap (C, truncate(relevantindex+1,C))
d2 = phi * inducedMap (A, truncate (relevantindex+1,A)) // inducedMap (B, truncate(relevantindex+1,B))
--basis (relevantindex+1, phi)
-- We find the induced map on truncations, which in large degree is the induced map on global sections.

str = map (kk, R)

transfer1 = (str ** d1) ** E^{relevantindex+1}
transfer2 = (str ** d2) ** E^{relevantindex+1}

-- We need to coerce an integer matrix over R to one over E. This is a hugely silly way of doing that

e1 = extend (T1, T2, transfer1)
e2 = extend (T2, T3, transfer2)

--Now we extend the morphisms on global sections to one on the whole chain complex.

T1cover = cone e2

coneQuasi = extend (T1, naiveTruncation (T1cover, -high_0, -low_0), e1_-relevantindex)

--The quasi-isomorphism between T1 and the cone over e2

coneInverseProjection = extend (T2, T1cover, id_(T2_-relevantindex))
coneInverseInclusion = extend (naiveTruncation (T1cover[1], -relevantindex, 10), T3, gens (ker coneInverseProjection[1])_-relevantindex)[-1]

--These two maps are the two maps which are not chain complex morphisms

prT1 = inducedMap (HH^5 T1, T1_-5)
prT1cover = inducedMap (HH^5 T1cover, T1cover_-5)

coneQuasiInverse = map (naiveTruncation (T1cover, -high_0, -low_0), naiveTruncation (T1, -high_0, -low_0) ,i -> (extend (naiveTruncation (T1cover, -high_0, -low_0), T1, ((inverse (HH^5 coneQuasi)) * prT1) // prT1cover))_i)

-- This is a quasi-isomorphism in middle degree, we lose some information when we use truncation. In fact even prune HH (id_T1) is not the action of the identity on cohomology because of truncation issues!

coneInclusion = map (naiveTruncation(T1cover, -relevantindex, 10), naiveTruncation (T2,-relevantindex,10), i-> T1cover_i_(positions (toList(0..rank source coneInverseProjection_(i) -1), j -> not((coneInverseProjection_(i))_{j} == 0))))
coneProjection = map (T3[-1], naiveTruncation (T1cover, -relevantindex, 10), i-> T1cover_i^(positions (toList(0..rank target coneInverseInclusion_(i) -1), j -> not((coneInverseInclusion_(i))^{j} == 0))))
boundaryMap = coneProjection * coneQuasiInverse

--One should expect that e1 is the same as coneInclusion composed with coneQuasi, but in general it seems not to be. They agree on cohomology, so it's a problem in the extend method used to produce them. The most canonical objects to work with are e1, e2 and boundaryMap

prune HH (e2[-1] * boundaryMap)
prune HH (e1 * e2)
prune HH (boundaryMap * e1)

-- The composite e1*e2 will not be zero, even on cohomology, since it should know the boundary map. The other composites however are zero.