restart
loadPackage "TateOnProducts"
loadPackage "Complexes"
kk = ZZ/101

relevantindex = 7
--Use 6 if using the old ChainComplex extend function
low = {0}
high = {6}

S = kk [x_0..x_2]
P1234 = ideal (x_0*x_1 - x_0 * x_2, x_0 * x_2 -x_1 * x_2)
P5 = ideal random (S^1, S^{2:-1})

I = saturate intersect (P1234, P5)

--Set up all the relevant packages and rings. The constant relevantindex is the choice of index to use when calculating the action on global sections. Low and high are the bounds to use when calculating the cohomology table.

(R,E) = productOfProjectiveSpaces ({2}, CoefficientField=>kk)

I = sub (I, vars R)

ExactSequence = res I

T1 = tateResolution (C = module I,low,high)
T2 = tateResolution (B=ExactSequence_1,low,high)
T3 = tateResolution (A=ExactSequence_2,low,high)

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

e1 = extend (T1[-relevantindex], T2[-relevantindex], transfer1)[relevantindex]
e2 = extend (T2[-relevantindex], T3[-relevantindex], transfer2)[relevantindex]

--Now we extend the morphisms on global sections to one on the whole chain complex.

NewT1 = complex T1
NewT2 = complex T2
NewT3 = complex T3

Newe1 = extend (NewT1, NewT2, transfer1)
Newe2 = extend (NewT2, NewT3, transfer2)

debugLevel=1

assert isComplexMorphism Newe1
assert isComplexMorphism Newe2

liftedmap = (transfer1 * NewT2.dd_-6 // NewT1.dd_-6)

liftedmap = (transfer1 * NewT2.dd_-7 // NewT1.dd_-7)

f1 = (NewT1.dd_-7 * liftedmap)  -- _{0..3}^{0..3}
f2 = (transfer1 * NewT2.dd_-7)  -- _{0..3}^{0..3}

NewT2.dd_7

NewT1
NewT2

prune HH_-7 NewT1
prune HH_-7 NewT2

Newe1_-8 == transfer1
Newe1_-8 * dd^NewT2_-7
g1 = (Newe1_-8 * dd^NewT2_-7) // dd^NewT1_-7
(Newe1_-8 * dd^NewT2_-7) % dd^NewT1_-7

(Newe1_-8 * dd^NewT2_-7) - dd^NewT1_-7 * g1 -- nonzero...
