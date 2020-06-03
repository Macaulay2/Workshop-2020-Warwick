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
