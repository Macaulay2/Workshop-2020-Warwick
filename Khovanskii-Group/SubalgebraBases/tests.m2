-- 1) Subring tests
-- 2) SubalgebraBases tests

-- 1) Subring tests

TEST ///
R = QQ[x1, x2, x3];
S = QQ[e1, e2, e3, y];
f = map(R, S, {x1 + x2 + x3, x1*x2 + x1*x3 + x2*x3, x1*x2*x3,
(x1 - x2)*(x1 - x3)*(x2 - x3)});
A = subring matrix f;
pR = presentationRing A;
assert (presentation A == matrix {{pR_0^2*pR_1^2-4*pR_0^3*pR_2-4*pR_1^3+18*pR_0*pR_1*pR_2-27*pR_2^2-pR_3^2}})
///
