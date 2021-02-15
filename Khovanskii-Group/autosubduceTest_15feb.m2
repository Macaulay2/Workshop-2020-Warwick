restart
needsPackage "SubalgebraBases"
-- test 3 currently fails 
-- test 3 fails because 'sagbi' goes not return the correct generators
-- 'sagbi' fails because 'autosubduce' removes the term -1 from first generator
-- 'autosubduce' fails because 'reducedGens' has a+b+c without -1
-- 'reducedGens = {a+b+c, ...}' because 'subduction(notS, a+b+c-1) = a+b+c'
-- problem seems to be with 'rawSubduction'
 
kk = ZZ/101
R = kk[a,b,c]
F = matrix{{a+b+c-1, a^2+b^2+c^2-a, a^3+b^3+c^3-b}}
ans = matrix {{a+b+c-1, a*b+a*c+b*c+50*b+50*c, a*b*c+50*b^2+50*b*c+50*c^2-9*b+25*c}}

A = subring F
B = autosubduce A
peek B
--------

F12 = matrix{{F_(0,1), F_(0,2)}}
notS = subring F12
s = a+b+c-1
subduction(notS, s)

peek notS#"PresRing"

----------
-- checking subduction we that rawSubduction returns something strange:

-- INPUTS:
-- numblocks = 1
-- fMat = matrix{{p_0+p_1+p_2-1}}
-- F = map(kk[p_0..p_4],kk[p_0..p_4],{p_0, p_1, p_2, p_0^2+p_1^2+p_2^2-p_0, p_0^3+p_1^3+p_2^3-p_1})
-- syz_ideal = ideal(-p_0^2+p_3,-p_0^3+p_4)
-- J = gb syz_ideal
-- result = rawSubduction(numblock, raw fMat, raw F, raw J)

-- OUTPUT:
-- result == p_0+p_1+p_2 <-- equal to a+b+c when promoted
