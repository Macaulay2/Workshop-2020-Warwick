needsPackage "SubalgebraBases"

-*
IN: a List L of algebra generators
OUT: true/false: do the generators form a sagbi basis (wrt order on ambient ring)
*-
isSagbi = method()
isSagbi Subring := S -> (
    presS := S#"PresRing";
    IA := presS#"SyzygyIdeal";
    GBIA := gens gb IA;
    monomialSyzygies := selectInSubring(1, GBIA);
    remainders := compress subduction(S, presS#"FullSub" monomialSyzygies);
    numcols remainders == 0
    )
isSagbi List := L -> isSagbi subring L
isSagbi Matrix := M -> isSagbi subring M
end--
restart
needs "isSagbi.m2"
R = QQ[x,y,z]
L = {y*(x-1), y*x^2, y*(x^3+x^2+x+1), y^2} 
isSagbi L
isSagbi matrix{L}
M = {x+y+z,x*y+x*z+y*z, x*y*z, (x-y)*(x-z)*(y-z)}
isSagbi M
N = {x+y+z,x*y+x*z+y*z, x*y*z}
isSagbi N
