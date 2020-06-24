-- workflow for installing this version of the package
-- run this as < needs "install.m2" >  from within a running instance of M2
uninstallPackage "SubalgebraBases"
restart
-- !!! assumes your M2 session starts inside of ".../-2020-Warwick/Khovanskii-Group/"
pathToPackage = "./SubalgebraBases.m2"
installPackage(
    "SubalgebraBases",
    FileName=>pathToPackage,
    RemakeAllDocumentation => not instance(newDoc, Symbol),
    RerunExamples => true
    )
check "SubalgebraBases"
help SubalgebraBases
-- viewHelp SubalgebraBases -- help in browser

end--
uninstallPackage "SubalgebraBases"
restart
installPackage "SubalgebraBases"



-- run next two lines to remake documentation: skipping the first keeps the same documentation
restart
newDoc = 2020
needs "install.m2"

R = QQ[x1, x2, x3];
S = QQ[e1, e2, e3, y];
needsPackage "SubalgebraBases"
viewHelp "SubalgebraBases"
A = subring {
    x1 + x2 + x3, 
    x1*x2 + x1*x3 + x2*x3,
    x1*x2*x3,
    (x1 - x2)*(x1 - x3)*(x2 - x3)
    }
d=1
f = x1^(d+1)*x2^d+x2^(d+1)*x3^d+x3^(d+1)*x1^d
f//A
f%A
member(f, A)



--- Jun 24: testing rawMonoidNumberOfBlocks
restart
needsPackage "SubalgebraBases"
debug Core
R=QQ[a,b,c,MonomialOrder=>{Weights=>{1,0,0}}]
rawMonoidNumberOfBlocks raw monoid R
MO=(options R).MonomialOrder
S=QQ[gens R|{x,y,z},MonomialOrder=>MO|{Eliminate 3}]
rawMonoidNumberOfBlocks raw monoid S
MOS=(options S).MonomialOrder
peek S
