restart
needsPackage "SubalgebraBases"

(c, d) = (6, 4) -- these can be anything
R=QQ[t,x_1..x_(c+d),MonomialOrder=>{Weights=>{1}|toList(0:(c+d)),Lex}]
M=matrix{
    toList(x_1..x_c),
    toList(x_(1+d)..x_(c+d))
    }
A=subring(drop(gens R,1) | apply(subsets(c,2), s -> t* det M_s))
gens A

subalgebraBasis(A, Limit=>10, PrintLevel=>1)
subalgebraBasis(A, Limit=>9, PrintLevel=>1)
