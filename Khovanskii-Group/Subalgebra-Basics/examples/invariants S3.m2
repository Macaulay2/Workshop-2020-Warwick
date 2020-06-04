invariantsSn = (n) -> (
    -- ring of invariants of S_n
    R = ZZ/101[x_0 .. x_(n-1)]; 
    map(R^1, n, (j,i) -> sum apply(toList(x_0 .. x_(n-1)), x->x^(i+1))))

F = invariantsSn 3
ans = matrix {{x_0+x_1+x_2, x_0*x_1+x_0*x_2+x_1*x_2, x_0*x_1*x_2}}

S = subring ans ;
sB = time subalgebraBasis S;


