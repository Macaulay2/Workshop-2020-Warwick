L=[l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12]

m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);

cvx_begin
    variable x(n)
    minimize( norm(A*x-b) )
cvx_end

cvx_begin
    variable Z(n,n) hermitian toeplitz
    dual variable Q
    minimize( norm( Z - P, 'fro' ) )
    Z == hermitian_semidefinite( n ) : Q;
cvx_end

cvx_begin sdp
    variable Z(n,n) hermitian toeplitz
    dual variable Q
    minimize( norm( Z - P, 'fro' ) )
    Z >= 0 : Q;
cvx_end