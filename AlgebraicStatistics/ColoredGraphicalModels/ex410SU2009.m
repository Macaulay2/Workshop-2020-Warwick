% Example 4.10 from Sturmfels-Uhler (2009)
% "Multivariate Gaussians, Semidefinite Matrix Completion, 
% and Convex Algebraic Geometry"

% problem size
n = 5;

% create and solve the problem
cvx_begin sdp
  % A is a PSD symmetric matrix (n-by-n)
  variable K(n,n) symmetric;
  K >= 0;

  % constrained matrix entries.
  K(1,3) == 0;
  K(2,4) == 0;
  trace(K)==1;  
  
  % find the solution to the problem
  maximize( log_det( A ) )
  % maximize( det_rootn( A ) )
cvx_end

% display solution
disp(['Matrix A with maximum determinant (' num2str(det(A)) ') is:'])
A
disp(['Its eigenvalues are:'])
eigs = eig(A)