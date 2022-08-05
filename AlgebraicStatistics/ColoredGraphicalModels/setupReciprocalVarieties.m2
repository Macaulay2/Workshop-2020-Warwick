newPackage(
        "CatalecticantMatrices",
        Version => "1.0", 
        Date => "6 December 2020",
        Authors => {
            {Name => "Laura Brustenga i Moncusi", 
                Email => "brust@math.ku.dk", 
             Name => "Elisa Cazzador", 
                Email => "elisacaz@math.uio.no", 
             Name => "Roser Homs", 
                Email => "roser.homs@tum.de" 		
		}
            },
        Headline => "tools for catalecticant matrices"
        )

export {    
    "adjugate",
    "genericCatalecticantMatrix",
    "parametrizedImage",
    "symmetricEntries",
    "veronese",
    "toCat",
    "orthogonal"
    }

----------------------------------------------------------------------------------------
-- The present Macaulay2 file is an auxiliar package which is part of the             --
-- complementary material to the paper "Inverting catalecticants of ternary quartics" --                       
-- by Laura Brustenga i Moncusí, Elisa Cazzador and Roser Homs                        --
----------------------------------------------------------------------------------------


----------------------------------------------------------------------------------------
-- ADJUGATE MATRIX                                                                    --
----------------------------------------------------------------------------------------   
-- Compute adjugate of a symmetric matrix                                             --
-- Input:                                                                             --
--   * M, a matrix                                                                    --
-- Output:                                                                            --
--   * adjugate matrix of M                                                           --
----------------------------------------------------------------------------------------   

  adjugate = M -> ( 
    m := numcols M; 
    if m != numrows M then error "Matrix is not square";
    if entries M != entries transpose M then error "Matrix is not symmetric";
    adjugateM := for i to m-1 list ( 
	  for j to m-1 list (-1)^(i+j)*det(submatrix'(M, {j}, {i})) 
	  ); 
    matrix adjugateM
  ) 
 

----------------------------------------------------------------------------------------
-- GENERIC CATALECTICANT MATRIX                                                       --
----------------------------------------------------------------------------------------
-- Generic catalecticant matrix of (n+1)-ary forms of degree 2k with                  --
-- coefficients in the ring X (i.e. the entries of the matrix are the variables in X) --
-- Input:                                                                             --
--    * k, an integer                                                                 --
--    * n, an integer                                                                 --
--    * X, a ring  with indexed variables                                             --
-- Output:                                                                            --
--    * generic catalecticant matrix with entries in the ring X                       --
----------------------------------------------------------------------------------------

  genericCatalecticantMatrix = (k,n,X) -> (
    --projective dimension of space of catalecticants
    N:=#(gens X)-1;
    -- local variables for auxiliary ring R
    a:=symbol a;
    x:=symbol x;
    -- ring with variables a (coefficients of the form) and x (variables of the form)
    R:=QQ[a_0..a_N,x_0..x_n];
    -- list of variables a
    aList:=drop(gens R,-(n+1));
    -- list of variables x
    xList:=drop(gens R,N+1);
    -- basis of monomials of degree d=2k in variables x
    dBasis:=flatten entries basis(2*k,R,Variables=>xList);
    -- expression of a (n+1)-ary form of degree 2k 
    -- with variables x and coefficients a
    form:=sum(for i to #aList-1 list aList_i*dBasis_i);
    -- basis of monomials of degree k in variables x
    kBasis:=flatten entries basis(k,R,Variables=>xList);
    -- contraction (apolar action) of the form 
    -- with respect to operators given by monomials of degree k
    CM:=contract(matrix{kBasis},form);
    -- retrieve coefficients a of the contraction: catalecticant matrix
    (M,C):=coefficients(CM,Variables=>xList,Monomials=>kBasis);
    -- define ring map to return catalectican matrix to X
    f:=map(X,R,join(gens X,apply(xList,i->0)));
    -- retrieve catalecticant matrix in X, remove unnecessary information 
    matrix entries f(C)
  )

  

----------------------------------------------------------------------------------------
-- SYMMETRIC ENTRIES                                                                  --
----------------------------------------------------------------------------------------   
-- Extracts the entries of the upper triangular part of a symmetric matrix            --
-- Input:                                                                             --
--   * A, a matrix                                                                    --
--   * m, the size of A                                                               --
-- Output:                                                                            --
--   * list of entries of A in the upper triangular part                              --
----------------------------------------------------------------------------------------   
  
  symmetricEntries = (A, m) -> (
      flatten(
	for i to m-1 list(
	    for j from i to m-1 list A_(i,j)
	    )
      )
  )



----------------------------------------------------------------------------------------
-- PARAMETRIZED IMAGE                                                                 --
----------------------------------------------------------------------------------------
-- Parametrization for the image of a symmetric matrix via the restriction of the     --
-- inversion map to an LSSM                                                           --
-- Input:                                                                             --
--    * p, a matrix in the LSSM with QQ-entries                                       --
--    * L, the generic element in the LSSM                                            --
--    * S, the generic symmetric matrix of the same size of L                         --
-- Output:                                                                            --
--    * ideal of parametrization for the image of p via the inversion map             --
----------------------------------------------------------------------------------------

  parametrizedImage = (p, L, S) -> (
    -- ring with variables of the LSSM, the symmetric space and an auxiliary t
    t:=symbol t;
    R:=QQ[t]**(ring L)**(ring S);
    p':=sub(p,R);
    L':=sub(L,R);
    S':=sub(S,R);
    -- variables of symmetric space
    sVar:=support S';
    -- number of variables in the ring of S of symmetric matrices
    M:=#(sVar)-1;
    -- syze of matrices
    m:=numgens source L;
    -- rank of the point whose image we parametrize
    r:=rank p;
    -- line of full-rank points approaching p
    use R;
    line:=p'+t*L';
    -- phi(A + tL)
    imageLine:=adjugate line;
    -- divide all the entries by the maximum power ov t dividing them
    divImageLine:=matrix(
      for i to m-1 list (for j to m-1 list imageLine_(i,j) // 
        t^(m-r-1))
    );
    -- compute the limit setting t=0 
    limit:=sub(divImageLine, t=>0); 
    -- entries of the limit matrix
    entr:=symmetricEntries(limit,m);
    -- parametrization for phi(p)
    param:=ideal(
      for i to M list sVar_i-entr_i
    );
    -- return the parametrization in the ring of the product
    sub(param, ring(L)**ring(S))
  )




----------------------------------------------------------------------------------------
-- VERONESE                                                                           --
----------------------------------------------------------------------------------------
-- Coordinates for the d-uple embedding of a point                                    --
-- Input:                                                                             --
--    * P, a list of rational numbers, corresponding to the coordinates of a point    --
--    * d, the degree of the Veronese embedding                                       --
-- Output:                                                                            --
--    * a list of rational numbers, the coordinates of the d-uple embedding of P      --
----------------------------------------------------------------------------------------

  veronese = (P, d) -> (
    -- projective dimension of the space where P belongs to
    l:=length(P)-1;
    -- the ring of the projective space of P
    u:=symbol u;
    U:=QQ[u_0..u_l];
    -- list of degree-4 monomials in the variables of U
    dBasis:=first entries basis(d,U);
    -- substitute the coordinates of P in the monomials
    for monom in dBasis list(
      sub(monom, for i to l list u_i=>P_i)
      )
  )



----------------------------------------------------------------------------------------
-- TO CAT                                                                             --
----------------------------------------------------------------------------------------
-- Catalecticant matrix associated with the d-uple embedding of a point               --
-- Input:                                                                             --
--    * P, a list of rational numbers, corresponding to the coordinates of a point    --
--    * d, the degree of the Veronese embedding                                       --
-- Output:                                                                            --
--    * the catalecticant matrix associated with the d-uple embedding of P            --
----------------------------------------------------------------------------------------

  toCat = (P,d) -> (
    -- half of degree d, to build the generic catalecticant
    k:=d//2;
    -- projective dimension of the space where P belongs to
    n:=length(P)-1;
    -- ring of catalecticant matrices Cat(k,n)
    x:=symbol x;
    N:=binomial(n+2*k,2*k)-1;
    X:=QQ[x_0..x_N];
    -- generic catalecticant matrix
    cat:=genericCatalecticantMatrix(k,n,X);
    -- v_d(P), the d-uple embedding of P
    veron:=veronese(P,d); 
    -- substitute the coordinates of v_d(P) in the entries of catalecticant matrix
    sub(cat, for i to N list x_i=>veron_i)
    )


----------------------------------------------------------------------------------------
-- ORTHOGONAL                                                                         --
----------------------------------------------------------------------------------------
-- Orthogonal space to an LSSM                                                        --
-- Input:                                                                             --
--    * L, the generic matrix in the LSSM                                             --
--    * S, the generic symmetric matrix of the same size                              --
-- Output:                                                                            --
--    * the ideal of the orthogonal space to L in the ring of S                       --
----------------------------------------------------------------------------------------

  orthogonal = (L,S) -> (
    -- variables of L
    varsL:=support L;
    -- projective dimension of the space of L
    N:=#(varsL)-1;
    -- variables of S
    varsS:=support S;
    -- projective dimension of the sapce of S
    M:=#(varsS)-1;
    -- ring of the product, S-variables seen as coefficients
    R:=QQ[varsS][varsL];
    -- map  L and S in the ring R
    L':=sub(L,R);
    S':=sub(S,R);
    -- trace of the product between M and S
    traceProduct:=trace(L'*S');
    -- ideal of the orthogonal space to L
    use R;
    ortho:=trim ideal(for i to N list(
      coefficient(varsL_i, traceProduct)
      ));
    -- map the orthogonal in the ring of S
    sub (ortho, ring S)
    )
