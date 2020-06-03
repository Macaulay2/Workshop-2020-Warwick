
-- Bergman fan work


needsPackage "Matroids"
M = matroid({0,1,2,3},{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}) --you can modify the matroid to play with Bergman
isWellDefined M
rank(M)
viewHelp "Matroids"

P = latticeOfFlats M
C = chains(P)
chains(P,rank(M)-1)
#(M)


-- First version of code for matroids of rank 2, needs some more looping for maximal chains
-- of length bigger than 1

-- BergmanList produces a list of maximal length proper chains
BergmanList =  (M) -> (
      P := latticeOfFlats(M);
      C := chains(P,rank(M)-1);
      L := {};
      n := #M-1;
      N := {toList(0..n)}; 
      for l in C do(
    	  if l != {{}} then (
	      if l != N then (
	    	  L = append(L,l);
	    	  );
	      );
    	  );
      L
      )




-- BergmanFan returns a matrix whose columns are the sums of elementary vectors associated to the elements in the maximal chain
BergmanFan =  (M) -> (
      P := latticeOfFlats(M);
      C := chains(P,rank(M)-1); -- this might have problems
      L := {}; 
      n := #M-1;
      N := {toList(0..n)}; 
      for l in C do(
    	  if l != {{}} then (
	      if l != N then (
	    	  L = append(L,l);
	    	  );
	      );
    	  );
      m:= mutableMatrix(ZZ,#M,#L);
      i:=0 ;
      for l in L do (
	  A := l#0;
	  for a in A do (
	      m_(a,i) = 1;
	      );
	  i = i+1;
	  );
      m
      )

BergmanFan(M)
BergmanList(M)

toList E


-- BergmaneI gives the indicator vector eI (as a list) for a subset I
-- assuming that ground set is [[n]] = {0,1,...,n} and I is a subset of the ground set
-- we should include the validity test for matroid in the beginning
BergmaneI = (M,I) -> (
    L := {};
    E := toList  M.groundSet;
    n := #E;
    for i in E do(
        if member(i,I) then  L= append(L, 1) else L= append(L,0);
        );
    L)


BergmaneI(M,{0,1})


-- BergmanchaineC gives the generators of a cone corresponding to a proper chain of flats C of matroid M
-- here by a proper chain we mean a chain that does not consist of empty set or the ground set.
-- this calls and therefore depends upon BergmaneI
-- this does not check whether C has flats in it or not
BergmanchaineC = (M,C) -> (
    L := {};
    E := toList M.groundSet; -- later we will remove redundancies
    n = #E;
    for F in C do(
	L = append(L,BergmaneI(M, F));
	);
    MatrixExpression L)
    
    
    
BergmanchaineC(M, {{0},{0,1,2,3}})
       
	

