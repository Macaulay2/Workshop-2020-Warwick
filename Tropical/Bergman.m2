
-- Bergman fan work


needsPackage "Matroids"
M = matroid({0,1,2,3},{{0,1},{2,3}}) --you can modify the matroid to play with Bergman
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