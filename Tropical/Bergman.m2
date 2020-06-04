
-- Bergman fan work

needsPackage "Polyhedra"
needsPackage "Matroids"
M = matroid({0,1,2,3},{{0,1},{2,3}}) --you can modify the matroid to play with Bergman
rank(M)
viewHelp "Matroids"
isWellDefined M

M = matroid{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}
rank M
isWellDefined M

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



M = symbol M
-- BergmaneI returns the indicator vector (as a list) for a subset I of the ground set of a matroid M
-- We are assumin that the matroid is over the ground set [[n]] = {0,1,2,..,n}
-- does not include the test whether M is a valid matroid
BergmaneI = (M, I) -> (
    E := toList M.groundSet;
    n := #E;
    L = {};
    for i in E do(
	if member(i,I) then L =  append(L,1) else L = append(L,0);
	);   
    L
    )

BergmaneI(M,{1,2})

-- BergmanconeC returns the matrix of of generators of the cones corresponding to the chain of flats C
-- it does not check whether C is a chain of flat or not
-- This calles and therefore depends on BergmaneI
-- We will remove redundancies later
BergmanconeC  = (M, C) -> (
    E := toList M.groundSet;
    n := #E ;
    L := {};
    for F in C do(
	L = append(L, BergmaneI(M, F));
	); 
   transpose  matrix L
    )


-- BergmanFanI returns the fan
-- still needs to put some checks
-- this depends on 
BergmanFanI = (M) -> (
    E := toList M.groundSet;
    n := #E ;
    r := rank M;
     L := {};
    LM := latticeOfFlats M;
    redLM := dropElements(LM, {{}, E});
    rset := {toList(1..r-1)};
    redOrdcplx := maximalChains redLM;
    for C in redOrdcplx do(
	L = append(L, coneFromVData BergmanconeC(M,C));
	);
    fan L 
    )


U24 = uniformMatroid(2,4)
F = BergmanFanI U24
rays F
maxCones F



MK4 = matroid completeGraph 5
F = BergmanFanI MK4 -- this was taking some time on my machine
rays F
maxCones F
