newPackage(
        "StatGraphs",
        Version => "0.1", 
        Date => "3 August 2020",
        Authors => {{Name=> "Carlos Amendola", 
	   Email=> "carlos.amendola@tum.de",
	   HomePage=>"http://www.carlos-amendola.com/"},
       
	  {Name => "Luis David Garcia Puente",
	   Email => "lgarcia@shsu.edu",
	   HomePage => "http://www.shsu.edu/~ldg005"},
       
          {Name=> "Roser Homs Pons", 
	   Email=> "roserhp@gmail.com",
	   HomePage=>"https://personal-homepages.mis.mpg.de/homspons/index.html"},
       
          {Name=> "Olga Kuznetsova", 
	   Email=> "kuznetsova.olga@gmail.com",
	   HomePage=>"https://okuznetsova.com"},
       
          {Name=> "Harshit J Motwani", 
	   Email=> "harshitmotwani2015@gmail.com"}},
        Headline => "Graphs specific for algebraic statistics",
        DebuggingMode => true
        )

export {
    "HF"
    }

HF = method()
HF(ZZ, ZZ, Ideal) := List => (lodegree, hidegree, I) -> (
    for i from lodegree to hidegree list hilbertFunction(i, I)
    )
HF(ZZ, ZZ, Module) := List => (lodegree, hidegree, I) -> (
    for i from lodegree to hidegree list hilbertFunction(i, I)
    )

beginDocumentation()

doc ///
  Key
    MyExcellentPackage
  Headline
    shortcut function for the Hilbert function
  Description
    Text
      This package contains a function to display
      the Hilbert function of an ideal in several degrees at once.
      
      For example, we compute the first few elements of the Hilbert function 
      of the rational quartic in $P^3$.
    Example
      R = ZZ/101[a..d];
      I = monomialCurveIdeal(R, {1,3,4})
      HF(0,5,I)
      assert({1, 4, 9, 13, 17, 21} == HF(0,5,I))
  SeeAlso
    HF
    hilbertFunction
    monomialCurveIdeal
///

doc ///
  Key
    HF
    (HF,ZZ,ZZ,Ideal)
  Headline
    shortcut for the Hiilbert function in a range of degrees
  Usage
    L = HF(lo, hi, I)
  Inputs
    lo:ZZ
      the low degree of the degree range
    hi:ZZ
      the top degree of the degree range
    I:Ideal
      a homogeneous ideal in a polynomial ring $S$
  Outputs
    L:List
      the values of the hilbert function $dim(S/I)_i$ for $i = lo, \ldots, hi$
  Description
    Text
      This is a quick way to get the several values of the Hilbert function 
    Example
      R = ZZ/101[a..d];
      I = monomialCurveIdeal(R, {1,3,4})
      HF(0,5,I)
      assert({1, 4, 9, 13, 17, 21} == HF(0,5,I))
  SeeAlso
    hilbertFunction
///

TEST ///
  R = ZZ/101[a..d]
  I = monomialCurveIdeal(R, {1,3,4})
  HF(0,5,I)
  assert({1, 4, 9, 13, 17, 21} == HF(0,5,I))
///

end--

