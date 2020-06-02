newPackage(
        "MyExcellentPackage",
        Version => "0.1", 
        Date => "2 June 2020",
        Authors => {{Name => "Mike Stillman", 
                  Email => "", 
                  HomePage => ""}},
        Headline => "a function for displaying the Hilbert function",
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

restart
uninstallPackage "MyExcellentPackage"
restart
installPackage "MyExcellentPackage"
check MyExcellentPackage
viewHelp "MyExcellentPackage"
help MyExcellentPackage
help HF
restart
needsPackage "MyExcellentPackage"

-- eg:
path = prepend(path, "~/mike-M2/intro-to-packages")

R = ZZ/101[a..d]
I = monomialCurveIdeal(R, {1,3,4})
methods HF
help HF
HF(0,10,I)


TEST ///
-- test code and assertions here
-- may have as many TEST sections as needed
///

-- f11 works after you have done setupEmacs() inside M2.
-- BUT: for mac, and also perhaps ubuntu, system grabs f11.
--  system preferences.
