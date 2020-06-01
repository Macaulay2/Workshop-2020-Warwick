# History 

1. Around 1997, Mike Stillman and Harry Tsai created some code to compute SAGBI bases. Two versions of the SAGBI code were written, one just using the front-end of Macaulay2 and the other built into the C/C++ engine of Macaulay2. 
2. Around 2006, Mike Sillman created a package called SubalgebraBases to provide documentation and the front-end package for computing SAGBI bases. Sometime before 2018, a change to the way monomial orders are represented broke at least one of the two implementations of SAGBI bases. 
3. In the Summer of 2019 (at the David Cox retirement conference), Michael Burr, Tim Duff, and Mike Stillman fixed the error in the engine due to the change in monomial orders for computing SAGBI bases. 
4. In the Summer of 2020 (at the Cleveland M2 workshop), Michael Burr, Tim Duff, and Elise Walker worked on modernizing and cleaning the code. Some future projects were also planned. 

# Where's the Code?

1. [Version on the master branch](https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/SubalgebraBases.m2). All tests pass since the release of Macaulay2 1.15.
2. However, the code was still a bit opaque. A stripped-down, working rewrite is available on this repository. Try installing it!

# What are (Khovanskii bases / SAGBI bases / Canonical Subalgebra Bases)?

1. See [here](./Subalgebra-Basics/StuCh11.pdf) for a first answer to this question.

# Potential Projects 

A work in progress. If you have other ideas, include them in your response to the group survey.

## New to M2

1. Write and export additional functions (such as NOBodies and NOBody volume) to compute interesting aspects of SAGBI/Khovanskii bodies. This should tie into other packages. 

2. Add / clean up the documentation, add (additional) tests, and prepare the package for future publication. 

3. C

## Some experience with M2

1. At the Cleveland M2 workshop, we rewrote and commented the old code. There are two functions which we didn't comment. Figure out what these functions are trying to do and make sure that they are sane. Then, rewrite them to be simpler. Some experience with M2: At the Cleveland M2 workshop, we rewrote and commented the old code for the engine-based computations. Most of the top level code repeats what was in the engine-based computations (but there are some differences). Attempt to merge the two functions. 

## Expert at M2

1. Create a SubAlgebra type. Originally, we had created a SubRing type as a type of hashtable, but it really should be a type of ring. It may also be better to focus on SubAlgebras instead of SubRings (but both types could be written concurrently). SubAlgebra would represent a finitely generated subalgebra. It would store its base ring R, its ambient ring S, and its generators G. The subalgebra can be represented as a map from R^|G| -> S. Elements could be given by points in R^G. Equality can be checked in S. A SAGBI basis could be used for membership (although there is a pure Groebner/elimination approach too). A good plan should be created here as this could be useful to many other projects. 

2. Create a valuation type. There are many design choices that should be considered here. Should there be a single valuation type or should all common valuation constructions have their own types? A single valuation type would not pollute the namespace, but just becomes a container to run arbitrary code (for arbitrary valuations). Some code has been started in this direction, but much more work (and a design) remains. This should also tie in with the TropicalGeometry package.

3. Implement additional algorithms (eg. intrinsic Groebner, 

